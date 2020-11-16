#!/usr/bin/env python
# coding: utf-8
"""
    Script to run Findr on subsets (subsamples) of gene expression data.

PARAMETERS:

    - Number of segregants (samples) to use for this subset:
        size_of_subsets = 10
        
    - Number of subsets to draw for a given size:
        n_subsets = 4

References:

    [1] Ludl, A-A and Michoel, T (2020) Comparison between instrumental variable and mediation-based methods for reconstructing causal gene networks in yeast
    (submitted)
    arxiv: https://arxiv.org/abs/2010.07417


    [2] Lingfei Wang and Tom Michoel.
        Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data.
         PLOS Computational Biology, 13(8):1–26, 08 2017.

    [3] Lin S. Chen, Dipen P. Sangurdekar, and John D. Storey.
        trigger: Transcriptional Regulatory Inference from Genetics of Gene ExpRession,
         2007. R package version 1.16.0.

"""

vers='09'
#
import os
import numpy as np
import pandas as pan
import findr
import roman
import random

#
#   random RNG state/seed
#
random.seed()

#
#   PARAMETERS AND FOLDERS:
#
# Number of segregants to use for this subset
size_of_subsets = 1000
# number of subsets to draw
n_subsets = 8

# Output folders:
out_folder_top = '../../data/output/'

out_path = out_folder_top +   'subsampling_findr'+ vers +'/'
if not os.path.exists(out_path):
    os.makedirs(out_path)

this_out_path   =   out_path    +   'size_' +str(size_of_subsets)+'/'
if not os.path.exists(this_out_path):
    os.makedirs(this_out_path)
print('# Output paths have been created.')

#
#   MAIN function:
#
def main():
    """ Function doc """

    print('# 0.  Initialise findr.')
    l=findr.lib(loglv=12)

    print('# 1.  Load yeast data.')
    n_eqtls=None
    n_genes=None
    n_segregants=None

    in_path  = out_folder_top

    print('# 1.a  load genes-eqtls table.')
    gene_eqtl_table = pan.read_csv( in_path + 'strongest_eqtls_r_3columns.csv',
                                  #header=0,
                                   nrows=n_genes,
                                    delimiter=',')
    print(gene_eqtl_table.columns)
    print(gene_eqtl_table.shape)
    n_eqtls = gene_eqtl_table.shape[0]

    print('# 1.a.1  sort genes-eqtls table.')
    eqtls_genes_sorted_by_chr = roman_chromosome_sorting(gene_eqtl_table['pmarker'], gene_eqtl_table['gene'])

    print('# 1.b  load genotypes-experiments table.')
    genotypes = pan.read_csv( in_path + "genotypes_binary_strongest_eqtl.csv",
                                    nrows = n_segregants,
                                    delimiter=',')

    genotypes = genotypes.drop(genotypes.columns[0], axis=1)
    print(genotypes.columns)
    print(genotypes.index)
    print(genotypes.shape)

    print('# 1.b.1  sort columns of genotypes-experiments table.')
    eqtl_cols = eqtls_genes_sorted_by_chr['name']
    genotypes = genotypes[eqtl_cols]

    print('# check number of eqtls in genotypes')
    count_occuring_eq = sum( [ 1 for col in eqtl_cols if col in genotypes.columns] )
    print('  number of missing eqtls:',len(eqtl_cols) -count_occuring_eq)
    print(genotypes.shape)

    #
    #   Randomly draw n_subsets integers from the range of segregants
    #
    n_segregants_loaded = genotypes.shape[0]
    list_of_selected_segregants = random.sample( range(n_segregants_loaded), size_of_subsets )    
    list_of_selected_segregants = sorted(list_of_selected_segregants)
    print(list_of_selected_segregants)

    print('# 1.c  load gene expression table.')
    yeast_exp = pan.read_table(
                        in_path + "expression_statsmodels_linreg_residuals_01.csv",
                            nrows = n_segregants,
                            delimiter=',')
    yeast_exp = yeast_exp.drop(yeast_exp.columns[0:1], axis=1)
    print(yeast_exp.columns)
    print(yeast_exp.index)
    print(yeast_exp.shape)

    #   check if certain columns exist:
    sc_1=sorted([col for col in yeast_exp.columns if 'YAL' in col])
    print(sc_1)
    
    #   check how many gene names are missing in expr data
    print('# check number of genes in expression data')
    gene_cols = eqtls_genes_sorted_by_chr['gene']
    count_occuring_genes = sum( [ 1 for col in gene_cols if col in yeast_exp.columns] )
    print('  number of missing genes',len(gene_cols) -count_occuring_genes)
    missing_genes = [ col for col in gene_cols if col not in yeast_exp.columns]
    print('  set of missing genes', missing_genes)
    extra_genes = [ col for col in yeast_exp.columns if col not in list(gene_cols) ]

    print('  number of extra genes', len(extra_genes))
    print('  number of genes with eqtls', len(gene_cols))
    print('  sum of both', len(gene_cols) +len(extra_genes))

    print('# 1.c.1  sort columns of expression values.')
    gene_cols = eqtls_genes_sorted_by_chr['gene']
    yeast_exp_with_eqtls = yeast_exp[gene_cols]
    print('#   shape of yeast_exp_with_eqtls:',len(yeast_exp_with_eqtls.columns))


    yeast_exp_without_eqtls = yeast_exp[  extra_genes  ]
    print('#   shape of yeast_exp_without_eqtls:',len(yeast_exp_without_eqtls.columns))

    yeast_exp = pan.concat( [yeast_exp_with_eqtls,yeast_exp_without_eqtls] , sort=False, axis=1 )
    print('#   shape of yeast_exp:',len(yeast_exp.columns))
    expr_ordering_check = ( list(yeast_exp.columns) == list(yeast_exp_with_eqtls.columns) + list(yeast_exp_without_eqtls.columns) )
    print('#   check ordering of columns:', expr_ordering_check )

    print('### columns sanity check')
    print(yeast_exp.head())
    print(yeast_exp_with_eqtls.head())
    print(yeast_exp_without_eqtls.head())

    print('# 1.e  Convert matrices to findr format.')
    # Conversion
    #   the below converts from range [-1,0,1] to range [0,1,2]
    array_genotypes=( 0 < genotypes.to_numpy() ).astype('B')
    array_yeast_exp=yeast_exp.to_numpy().astype('float32')

    # transposes:
    array_genotypes = np.transpose(array_genotypes)
    array_yeast_exp = np.transpose(array_yeast_exp)

    print(array_genotypes.shape)
    print(array_yeast_exp.shape)

    #
    #  Loop over subsamples:
    #
    for i in range(n_subsets):
        i_outpath = this_out_path+str(i)
        if not os.path.exists(i_outpath):
            os.makedirs(i_outpath)
        #
        #   Randomly draw n_subsets integers from the range of segregants
        #
        n_segregants_loaded = genotypes.shape[0]
        list_of_selected_segregants = random.sample( range(n_segregants_loaded), size_of_subsets )    
        list_of_selected_segregants = sorted(list_of_selected_segregants)
        print(list_of_selected_segregants)

        #
        #   Save list_of_selected_segregants:
        #
        np.savetxt( i_outpath + '/selected_segregants.txt' ,list_of_selected_segregants, fmt='%d')

        a_subsample_of_genotypes = array_genotypes[:,list_of_selected_segregants]
        a_subsample_of_expr = array_yeast_exp[:,list_of_selected_segregants]

        print(a_subsample_of_genotypes.shape)
        print(a_subsample_of_expr.shape)

        print('# 2.  Run findr:')
        print('# 2.1  Run pij rank:')
        #  from the documentation: pijr=l.pij_rank(dt,dt2,nodiag=False)
        pijr=l.pij_rank(a_subsample_of_expr, a_subsample_of_expr, nodiag=True)
        print('#   saving file ...')
        np.savetxt( i_outpath + '/findr_pijrank' +'_yeast_'+vers ,pijr['p'])
        #np.savetxt( i_outpath + 'findr_pijrank' +'_yeast_'+vers+'.gz' ,pijr['p'])

        print('# 2.2.0  For pij gassist, define array of genes with eqtls:')
        #  define array of expression values corresponding to genotypes
        array_genes_with_eqtls = array_yeast_exp[:n_eqtls]
        print(array_genes_with_eqtls.shape)

        print('#  - b -  define Subsample of array of genes with eqtls:')
        a_subsample_of_genes_with_eqtls = array_genes_with_eqtls[:,list_of_selected_segregants]
        print(a_subsample_of_genes_with_eqtls.shape)

        print('# 2.2.1  Run pij gassist:')
        # pij_gassist
        #   This produces the probabilities for the novel causal inference test.
        #   Use  pij_gassist_trad for the traditional causal inference test.
        #   The first argument should be assorted to genotypes.
        #   This needs causal anchor information:
        #        eg eqtls where there should be one E for each A.
        #    Note: that appears to be the case in deed.
        #            select first columns and then check what can be done on the rest of the data.
        #
        print('genotypes shape', a_subsample_of_genotypes.shape)
        print('yeast exp shape', a_subsample_of_expr.shape)
        pijga=l.pij_gassist(a_subsample_of_genotypes, 
                             a_subsample_of_genes_with_eqtls,
                              a_subsample_of_expr, nodiag=True)
        print('#   saving file ...')
        np.savetxt( i_outpath + '/findr_pijga' +'_yeast_'+vers ,pijga['p'])
        #np.savetxt( i_outpath + '/findr_pijga' +'_yeast_'+vers+'.gz' ,pijga['p'])

        print('# 2.3  Run pijs gassist:')
        # pijs_gassist
        #   This produces ALL probabilities for causal inference tests as detailed in the article [ref 1].
        #   This needs causal anchor information: eg eqtls where there should be one E for each A.
        #      Note: that appears to be the case in deed.
        #            select first columns and then check what can be done on the rest of the data.
        #
        pijsga=l.pijs_gassist(a_subsample_of_genotypes, 
                             a_subsample_of_genes_with_eqtls,
                              a_subsample_of_expr, nodiag=True)
        pijsga.keys()
        count_ps=0
        print('# Saving pijs without diagonal to file:')
        for pi in list(pijsga.keys())[1:]:
            print('# saving '+pi)
            np.savetxt( i_outpath + '/findr_pijsga_'+pi+'_yeast_'+vers,pijsga[pi])
            #np.savetxt( i_outpath + '/findr_pijsga_'+pi+'_yeast_'+vers+'.gz',pijsga[pi])
            count_ps+=1
        print('# done saving '+str(count_ps)+' pijs to file.')
        
    print('# All done.')


def roman_chromosome_sorting(col_names, gene_names):
    """ Function doc:
        This function reorders the indices according to position along the chromosome.
        
        input:
            :col_names: (string) name of a marker (eg. SNP or eQTL).
            :gene_names: (string) name of a gene
            
        output:
            :dc: Data frame with name and gene ordered by position.
                Note: it is also saved to a file: sorted_column_names.csv
    """    
    df = pan.DataFrame()
    df['name'] = col_names
    df['gene'] = gene_names
    df['chr_label'] = df.apply (lambda row: row['name'].split(':')[0], axis=1)
    df['bases'] = df.apply (lambda row: row['name'].split('_')[1], axis=1)
    df['position'] = df.apply (lambda row: int( row['name'].replace(':','_').split('_')[1] ), axis=1)
    df['chromosome'] = df.apply (lambda row: roman.fromRoman( row['chr_label'].strip('chr') ), axis=1)
    dfs=df.sort_values(axis=0,by='chromosome position'.split() )
    print(dfs)
    dc = dfs['name gene'.split()]
    dc.to_csv('sorted_column_names.csv', header=1)
    return dc

#
#   RUN main
#
main()

# EOF

#!/usr/bin/env python
# coding: utf-8
#
# A python script to run Findr [1] on yeast data from [2].
#
#  References:
#
#   1. Lingfei Wang and Tom Michoel.
#        Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data.
#         PLOS Computational Biology, 13(8):1–26, 08 2017.
#          paper: https://doi.org/10.1371/journal.pcbi.1005703
#           repository: https://github.com/lingfeiwang/findr
#
#   2. Albert, F. W., Bloom, J. S., Siegel, J., Day, L., & Kruglyak, L. (2018).
#       Genetics of trans-regulatory variation in gene expression. Elife, 7, e35471. 
#        doi:10.7554/eLife.35471
#         link: https://doi.org/10.7554/elife.35471
#
vers='01'
#
import os
import numpy as np
import pandas as pan
import findr
import roman

def roman_chromosome_sorting(col_names, gene_names, out_path='../data/output/'):
    """ Function doc
        A function to sort the gene labels according to their chromosome numbering and position along the chromsome.
        Inputs: 
                col_names: list of strings,
                gene_names: list of strings,
                outpath: string.
        Output: dc: pandas DataFrame and written to csv
    """
    df = pan.DataFrame()
    df['name'] = col_names
    df['gene'] = gene_names
    #print(df)
    df['chr_label'] = df.apply (lambda row: row['name'].split(':')[0], axis=1)
    df['bases'] = df.apply (lambda row: row['name'].split('_')[1], axis=1)
    df['position'] = df.apply (lambda row: int( row['name'].replace(':','_').split('_')[1] ), axis=1)
    df['chromosome'] = df.apply (lambda row: roman.fromRoman( row['chr_label'].strip('chr') ), axis=1)
    #print(df)
    dfs=df.sort_values(axis=0,by='chromosome position'.split() )
    print(dfs)
    dc = dfs['name gene'.split()]
    dc.to_csv( out_path + 'sorted_column_names.csv', header=1)
    return dc


def main():
    """ Function doc
        Run findr P0: pij_rank, and causl tests: pij_gassist on yeast data and save results
    """
    print('# Run Analysis with Findr inference tests:')

    out_path = '../data/output/'
    in_path  = out_path
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    this_out_path  =  out_path  +  'findr_yeast/'
    if not os.path.exists(this_out_path):
        os.makedirs(this_out_path)
    print('# Output paths have been created.')

    print('# 0.  Initialise findr.')
    l=findr.lib(loglv=12)

    print('# 1.  Load yeast data.')
    n_eqtls=None
    n_genes=None
    n_segregants=None

    print('# 1.a  load genes-eqtls table.')
    gene_eqtl_table = pan.read_csv( out_path + 'strongest_eqtls_r_3columns.csv',
                                  #header=0,
                                   nrows=n_genes,
                                    delimiter=',')
    print(gene_eqtl_table.columns)
    print(gene_eqtl_table.shape)
    n_eqtls = gene_eqtl_table.shape[0]
    print('# 1.a.1  sort genes-eqtls table.')
    eqtls_genes_sorted_by_chr = roman_chromosome_sorting(gene_eqtl_table['pmarker'],
                                        gene_eqtl_table['gene'],
                                        out_path)


    print('# 1.b  load genotypes-experiments table.')
    genotypes = pan.read_csv( out_path + "genotypes_binary_strongest_eqtl.csv",
                                    nrows = n_segregants,
                                    delimiter=',')
    genotypes = genotypes.drop(genotypes.columns[0], axis=1)
    print(genotypes.columns)
    print(genotypes.shape)
    print('# 1.b.1  sort columns of genotypes-experiments table.')
    eqtl_cols = eqtls_genes_sorted_by_chr['name']
    genotypes = genotypes[eqtl_cols]

    print('# 1.b.2 check number of eqtls in genotypes')
    count_occuring_eq = sum( [ 1 for col in eqtl_cols if col in genotypes.columns] )
    print('#  number of missing eqtls:',len(eqtl_cols) -count_occuring_eq)
    print(genotypes.shape)

    print('# 1.c  load gene expression table.')
    yeast_exp = pan.read_table(
            out_path + "expression_statsmodels_linreg_residuals_01.csv",
                                    #usecols = range(1,n_eqtls+1),
                                    nrows = n_segregants,
                                    delimiter=',')

    yeast_exp = yeast_exp.drop(yeast_exp.columns[0:1], axis=1)
    print(yeast_exp.columns)
    print(yeast_exp.shape)

    #   check for YAL genes:
    #sc_1=sorted([col for col in yeast_exp.columns if 'YAL' in col])
    #print(sc_1)
    #   check how many gene names are missing in expr data
    print('# check number of genes in expression data')
    gene_cols = eqtls_genes_sorted_by_chr['gene']
    count_occuring_genes = sum( [ 1 for col in gene_cols if col in yeast_exp.columns] )
    print('#  number of missing genes',len(gene_cols) -count_occuring_genes)
    missing_genes = [ col for col in gene_cols if col not in yeast_exp.columns]
    print('#  set of missing genes', missing_genes)
    extra_genes = [ col for col in yeast_exp.columns if col not in list(gene_cols) ]

    print('#  number of extra genes', len(extra_genes))
    print('#  number of genes with eqtls', len(gene_cols))
    print('#  sum of both', len(gene_cols) +len(extra_genes))

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

    #print('### columns sanity check')
    #print(yeast_exp.head())
    #print(yeast_exp_with_eqtls.head())
    #print(yeast_exp_without_eqtls.head())
    
    #
    #   save columns in order
    #
    yeast_exp.columns.to_frame().reset_index(drop=True).to_csv( this_out_path + 'columns_in_order_for_findr.csv', header=1)


    print('# 1.d  Convert matrices to findr format.')
    # Conversion
    #   the below converts from range [-1,0,1] to range [0,1,2]
    array_genotypes=( 0 < genotypes.to_numpy() ).astype('B')
    array_yeast_exp=yeast_exp.to_numpy().astype('float32')

    # transposes
    array_genotypes = np.transpose(array_genotypes)
    array_yeast_exp = np.transpose(array_yeast_exp)

    print('# 2.  Run findr:')
    print('# 2.1  Run pij rank: P0')
    # pij_rank: 
    #   This is the coexpression (correlation based) test P0.
    #   Note that the resulting matrix is not strictly symmetric,
    #    because it is estimated from the observed distribution
    #     of LLR test statistics for each A gene separately.
    #
    pijr=l.pij_rank(array_yeast_exp, array_yeast_exp, nodiag=True)
    np.savetxt( this_out_path + 'findr_pijrank_yeast_01', pijr['p'])

    print('# 2.2.0  For pij gassist, define array:')
    #  define array of expression values corresponding to genotypes
    array_genes_with_eqtls = array_yeast_exp[:n_eqtls]
    print(array_genes_with_eqtls.shape)

    print('# 2.2.1  Run pij gassist: P (Findr test)')
    # pij_gassist:
    #   This produces the probabilities for the novel causal inference test.
    #   The first argument should be assorted to genotypes.
    #   This needs causal anchor information: eg eqtls where there should be one E for each A gene.
    #
    print('genotypes shape',array_genotypes.shape)
    print('yeast exp shape',array_yeast_exp.shape)

    pijga=l.pij_gassist(array_genotypes, 
                         array_genes_with_eqtls,
                          array_yeast_exp, nodiag=True)

    np.savetxt( this_out_path + '/findr_pijga_yeast_01',pijga['p'])

    print('# 2.3  Run pijs gassist: P2, P3, P4, P5')
    # pijs_gassist:
    #   This produces ALL probabilities for causal inference tests as detailed in the article [ref 1].
    #   This needs causal anchor information: eg eqtls where there should be one E for each A gene.
    #
    pijsga=l.pijs_gassist(array_genotypes,
                           array_genes_with_eqtls,
                            array_yeast_exp, nodiag=True)
    pijsga.keys()

    count_ps=0
    print('#  Saving pijs without diagonal to file:')
    for pi in list(pijsga.keys())[1:]:
        print('#   saving '+pi)
        np.savetxt( this_out_path + '/findr_pijsga_'+pi+'_yeast_01', pijsga[pi])
        count_ps+=1
    print('#  done saving '+str(count_ps)+' pijs to file.')

    print('# Done with Findr.')


# Run main:
main()

# EOF

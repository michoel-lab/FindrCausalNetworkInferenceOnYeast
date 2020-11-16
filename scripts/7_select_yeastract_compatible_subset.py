"""
#  
#  Copyright 2020 AL at UIB.
#
#   Select subset of yeastract matrices to compare with findr data.
#

    NOTE: This should be run after running Findr.

"""
import os
import pandas as pan

import index_labels

vers='01'

def main():
    """ Function doc 
    
        Method to select subset of yeastract matrices for comparison
         and further analysis of causal matrices inferred with findr.
        
        Load the following files:
        
        1.  Yeast gene library from Ensembl library:

            Files:   cleaned_genes_ensembl83_yeast_step1.txt
                       and
                     cleaned_genes_ensembl83_yeast_step2.txt
            Header:
            #!genome-build SGD R64-1-1
            #!genome-version R64-1-1
            #!genome-date 2011-09
            #!genome-build-accession GCA_000146045.2
            #!genebuild-last-updated 2011-12
        
        
        2.  Labels for causal network matrix inferred with findr  [Wang, Michoel (2017)]
            from data from Albert et al. elife ( 2018, 35471v1) paper.
        
            Files:  columns_in_order_for_findr_01.csv
                      and
                    findr_pijsga_p5_yeast_01

            Note: Rows are in the same order as columns,
            but only the eqtls are in rows so we need 
            the number of rows,
            eg. from the shape of the p matrix (eg. p5).
        
        3.  Yeastract regulatory network from:
            http://www.yeastract.com/formregmatrix.php
            
            Format is : 
                rows = TF,
                columns = genes.
            
            6 files (.csv.gz) in 2 groups of 3 :
            
                a) files with DNA binding only, DNA +/and Expr:
                
                    onlyDNABinding_RegulationMatrix_Documented_*
                    DNABinding_AND_Expr_TF_act_OR_inh_RegulationMatrix_Documented_*
                    DNABinding_PLUS_Expr_TF_act_OR_inh_RegulationMatrix_Documented_*
                
                b) files with 3 Expr only variants (act, inh, act OR inh).
                
                    Expr_TF_act_RegulationMatrix_Documented_*
                    Expr_TF_inh_RegulationMatrix_Documented_*
                    Expr_TF_act_OR_inh_RegulationMatrix_Documented_*
                        
        --- --- ---
        
        References:
        
        Genetics of trans-regulatory variation in gene expression
        FW Albert, JS Bloom, J Siegel, L Day, L Kruglyak
        Elife (2018)
            https://elifesciences.org/articles/35471v1.pdf
        
        Findr paper:
          Wang, Michoel (2017) PLoS Comput Biol 13(8): e1005703.
          https://doi.org/10.1371/journal.pcbi.1005703

        http://www.yeastract.com/
          YEASTRACT+: a portal for cross-species comparative genomics of transcription regulation in yeasts.
          Nucleic Acids Research, 48(D1):D642-D649
          (doi:10.1093/nar/gkz859) 
          P.T. Monteiro, J. Oliveira, P. Pais, M. Antunes, M. Palma, M. Cavalheiro, M. Galocha, C.P. Godinho, L.C. Martins, N. Bourbon, M.N. Mota, R.A. Ribeiro, R.Viana, I. Sá-Correia, M.C. Teixeira (2020)

        Ensembl library:
            ftp://ftp.ensembl.org/pub/release-83/gff3/saccharomyces_cerevisiae/
            http://www.ensembl.org/info/website/archives/index.html
            http://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index?db=core
            
    """
    path_input_data = '../data/input/'
    path_output_data = '../data/output/'
    
    out_vers_path = path_output_data + 'yeastract/'+ vers +'/'
    if not os.path.exists(out_vers_path):
        os.makedirs(out_vers_path)
    
    print('#  Selecting compatible subset of yeast regulatory network from yeastract with respect to elife paper.')
    print('#  1. Load ensembl yeast gene library:')
    #
    #   1.  Load ensembl list as dictionnaries.
    #       These files have been preprocessed with sed,
    #       They only include lines with a Name= entry.
    path_ensembl = path_output_data 

    #   The file used below only contains those entries with Name and ID.
    ensembl_file_full = 'cleaned_genes_ensembl83_yeast_step1.txt'
    df_full = pan.read_csv( path_ensembl + ensembl_file_full, skiprows=0,
                        names='chromosome position gene_str'.split(), 
                        usecols=[0,3,8],
                        comment='#',
                        sep='\t| ', engine='python')

    df_full['id'] = df_full['gene_str'].str.split(':').str[1]
    full_gene_list = df_full['id'].to_list()

    #   The file used below only contains those entries with Name and ID.
    ensembl_file_idn = 'cleaned_genes_ensembl83_yeast_step2.txt'
    df = pan.read_csv( path_ensembl + ensembl_file_idn, skiprows=0,
                        names='chromosome position gene_str name_str'.split(), 
                        usecols=[0,3,8,9],
                        comment='#',
                        sep='\t| ', engine='python')

    print('#  data frame as loaded')
    print(df)

    #   Create new column name and id
    df['name'] = df['name_str'].str.split('=').str[1]
    df['id'] = df['gene_str'].str.split(':').str[1]
    #   Sort it
    dfs=df.sort_values(axis=0,by='chromosome position'.split() )
    #   Make gene_name_dict
    print('#  data frame after column operations and sorting')
    print(dfs)
    df['name'] = df['name'].str.upper()
    df['id'] = df['id'].str.upper()
    dfi = dfs.set_index('name')
    gene_name_dict = dfi['id'].to_dict()
    print('#  number of genes:', len(gene_name_dict))

    print('#  2. Load elife gene list:')
    #
    #   2.  Load elife (Albert et al.) list of genes and eqtls
    #
    findr_data_path = path_output_data + 'findr_yeast/'
    file_columns_findr_order = 'columns_in_order_for_findr.csv'
    cols_findr_ordering = pan.read_csv( findr_data_path + file_columns_findr_order,
                       header=0)
    #
    #   Note: Rows are in the same order as columns,
    #         but only the eqtls are in rows 
    #         so we need the number of rows, eg. from the shape of the p matrix.
    #
    findr_index_labels = index_labels()
    #findr_index_labels.set_columns( cols_findr_ordering['name'] )
    findr_index_labels.set_columns( cols_findr_ordering['0'] )
    
    sub_folders = ''
    file_p5_findr = 'findr_pijsga_p5_yeast_01'
    p_matrix = pan.read_csv( findr_data_path + sub_folders + file_p5_findr,
                       header=None, sep=' ' )
    n_eqtls, n_genes = p_matrix.shape
    print('#    number of eqtls and genes',n_eqtls, n_genes)
    rows_findr_ordering = cols_findr_ordering[ :n_eqtls ]
    findr_index_labels.set_rows( rows_findr_ordering['0'] )

    print('# list of p5 rows')
    print(type(rows_findr_ordering['0']))

    print('#  3. Load yeastract matrices:')
    #
    # 3.  Load yeastract matrix
    #
    filepath_yt = path_input_data + 'yeastract/'
    yeastract_data = load_yeastract(filepath_yt, data_sel='all')
    
    print('#  4.  Select shared genes.')
    #
    #   4.  Select shared genes.
    #
    for ( yt_name, t_yeastract ) in yeastract_data:
        print('#  - Select common genes for ', yt_name)
        #
        #   a.  Map gene labels to elife version nomenclature using ensembl library.
        compatible_t_yeastract  =  relabel_matrix(  t_yeastract, gene_name_dict, full_gene_list )
        outfile = yt_name + '_yeastract_matrix_relabeled_for_findr.csv'
        compatible_t_yeastract.to_csv( out_vers_path + outfile, header=True, index=True, mode='w', sep=',' )
        #
        #   b.  Save output
        outfile = yt_name + '_yeastract_matrix_for_findr.csv'
        #   Save as csv (outfile, relabelled_matrix)
        yt_submat_common = select_common_rows_and_columns( compatible_t_yeastract, findr_index_labels )
        yt_submat_common.to_csv( out_vers_path + outfile, header=True, index=True, mode='w', sep=',' )
        print('#    common matrix shape',yt_submat_common.shape)
###

def  relabel_matrix(matrix, index_name_dict, full_gene_list):
    """ Function doc 
    
        Use dictionnary to relabel rows and columns of matrix.
        Intended to be used for instance with genes.
        
        Input:
            matrix:     (regulatory) matrix with row and column labels of genes or transcription factors (TF). 
            index_name_dict:    dictionnary to map TF names to gene names (obtained from ensembl library).
        
        Output:
            sub_matrix: a submatrix of the input matrix, that contains 
                        the rows and columns with labels
                        that are compatible with the index_name_dict.
    """
    old_cols = matrix.columns
    old_rows = matrix.index.values
    print('old_cols, n=',len(old_cols) )
    print(old_cols)
    print('old_rows, n=', len(old_rows))
    print(old_rows)
    print('#  parsing columns:')
    col_subset, new_cols_dict = select_compatible_labels(old_cols, index_name_dict, full_gene_list)
    print('#  parsing rows:')
    row_subset, new_rows_dict = select_compatible_labels(old_rows, index_name_dict, full_gene_list)
    sub_matrix = matrix.loc[ row_subset, col_subset ]
    sub_matrix.rename( columns = new_cols_dict, inplace=True)
    sub_matrix.rename( index = new_rows_dict, inplace=True)
    print('  number of compatible columns',len(col_subset))
    print('  number of compatible rows',len(row_subset))
    return sub_matrix


def  select_compatible_labels( old_labels, label_dict, full_label_list ):
    """ Function doc 
    
        Parse all old_labels to find those compatible with the label_dict.
        Incompatible labels are printed in output for logging.
        
        Input: 
            old_labels:     labels of row or column.
            label_dict:     dictionnary to map TF names to gene names (obtained from ensembl library).
    
        Output:
            label_subset:   subset of compatible labels.
            new_labels_dict:    dictionnary of compatible TF names to gene names.
    """
    label_subset=[]
    new_labels_dict={}
    
    for la in  old_labels:
        # The below line tested if the dash character was causing labels to be misclassified.
        li=la.replace(',','%2')
        
        label_x = li.upper()
        if li[-1]=='p':
            #   Note: this is for transcription factors (TF) from yeastract, which end in letter p.
            upper_li=li[:-1]
            label_x = upper_li.upper()
            
        if label_x in label_dict.keys() :
            label_subset.append( li )
            new_labels_dict[li] = label_dict[ label_x]
        
        elif label_x in full_label_list :
            label_subset.append( li )
            new_labels_dict[li] = label_x
            
        else:
            #   print output for log.
            print(' This label is incompatible:', li)            
    return label_subset, new_labels_dict
    

        

def select_common_rows_and_columns( a, selected_indices ):
    """ Function doc 
    
        Select rows and columns in an array. 
    
        Input:
                a:      pandas DataFrame,
                selected_indices:   index_labels object that contains the selected indices.
                
        Output:
                a_common:   submatrix with common rows and columns.
    """
    common_columns = list(set(a.columns).intersection(set(selected_indices.get_columns())))
    common_rows = list(set(a.index.values).intersection(set(selected_indices.get_rows())))
    print('  common columns', common_columns)
    print('  common rows', common_rows)
    a_common = a.loc[ common_rows, common_columns ]
    return a_common


def  load_yeastract(data_path,data_sel='all'):
    """ Function doc

    First load DNA binding only, Expr only, then DNA and Expr, DNA + Expr

    The optinal variable data_sel controls
        data_sel: string, possible values: all, DNA_only, Expr_only .

    Input:
            data_path: string,
            data_sel: string, states the selection to be chosen.
    Output:
            z: list containing tuples (name, array).
    """
    from glob import glob
    suffix = '.csv.gz'
    full_yt_file_list = glob(data_path+'*' + suffix)
    if data_sel=='all':
        yt_file_list = full_yt_file_list
        #pass
    elif data_sel=='DNA_only':
        yt_file_list = [ x for x in full_yt_file_list if 'DNA' in x  ]
    elif data_sel=='Expr_only':
        yt_file_list = [ x for x in full_yt_file_list if 'Expr' in x  ]
    if data_sel=='four':
        yt_file_list = [ x for x in full_yt_file_list if 'DNA' in x  ]
        yt_file_list += glob(data_path+ 'Expr_TF_act_OR_inh*' + suffix)
    yt_file_list = list(set(yt_file_list))
    print('#    number of yeastract files to load:', len(yt_file_list) )
    print('#    yeastract files to load:')
    print(yt_file_list)
    z=[]
    for fn in yt_file_list:
        print('  load  file : ' +fn)
        m = pan.read_csv(  fn,
                           header=0,
                           index_col=0,
                            delimiter=';')
                            #delimiter=',')
        m_name = fn.split('/')[-1].split('_yeastract')[0]
        z.append( (m_name,m) )
    return z

#
#   RUN main:
#
main()

# EOF

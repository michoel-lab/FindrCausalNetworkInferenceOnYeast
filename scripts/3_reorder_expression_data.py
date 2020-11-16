#!/usr/bin/env python
# coding: utf-8
#
#  Reorder the gene expression data matrices.
#   This puts genes in the same order as in the eQTL list.
#
# 0. Imports
#
import os
import pandas as pan
import numpy as np

in_path = '../data/input/'
out_path = '../data/output/'

if not os.path.exists(out_path):
    os.makedirs(out_path)

print('# Reorder the gene expression data matrices:')
print('#    Load data:')

#  1. Load expression data:
expr_values = pan.read_csv( in_path + "SI_Data_01_expressionValues.txt.zip",sep='\t')
expr_values.shape
expr_values.columns

#  2. Load ordering from eqtls ...
ser = pan.read_csv( out_path + 'strongest_eqtls_r.csv')
#  Sort them by 'alphabetical' order of eqtls
sorted_ser = ser.sort_values('pmarker')
#  Check statistics
#    a) for eqtls
eqtl_labels = sorted_ser.pmarker
#   below is the statistic of occurences of eqtls:
eqtl_labels.value_counts().value_counts()

#  Check statistics
#    b) for genges
gene_labels = sorted_ser.gene
#   below is the statistic of occurences of eqtls:
gene_labels.value_counts().value_counts()
#   create lists to map eqtl to gene names.
list_genes = []
list_eqtls = []
k=0
for line in sorted_ser.iterrows():
    list_genes.append( line[1].gene )
    list_eqtls.append( line[1].pmarker )
    k+=1

print('#    Reordering:')
#
#   3. Build data frame from expression data
#
sel_exp=pan.DataFrame()
#for col_name in eqtl_labels:
for i in range(len(list_genes)):
    sel_exp.insert( i, list_genes[i], expr_values[list_genes[i]] )
sel_exp.shape

# find and add the other expression data 
# and add a table that keeps the ordering information.
other_gene_labels=set(expr_values.columns)-set(list_genes)

numbers_of_genes = [len(other_gene_labels), len(expr_values.columns), len(list_genes)]
print('       numbers_of_genes: ', numbers_of_genes)
print('       genes with eQTLs: ', len(set(list_genes)))
print('       all genes: ', numbers_of_genes[0]+numbers_of_genes[2])

list_other_gene_labels=list(other_gene_labels)
n_eqtls = len(list_genes)
for i in range(len(list_other_gene_labels)):
    sel_exp.insert( n_eqtls+i, list_other_gene_labels[i],
                   expr_values[list_other_gene_labels[i]] )
sel_exp.shape

check_columns = (sel_exp.columns == list_genes + list_other_gene_labels).all()
print( '    Check if all expected columns are present:', check_columns)
print( '    number of columns selected:', len(sel_exp.columns) )

print('#    Writing output:')
#   4. Write output:
sel_exp.to_csv( out_path + 'expression_reordered_01.csv')
np.savetxt( out_path + 'list_genes_in_order_01.csv', sel_exp.columns.to_numpy(), fmt='%s')

#
#   save columns in order
#
sel_exp.columns.to_frame().reset_index(drop=True).to_csv( out_path + 'columns_in_order_for_findr_b.csv', header=1)

print( '# Done.' )
# EOF

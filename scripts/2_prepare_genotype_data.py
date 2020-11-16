#!/usr/bin/env python
# coding: utf-8
#
#  Preprocessing the genotype data matrices with genes
#   in the same order as in the eqtl list.
#
#   NOTE: need to prepare eQTL list first !
#
# 0. imports first
import pandas as pan
import numpy as np

print( '# Preprocessing genotype matrices for Findr:' )
#
#   0. output path
#
out_path = "../data/output/"

print( '#   Load genotype matrices:' )
#
#  1. Load data:
path_expr = "../data/input/"

genotypes = pan.read_csv( path_expr + "SI_Data_03_genotypes.txt.zip",sep='\t')
print( genotypes.shape )
print( genotypes.columns )
#
#   load ordering from eqtls ...
path_eqtls = out_path

ser = pan.read_csv( path_eqtls + 'strongest_eqtls_r.csv')
sorted_ser = ser.sort_values('pmarker')
eqtl_labels = sorted_ser.pmarker
#   below is the statistic of occurences of eqtls:
print( eqtl_labels.value_counts().value_counts() )
print( eqtl_labels )

selected_genotypes = genotypes.reindex( eqtl_labels.values )
print( selected_genotypes.columns.values )
print( selected_genotypes.columns.value_counts().isin([1]).all() )

print( sum(selected_genotypes.columns.value_counts()) )
print( sum(genotypes.columns.value_counts()) )


print( '#   Sort genotype matrices:' )
#
#   2. Sort the data according to label ordering (chromosome and genes).
#
sel_gt=pan.DataFrame()
i=0
for col_name in eqtl_labels:
    sel_gt.insert( i, col_name, genotypes[col_name],
                  allow_duplicates=True )
    i=i+1
print( sel_gt.shape )
print( sel_gt.columns.value_counts().value_counts() )
print( (sel_gt.columns.values == eqtl_labels).all() )

print( '#   Write output:' )
#   transform format of genotypes to binary:
sel_gt = (sel_gt > 0).astype(int)
#
#   3. Write output
#       save genotypes to file
sel_gt.to_csv( out_path + 'genotypes_binary_strongest_eqtl.csv' )
#print(sel_gt)

print( '# Done.' )
# EOF

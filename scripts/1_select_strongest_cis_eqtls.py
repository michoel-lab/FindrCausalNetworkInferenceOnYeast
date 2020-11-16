#!/usr/bin/env python
# coding: utf-8
#
# Python script to load and process yeast eQTLs from Albert et al. (2018):
#
#   Albert, F. W., Bloom, J. S., Siegel, J., Day, L., & Kruglyak, L. (2018).
#       Genetics of trans-regulatory variation in gene expression. Elife, 7, e35471. doi:10.7554/eLife.35471
#
#   paper:  https://doi.org/10.7554/elife.35471
#   data:   https://figshare.com/s/83bddc1ddf3f97108ad4
#
# 0. imports 
#
import pandas as pan
import numpy as np

out_path = '../data/output/'

print('# Load data for yeast eQTLs from Albert et al. (2018):')

#  1. Load data:
df_eqtl_info = pan.read_excel("../data/input/SI_Data_04_eQTL.xlsx")

df_eqtl_info.columns
df_eqtl_info.shape
df_eqtl_info.describe()

#  2. Select cis eQTLs
cis_indices = df_eqtl_info['cis']==1
df_cis=df_eqtl_info[cis_indices]

df_cis.describe()
df_cis.shape

#  3.  Sort cis eQTLs by r value (correlation coeff)
df_cis_sorted = df_cis.sort_values('r')

#   sort cis eQTLs by r absolute value (correlation coeff)
r_sort_indices  = df_cis['r'].abs().sort_values().index
df_cis_sorted_r = df_cis.reindex( r_sort_indices )
df_cis_sorted_r.to_csv(out_path+'cis_eqtls_sorted_by_r.csv')
#   write output to csv file:
df_cis_sorted_r.to_csv(out_path+'cis_eqtls_sorted_by_r_threecolumns.csv', columns=['gene','pmarker','r'])

r_sort_indices  = df_cis['r'].abs().sort_values().index
df_cis_sorted_r = df_cis.reindex( r_sort_indices )

#
#   check counts of values
df_cis_sorted_r['gene'].value_counts().value_counts()
#
#  results should be:
#   1    2800
#   2      83
#   3       1
#
#  Therefore there should be 2884 gene entries with one strongest eqtl.
#
#   histogram of gene repetitions in eqtl file
df_eqtl_info['gene'].value_counts().value_counts()

print('#    Select strongest eqtls by r coeff')
#
#  4. Select strongest eqtls by r coeff
#
df=df_cis_sorted_r
df['r_abs'] = df['r'].abs()
grouped = df.groupby('gene')
max_indices = grouped['r_abs'].idxmax()
#   Note: These max_indices are the indices in the original data, that is df_eqtl_info.
sel_max = df_eqtl_info.reindex( max_indices )
#   write output csv files:
sel_max.to_csv(out_path+'strongest_eqtls_r.csv')
sel_max.to_csv(out_path+'strongest_eqtls_r_3columns.csv',columns=['gene','pmarker','r'])
sel_max.to_csv(out_path+'strongest_eqtls_r_list.csv',columns=['gene','pmarker'])

# below is the statistic of occurences of eqtls:
sel_max['pmarker'].value_counts().value_counts()

sel_max.describe()
sel_max['gene'].value_counts().value_counts()

sel_max.columns

#   check absolute values of r:
da = sel_max
da['r_abs'] = da['r'].abs()
da.describe()
da['r_abs'].hist()

ax=da['r'].hist(bins=64)
ax.set_xlabel("r")
ax.set_ylabel("Counts")

print('# Done.')

# EOF

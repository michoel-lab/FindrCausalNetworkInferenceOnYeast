#!/bin/sh

echo 'Run Full Pipeline of Findr Analysis on Yeast Data:'

echo ' Beware this could take some time !'

log_file=findr_log_v1

  echo '# 1. Select strongest cis-eQTLs' | tee -a $log_file
  python  1_select_strongest_cis_eqtls.py  | tee -a $log_file

  echo '# 2. Prepare genotypes matrix' | tee -a $log_file
  python  2_prepare_genotype_data.py  | tee -a $log_file

  echo '# 3. Prepare expression data' | tee -a $log_file
  python 3_reorder_expression_data.py | tee -a $log_file

  echo '# 4. Prepare genotypes matrix' | tee -a $log_file
  python 4_covariate_regression_on_expression_data.py | tee -a $log_file

  echo '# 5. Preprocessing of Ensembl annotations for yeast genome' | tee -a $log_file
  sh 5_preprocessing_ensembl_data.sh    | tee -a $log_file

  echo '# 6. Prepare genotypes matrix' | tee -a $log_file
  python 6_run_findr.py  | tee -a $log_file

  echo '# 7. Select compatible submatrix from yeastract data' | tee -a $log_file
  python  ./7_select_yeastract_compatible_subset.py  | tee -a $log_file



echo 'Done'

# FindrCausalInferenceOnYeast [![DOI](https://zenodo.org/badge/313364218.svg)](https://zenodo.org/badge/latestdoi/313364218)

This repository contains python code to reproduce the causal networks inferred from differential gene expression in yeast as published in [1].
The code relies on Findr which should be installed according to the
instructions provided in the corresponding repository [2].
Releases are archived on zenodo: doi:[!](https://zenodo.org/badge/latestdoi/313364218)



## Results:

The regulatory relationships between yeast genes as inferred using this pipeline (for tests P2, P2P3, P2P5 and P) and published in [1] are made made available as gzipped csv files in (`data/predicted_networks`). The columns are in the following format:
   `regulator (name), target (name), weight (posterior probability)`.

## Data:

The following files are needed and should be put in the `data/input` folder:

### a. Expression data and genotypes from [3]:

The data from the supporting information of [3] at [ https://figshare.com/s/83bddc1ddf3f97108ad4 ].
The following files are used:
Data | Filename
-----|---------
expression data | SI_Data_01_expressionValues.txt.zip
covariates for expression data | SI_Data_02_covariates.xlsx 
genotypes | SI_Data_03_genotypes.txt.zip
eQTLs | SI_Data_04_eQTL.xlsx


### b. YEASTRACT ground truth data to compute precision and recall [4]:

Regulation Matrices can be  obtained from [ http://www.yeastract.com/formregmatrix.php ].
We retrieved the full ground-truth matrices containing all reported interactions of the following types from the YEASTRACT website: DNA binding evidence was used as the “Binding”, expression evidence including TFs acting as activators and those acting as inhibitors was used as the “Expression”, DNA binding and expression evidence was used as the “Binding & Expression”. Self regulation was removed from all ground truths. The matrices we retrieved are available as gzipped csv files in `data/input/yeastract`.

### c. Gene annotations from Ensembl [5]:

We use a file listing all genes, pseudogenes, etc. from Ensemble release 83: `Saccharomyces_cerevisiae.R64-1-1.83.gff3` .
The file should be processed with sed commands given in  `sed_processing_gff.sh` .
The result is a file where columns are separated by spaces, it contains gene name, start, end and a few more annotations.

## Required python packages:

This pipeline has been tested in python version 3.7.4.
The scripts requires Findr and the following packages:
   - numpy, pandas, 
   - statsmodels,
   - roman: to convert roman numerals to integers,
   - matplotlib and seaborn.


## Steps to run the analysis:

The scripts to run the analysis with Findr and to obtain binary causal networks for FDR thresholds given in [1].
The scripts should be run in the order they are numbered, the shell script "run_all.sh"  can
run them all, however *this may take a while*. **Therefore we recommend to run them in order**:
   - 1_select_strongest_cis_eqtls.py
   - 2_prepare_genotype_data.py
   - 3_reorder_expression_data.py
   - 4_covariate_regression_on_expression_data.py
   - 5_preprocessing_ensembl_data.sh
   - 6_run_findr.py
   - 7_select_yeastract_compatible_subset.py

The script in "subsamples" can be used to run Findr on randomly selected subsamples of the yeast data:
   - run_findr_subsampling_v9.py


## References:

1. Ludl, A-A and Michoel, T (2020) Comparison between instrumental variable and mediation-based methods for reconstructing causal gene networks in yeast
   (accepted)
   - arxiv: https://arxiv.org/abs/2010.07417
   - biorxiv: https://biorxiv.org/cgi/content/short/2020.10.13.337501v1

2. Findr paper:
    Wang, L and  Michoel, T (2017) PLoS Comput Biol 13(8): e1005703.
   - paper: https://doi.org/10.1371/journal.pcbi.1005703
   - source code: https://github.com/lingfeiwang/findr
   - python package: https://github.com/lingfeiwang/findr-python

3. Albert, F. W., Bloom, J. S., Siegel, J., Day, L., & Kruglyak, L. (2018). Genetics of trans-regulatory variation in gene expression. Elife, 7, e35471. doi:10.7554/eLife.35471
   - paper: https://doi.org/10.7554/elife.35471
   - data: https://figshare.com/s/83bddc1ddf3f97108ad4

4. Yeastract regulatory network:
    http://www.yeastract.com/formregmatrix.php

    YEASTRACT+: a portal for cross-species comparative genomics of transcription regulation in yeasts.
    Nucleic Acids Research, 48(D1):D642-D649   (doi:10.1093/nar/gkz859) 
    P.T. Monteiro, J. Oliveira, P. Pais, M. Antunes, M. Palma, M. Cavalheiro, M. Galocha, C.P. Godinho, L.C. Martins, N. Bourbon, M.N. Mota, R.A. Ribeiro, R.Viana, I. Sá-Correia, M.C. Teixeira (2020)
   - paper: https://doi.org/10.1093/nar/gkz859
    
5. Ensembl library for yeast (S. cerevisiae):
   - Ensemble release 83 (gff3 file): `ftp://ftp.ensembl.org/pub/release-83/gff3/saccharomyces_cerevisiae/`
   - Saccharomyces cerevisiae: http://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index?db=core
   - Ensembl Archives: http://www.ensembl.org/info/website/archives/index.html

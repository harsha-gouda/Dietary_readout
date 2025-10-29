# Dietary_readout from MS/MS food library
Classification based food biomarker extraction using a multilayer approach. Reference MS spectral matching of unknown samples with food library for dietary readout.

All the R code discribed in the folder "Code" has been tested on the windows11 laptop using the R version 4.4.1 using Rstudio. To obtain dietary scores for metabolomics data of your interest, please go to GNPS2 supported https://foodreadouts.gnps2.org/ and follow the SOP using GNPS workflows. 

Spectral files (mgf and csv) are deposited in Zenodo ([10.5281/zenodo.17469602](https://zenodo.org/records/17469602)) as files are large.

SOP for obtaining dietary readout from your datasets:

1. Run library matching workflow (GNPS2) with your mgf file of intrest matching against the food MS/MS library (500_foods_SPECTRUMID2.mgf)
2. Upload Library matching taskID and quantfile output into GNPS2 supported dietary readout web-based application (https://foodreadouts.gnps2.org/) and download dietary score results
3. Anlayze dietary score results in the statisitical analysis software of your interest or within the web-based application

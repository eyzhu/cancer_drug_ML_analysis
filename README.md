# Cancer drug machine learning analysis 

<img align="center" src="images/ML_fig1.jpg">

## Table of Contents

* [Overview](#overview)
* [Drug analysis](#drug-analysis)
* [Additional analysis](#additional-analysis)
* [Contributors](#contributors)

## Overview
* Raw scripts/Pipeline for analysis performed in manuscript __*Machine learning approach informs biology of cancer drug response*__.
* Analysis can be replicated by script provided in folders named after each drug.
* R version 3.5.1 was used for the analysis.

## Drug analysis
1. First download array expression data and place this file in the common_files folder. As of 4/17/2022, There's no live link to the file I used, so I uploaded it to dropbox. The file was orginally obtained from ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/sanger1018_brainarray_ensemblgene_rma.txt.gz. 
	*  https://www.dropbox.com/s/63d664sknfh8iv3/sanger1018_brainarray_ensemblgene_rma.txt?dl=0
2. Navigate to a numbered folder for a specific drug analysis, and run findFeatures.R to perform the analysis done in the paper.
	* The findFeatures.R script is broken down into sections as depicted in the schematic above. Section A performs the machine learning analysis, Section B performs the analysis using simple methods, and section C generates figures for the machine learning analysis.
	* Note that strAdjMat and importgenes variables (output of steps 2 and 3) are preloaded in the script to obtain the same results as the original analysis. These variables will slightly change depending on the version of string DB and permutations performed by Boruta.

## Contributors 
* Eliot Zhu
* https://github.com/eyzhu

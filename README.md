# TiSAn
Tissue-specific variant annotation

This project aims at providing a functional annotation tool for genetic variation, using prior on different tissues.
The general framework proposed with TiSAn will help researchers creating a predictive model for their tissue of interest.
Here, two models are made available as examples: the human brain and heart.

## Installation 
Run the following comand in R:

`devtools::install_github("kevinVervier/TiSAn")`

## Dependencies (R packages)
All packages should be installed automatically at the same time that TiSAN package.

## Databases
User can find human brain and heart databases, as gzipped .bed files (+index) at https://drive.google.com/open?id=0B8FEZgEQ9qH1WHRmTFVfZVd5Rjg

## Examples
To automatically annotate variants with TiSAn scores, we propose to use 'vcfanno' tool. Please note that most of the annotation tools will work with a datbase in bed format.
The following command calls vcfanno on a $VCF file containing positions to be annotated, and $CONFIG contains several parameters for vcfanno.

`vcfanno $CONFIG $VCF` 

Here is an example of CONFIG file:

`[[annotation]]`

`file="TiSAn_Brain.bed.gz"`

`columns=[4]`

`names=["TiSBrain"]`

`ops=["max"]`



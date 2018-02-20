This folder contains code to run TiSAn-build, the companion app for tissue-specific feature extraction.

# Panel1: expression Quantitative Loci (eQTL)

* In this panel, the user needs to provide a database of known tissue-specific eQTLs (e.g., GTEx consortium).
* Once a directory has been selected, the program will go through all the different tissues available in the database.
* A checkbox list will appear, and the user can select all the tissue types he considers as 'of interest' for a set of positive annotations.
* By defaut, the remaining tissues will be held as 'non-tissue' eQTLs.
* Once the tissue(s) of interest are selected in the checkbox list, the user can click on the 'eQTL tissue selection done' button.
* After a couple of minutes, two eQTL databases (tissue ('tissue_eqtl_gr.Rdata') and non-tissue ('non_tissue_eqtl_gr.Rdata')) should appear in the working directory.
* NB: if a file named 'tissue_eqtl_gr.Rdata' or 'non_tissue_eqtl_gr.Rdata' already exists, the program will not overwrite it.

# Panel2: differentially methylated regions (DMR)

* In this panel, the user needs to provide a database of known tissue-specific methylation profiles (e.g., RoadMap Epigenomics consortium).
* The program expects one single file containing all the methylation levels measured in different cell lines at different loci, and a metadata file contianing the information regarding the different sampled tissues.
* Once these two files have been provided, the program will go through all the different tissues available in the database.
* A checkbox list will appear, and the user can select all the tissue types he considers as 'of interest' for a set of positive annotations.
* By defaut, the remaining tissues will be held as 'non-tissue' DMRs.
* Once the tissue(s) of interest are selected in the checkbox list, the user can click on the 'Methylation tissue selection done' button.
* After a couple of minutes, a DMR database ('rme_gr.Rdata') should appear in the working directory.
* Here, it contains, for thousands of loci, the average DNA methylation level found in tissue and non-tissue cell lines.
* NB: if a file named 'rme_gr.Rdata' already exists, the program will not overwrite it.

# Panel3: Literature Mining (Pubmed)

* In this panel, the user needs to provide a database of known genes in their corresponding publications (e.g., gene2pubmed).
* The program expects one single file containing the pairs of genes and publications.
* The program also expects key words related to the tissue of interest (e.g., 'brain').
* Once it has been provided, the user can click on the 'tissue-specific keywords provided' button and the program will go through all the publications and search for keywords (it might take half an hour).
* Once the tissue-related genes are identified, they will be stored in 'tissue_gene_gr.Rdata'.
* The remaining genes will be stored in 'other_gene_gr.Rdata'.
* NB: if a file named 'tissue_gene_gr.Rdata' or 'other_gene_gr.Rdata' already exists, the program will not overwrite it.


# Panel4: Training set composition (LincSNP)

* In this panel, the user needs to provide a database of known disease-associated long intergenic non-coding RNAs (lincRNA) (e.g., LincSNP).
* The program expects one directory with one file per chromosome containing the pairs of SNPs and diseases.
* Once this directory has been provided, the program will go through all the different diseases available in the database.
* A checkbox list will appear, and the user can select all the tissue-related diseases he considers as 'of interest' for a set of positive annotations.
* Once it has been provided, the user can click on the 'lincSNP disease selection done' button and the program will go through all the known associations (might take 5 minutes).
* Once the disease-tissue SNPs are identified, they will be stored in 'train_lincsnp_pos.txt'.
* The remaining loci will be stored in 'train_lincsnp_neg.txt'.
* NB: if a file named 'train_lincsnp_pos.txt' or 'train_lincsnp_neg.txt' already exists, the program will not overwrite it.

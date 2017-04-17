---
title: "Create tissue-specific databases"
author: "Kevin Vervier"
date: "April 17, 2017"
output: html_document
---

In this vignette, we detail how tissue-specific databases were created using publicly available data sets. Results from this vignette are used in sources, like vignette 'estimate_distributions'.


# How to build tissue-related genes list using PubMed

Require to download gene to pubmed ID database (ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pumed.gz our version is from May 2016) and unzip it

Load Pubmed data mapping gene ID to publication ID
```{r}
pubmed = read.table('gene2pubmed')
# filter non-human organisms
pubmed = pubmed[pubmed$V1 == '9606',]
```
Query tissue, with a query in both titles and abstracts
```{r}
id=getIDs("your_tissue[Title] OR your_tissue[abstract]")
idx = which(pubmed$V3 %in% id) # map PubMed query to gene2pubmed table
genes.tissue = pubmed$V2[idx] # keep genes found in query
```

Create a table with gene ID, gen symbol and number of citations (decreasing order)

```{r}
library(org.Hs.eg.db)
genes.tissue.db = cbind(names(sort(table(as.character(genes.tissue)),decreasing=TRUE)[sort(table(as.character(genes.tissue)),decreasing=TRUE) > 4]),
                       annotate::getSYMBOL(names(sort(table(as.character(genes.tissue)),decreasing=TRUE)[sort(table(as.character(genes.tissue)),decreasing=TRUE) > 4]), data='org.Hs.eg')
                       ,sort(table(as.character(genes.tissue)),decreasing=TRUE)[sort(table(as.character(genes.tissue)),decreasing=TRUE) > 4])
colnames(genes.tissue.db) = c('GeneId','GeneSymbol','Citations')
```

Get corresponding gene locations

```{r}
library('biomaRt')
genes <- genes.tissue.db[,'GeneSymbol']
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37) #keep genome reference uniform !
genes_pos <- getBM(attributes=c('ensembl_gene_id','hgnc_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'hgnc_symbol', values = unique(genes), mart = ensembl)
genes_pos = genes_pos[order(genes_pos$chromosome_name,genes_pos$start_position),]

#filter all non-standard chromosome names
filt.gene.list= genes_pos[which(genes_pos$chromosome_name%in%c(1:22,'X','Y')),]
# realign with a match needed because not the same dimension
idx = match(genes.tissue.db[,'GeneSymbol'],filt.gene.list$hgnc_symbol)

if(length(which(is.na(idx)))>0){ 
  db = filt.gene.list[idx[-which(is.na(idx))],]
  db = cbind(db,genes.tissue.db[-which(is.na(idx)),'Citations'])
}else{
  db = filt.gene.list[idx,]
  db = cbind(db,genes.tissue.db[,'Citations'])
}
colnames(db)[7] = 'Citations'

# tissue gene database
write.table(db,file='your_tissue_genes.db',row.names = F,quote = F)
```

### How to build the non-tissue gene database

Here, use any database containing coding genes for hg19 (we used BRAINSPAN project list of genes: http://brainspan.org/api/v2/well_known_file_download/267666527)
```{r}
gene.info = read.csv('coding_genes_database',header=T)
# filter tissue genes
sub.db = gene.info[-which(gene.info$gene_symbol %in% db[,'hgnc_symbol']),]
# get their positions
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37) #keep genome reference uniform !
genes_pos <- getBM(attributes=c('ensembl_gene_id','hgnc_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'hgnc_symbol', values = unique(sub.db$gene_symbol), mart = ensembl)
genes_pos = genes_pos[order(genes_pos$chromosome_name,genes_pos$start_position),]
write.table(genes_pos,file='your_background_genes.db',row.names = F,quote = F)
```


# Create databases for tissue and non-tissue RoadMap Epigenomics(RME) methylation levels

Raw data can be downloaded at http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/DMRs/WGBS_DMRs_v2.tsv.gz

Code for splitting the RME consildated epigenomces file in two parts
```{r}
 system('cut -f 1,2,3,add_columns_corresponding_to_tissue WGBS_DMRs_v2.tsv >> WGBS_DMRs_v2_tissue.tsv') # column examples for Heart: 25,31,32,37,38
 system('cut -f 1-3,add_columns_corresponding_to_nontissue /WGBS_DMRs_v2.tsv >> WGBS_DMRs_v2_non_tissue.tsv') # non-heart columns: 24,26-30,33-36,39-43
```

Code for applying mean function to RME data
```{r}
#load tissue db
db.tissue = read.table(file='WGBS_DMRs_v2_tissue.tsv',header=T)
tmp = cbind(db.tissue[,1:3],apply(db.tissue,1,function(x)mean(as.numeric(x[-c(1:3)]),na.rm=T)))
tmp[which(is.nan(tmp[,4])),4] = NA
colnames(tmp) = c('chr','start','end','average_methylation')
write.table(tmp,file='WGBS_DMRs_v2_tissue_average.tsv',quote=F,row.names = F)
#load non-tissue db
db.bgd = read.table(file='WGBS_DMRs_v2_non_tissue.tsv',header=T)
tmp = cbind(db.bgd[,1:3],apply(db.bgd,1,function(x)mean(as.numeric(x[-c(1:3)]),na.rm=T)))
tmp[which(is.nan(tmp[,4])),4] = NA
colnames(tmp) = c('chr','start','end','average_methylation')
write.table(tmp,file='WGBS_DMRs_v2_non_tissue_average.tsv',quote=F,row.names = F)
```

# Create tissue and non-tissue GTEx eQTLs database

The raw data (V6) can be accessed through GTEx portal (http://gtexportal.org/home/datasets/GTEx_Analysis_V6_eQTLs.tar.gz)

```{r}
tissue.files= list.files('path_to_unzip_GTEx_files',pattern = 'your_tissue_name') #e.g., 'Heart'

#store SNPs info for tissue eQTL
snp.matrix = NULL
for(file in tissue.files){
  x = read.delim(file,header=T,as.is = T)
  snp.matrix = rbind(snp.matrix,cbind(x$rs_id_dbSNP142_GRCh37p13,x$snp_chrom,x$snp_pos,x$ref,x$alt,x$gene_name,x$tss_position,x$tss_distance,file))
}
colnames(snp.matrix) = c('rs_id_dbSNP142_GRCh37p13','snp_chrom','snp_pos','ref','alt','gene','tss_position','tss_distance','source')
write.table(snp.matrix,file='your_tissue_eQTL.db',quote=F,row.names = F)

#non-tissue eQTLs
files = list.files('path_to_unzip_GTEx_files',pattern = '.snpgenes')
# filter tissue files
files = files[-which(files %in% tissue.files)]
#store SNPs info for non-tissue eQTL
snp.matrix = NULL
for(file in files){
  x = read.delim(en.dir %&% file,header=T,as.is = T)
  snp.matrix = rbind(snp.matrix,cbind(x$rs_id_dbSNP142_GRCh37p13,x$snp_chrom,x$snp_pos,x$ref,x$alt,x$gene_name,x$tss_position,x$tss_distance,file))
}
colnames(snp.matrix) = c('rs_id_dbSNP142_GRCh37p13','snp_chrom','snp_pos','ref','alt','gene','tss_position','tss_distance','source')

write.table(snp.matrix,file='your_background_eQTL.db',quote=F,row.names = F)
```

# Create fetal heart enhancers database

```{r}
# load the catalog
db = read.csv('their_fetal_enhancer.db',header=TRUE) # example found here http://www.nature.com/article-assets/npg/ncomms/2016/161005/ncomms12923/extref/ncomms12923-s3.xlsx
# number of enhancers with PreNatal score higher than postNatal
length(which(db$score_Prenatal > db$score_Postnatal)) # 27,375
idx = which(db$score_Prenatal > db$score_Postnatal)

db = db[idx,]

# remove chr X and Y
db = db[-which(db$chromosome%in%c('chrX','chrY')),] # 349 non autosomal --> 27,026 candidate regions

# score > 0.2
#length(which(db$score_All > 0.2)) # 5,919 regions
length(which(db$score_Prenatal > 0.2)) # 8,991 regions

db = db[which(db$score_Prenatal > 0.2),]

#table(db$overlaps.gene.body.) # 2,927 that do not overlap gene body

############
# formatting 
db = db[,c(1,2,3,6,8,9,10)]
db$chromosome = gsub('chr','',db$chromosome)
db = db[order(db$chromosome,db$start,decreasing = FALSE),]

write.csv(db,file='your_fetal_enhancer.db',row.names = FALSE)
```

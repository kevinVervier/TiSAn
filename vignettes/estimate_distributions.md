---
title: "Estimate proximity distributions"
author: "Kevin Vervier"
date: "April 17, 2017"
output: html_document
---

In this vignette, we provide source to estimate distribution parameters we used as features in TiSAn models. Those parameters can be reported in the corresponding scripts, like 'R/RME_methyl.R'.

# ENCODE/RoadMap Epigenomics methylated regions

Estimate distribution for proximity with closest tissue and non-tissue RoadMap Epigenomics(RME) region.

Generate a random sample of 10,000 genomic positions
```{r}
bed = circlize::generateRandomBed(nr=10000)
chr = substr(bed$chr,4,6)
position = bed$start
```

Load tissue methylation values (from 'create_databases' vignette)
```{r}
db.tissue = read.table('your_tissue_RME.db',header=TRUE)
```
For every chr-pos pair, get the ones inside a methylated region (no check on chromosome)
```{r}
list.pos = lapply(1:length(chr), function(i) which(db.tissue$start <= position[i] & db.tissue$end >= position[i]))
```
For every chromosome check which methylated regions belong to it
```{r}
list.chr = lapply(c(1:24,'X','Y'), function(i) which(db.tissue$chr == i))
names(list.chr) = c(1:24,'X','Y')
```
Find which chr-pos pairs are inside a methylated region
```{r}
bool = sapply(1:length(chr), function(i) as.integer(as.logical(length(intersect(db.tissue$chr[list.pos[[i]]],chr[i])))))
```
For chr-pos pairs not in a methylated region, find the closest region
```{r}
options(warn=-1) # remove warnings for potentially infinite values
tmp = sapply(c(1:length(chr))[!bool], function(i) min(min(abs(position[i] - db.tissue$start[list.chr[[chr[i]]]])), min(abs(position[i] - db.tissue$end[list.chr[[chr[i]]]])))) # Inf value means Chr Y  dist = rep(0,length(position))
options(warn=0)
#init distance vector
dist = rep(0,length(position))
# raw bp distance
if(length(tmp) > 0) dist[!bool] = tmp 
```
Fit Weibull distribution after removing distribution tail
```{r}
tmp = dist[which(dist < 10000)]
fit <- MASS::fitdistr(tmp + 1 , "weibull")$estimate
# shape and scale need to be pasted in 'RME_methyl.R'
```  

# Proximity with tissue and non-tissue genes
Estimate distance to closest tissue gene distribution.
Load databases generated in 'create_databases' vignette
```{r}
tissue.gene.list = read.table('your_tissue_genes.db',header=T)
bgd.gene.list = read.table('your_background_genes.db',header=T)
```
Generate a sample of 10,000 random genomic positions
```{r}
bed = circlize::generateRandomBed(nr=10000)
chr = substr(bed$chr,4,6)
position = bed$start
```
For each position, find which tissue gene is the closest (from gene start or end)
```{r}
tissue.idx = sapply(1:length(chr), function(i) which(tissue.gene.list$chromosome_name == chr[i]))
x = sapply(1:length(chr),function(i) min(abs(c(tissue.gene.list$end_position[tissue.idx[[i]]],tissue.gene.list$start_position[tissue.idx[[i]]] ) - position[i])))
````
Filter distribution tail (e.g., positions far from any gene)
```{r}
tmp = x[x < 10000000]
```
Fit a Weibull distribution and get scale and shape parameters
```{r}
MASS::fitdistr(tmp + 1 , "weibull")$estimate # parameters shape and scale go in getClosestGene.R in tissue genes section
```

For each position, find which non-tissue gene is the closest (from gene start or end)
```{r}
bgd.idx = sapply(1:length(chr), function(i) which(bgd.gene.list$chromosome_name == chr[i]))
x = sapply(1:length(chr),function(i) min(abs(c(bgd.gene.list$end_position[bgd.idx[[i]]],bgd.gene.list$start_position[bgd.idx[[i]]] ) - position[i])))
```
Filter distribution tail (e.g., positions far from any gene)
```{r}
tmp = x[x < 10000000]
```
Fit a Weibull distribution and get scale and shape parameters
```{r}
MASS::fitdistr(tmp + 1 , "weibull")$estimate # shape and scale go in getClosestGene.R in non-tissue genes section
```

# Gene-Tissue Expression (GTEx) expression quantitative trait loci (eQTL)

Estimate distance to closest GTEx eQTLs for tissue and non-tissue databases

Get Weibull coefficients for tissue eQTLs

```{r}
file = 'your_tissue_eQTL.db' #here is the path to the tissue-related GTEx eQTL database (created by vignette 'create_databases')
x=read.table(file,header=T)
# get the closest distance between a gene and an eQTL
tmp = tapply(x$tss_distance,x$gene,function(x) min(abs(x)))
# add 1 because Weibull distribution is not defined for 0 values
tmp = tmp + 1
# estimate distribution parameters
MASS::fitdistr(tmp, "weibull")$estimate #
```

Get Weibull coefficients for non-tissue eQTLs
```{r}
file = 'your_background_eQTL.db'#here is the path to the non-tissue related GTEx eQTL database (created by vignette 'create_databases')
x=read.table(file,header=T)
# get the closest distance between a gene and an eQTL
tmp = tapply(x$tss_distance,x$gene,function(x) min(abs(x)))
# add 1 because Weibull distribution is not defined for 0 values
tmp = tmp + 1
# estimate distribution parameters
MASS::fitdistr(tmp4, "weibull")$estimate #
```
Those parameters need to be reported in R/GTEx_eQTL.R script in shape.tissue, scale.tissue, shape.other, scale.other

# Developmental brain methylation loci

Estimate distance to developmental brain methylation loci
Load data for dDMRs and dDMPs (Speirs et al. 2015)
```{r}
db.regions = read.csv(file='Bray2015.csv',header=T) 
db.pos = read.csv(file='dDMPs.csv',header=T)
```
Generate a random sample of 10,000 genomics positions
```{r}
bed = circlize::generateRandomBed(nr=10000)
chr = substr(bed$chr,4,6)
position = bed$start
```

Estimate distance to dDM single positions

```{r}
x = sapply(1:length(chr), function(i) min(abs(db.pos$Position[which(db.pos$Chr == chr[i])] - position[i])))
signes = sapply(1:length(chr), function(i) which.min(abs(db.pos$Position[which(db.pos$Chr == chr[i])] - position[i])))
y = unlist(sapply(1:length(chr), function(i) db.pos$Position[which(db.pos$Chr == chr[i])[signes[[i]]]]-position[i]))

tmp = x[x < 1000000]
MASS::fitdistr(tmp + 1 , "weibull")$estimate # those parameters are the position-related parameters in 'dDMPs_dist.R'
```

Estimate distance to dDM regions
```{r}
bool.reg =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.regions$chrom == chr[i] & db.regions$start <= position[i] & db.regions$end >= position[i]))))) 
x = sapply(c(1:length(chr))[!bool.reg], function(i) min( min(abs(db.regions$start[which(db.regions$chrom == chr[i])] - position[i])), min(abs(position[i] - db.regions$end[which(db.regions$chrom == chr[i])]))))
tmp = x[x < 10000000]
MASS::fitdistr(tmp + 1 , "weibull")$estimate # those parameters are region-related in 'dDMPs_dist.R'
```

# Fetal heart developmental enhancer site 

```{r}
#load db
db.regions = read.csv(file='your_fetal_enhancer.db',header=T)

#random sample
bed = circlize::generateRandomBed(nr=10000)
chr = substr(bed$chr,4,6)
position = bed$start

#compute region dist
bool.reg =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.regions$chromosome == chr[i] & db.regions$start <= position[i] & db.regions$end >= position[i]))))) 
x = sapply(c(1:length(chr))[!bool.reg], function(i) min( min(abs(db.regions$start[which(db.regions$chromosome == chr[i])] - position[i])), min(abs(position[i] - db.regions$end[which(db.regions$chromosome == chr[i])]))))
tmp = x[x < 10000000]
fit <- fitdistr(tmp + 1 , "weibull")$estimate # those parameters are the position-related parameters in 'fetal_cardiac_enhancers_dist.R'

```

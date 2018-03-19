---
title: 'Training TiSAn models: the pipeline'
author: "Kevin Vervier"
date: "March 2, 2018"
output: html_document
---

In this vignette, we detail how tissue-specific models are trained, from collecting training data and annotations to make genome-wide predictions.

# How to select training set loci

First thing first, the user is asked to define which genomic positions are related to the tissue of interest.
We recommend to have at least 10,000 positive and 10,000 negative examples.
Many databases could be used to gather those loci. In the current manuscript, we focused on two kinds: genotype arrays and non-coding RNA.

* Genotype arrays designed by large consortia are usually screening variants with strong disease associations (e.g., PsychArray for psychiatric disorders, or MetaboChip for cardiovascular diseases). 
* Large intergenic non-coding RNAs catalogs, such as LincSNP, provide important resources for loci out of the genic regions.

Users could refer to the TiSAn-build application (last panel: lincSNP) to extract positions based on a list of tissue-related disorders/traits. For additional details on the process, users could also check the vignette 'create_training_set'.
Users could have additional sources for training data (e.g., genotype array), and simply need to name them with 'pos' and 'neg' in their names.
After this step, a training set is saved, with at least column 1 being the chromosome number and column 2 is the variant location.

In this section we load the entire set of training positions across the different training files.
```{r}
data.loc='../TiSAn-build' # default location of the training files is in TiSAn-build (can be changed).
#get all positive training files
pos.files <- list.files(path=data.loc,pattern = 'train.*pos',full.names = TRUE)
#get all negative training files
neg.files <- list.files(path=data.loc,pattern = 'train.*pos',full.names = TRUE)

# merge and create a label tag (0 if negative example and 1 if positive example)
X = NULL
Y = NULL
for(pos.file in pos.files){
  tmp <- read.delim(pos.file)[,1:2]
  X = rbind(X,tmp)
  Y = c(Y,rep(1,nrow(tmp)))
}
for(neg.file in neg.files){
  tmp <- read.delim(neg.file)[,1:2]
  X = rbind(X,tmp)
  Y = c(Y,rep(0,nrow(tmp)))
}
save(X,Y,file='training_set.Rdata')
```

# How to prepare feature extraction (distance to loci)

In this section, we estimate the parameters of the Weibull distance for each single annotation (e.g., eQTL, methylation, developmental, ...). For additional details on the process, users could also check the vignette 'estimate_distributions'.
We assume that all annotation databases are stored in files with 'gr' in their name.
After this step, you should get a Rdata file with all the estimated parameters saved in it.
```{r}
# need MASS package for distribution fitting
if(!require(MASS)) install.packages('MASS')
# need circlize package to sample genomic positions
if(!require(circlize)) install.packages('circlize')
# need GenomicRanges package for Granges object
if(!require(GenomicRanges)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}
# Generate a random sample of genomic positions
bed = circlize::generateRandomBed(nr=10000)
chr = substr(bed$chr,4,6)
position = bed$start
random_gr = GRanges(seqnames = Rle(chr),
                    ranges = IRanges(start=position, end=position))

data.loc='../TiSAn-build' # default location of the training files is in TiSAn-build (can be changed).
#get all positive training files
anno.files <- list.files(path=data.loc,pattern = 'gr',full.names = TRUE)
# init a database to store parameters
param_db <- NULL
# loop over each annotation and estimate Weibull distribution parameters
for(f in anno.files){
  # load file
  bar <- load(f)
  #in case of multiple db are stored in the same object (not recommended)
  for(obj in bar){
    # compute the distance to the nearest element for each random locus
    dist <- distanceToNearest(random_gr,get(bar))
    # only keep distance value
    dist <- dist@elementMetadata$distance
    # remove tail values
    tmp = dist[which(dist < 10000)]
    fit <- MASS::fitdistr(tmp + 1 , "weibull")$estimate
    param_db <- rbind(param_db,c(obj,fit[1],fit[2]))
  }
}
colnames(param_db) = c('anno','shape','scale')
# save all the distribution parameters
save(param_db,file='weibull_param.Rdata')
```


# How to extract features for the training set

In this section, we annotate each training locus with the different database signals, and format it as a feature representation, used as such by most of the machine learning algorithms.

```{r}
# load training positions
load('training_set.Rdata')
# source the feature extraction functions
source('../R/feature_extraction.R')

feat = feature_extraction(chr=X[,1],position=X[,2])
# save the extracted features
save(X,Y,feat,file='training_set.Rdata')

```
NB: in this example, we omitted the extraction of the 'compositional' features (kmers, centromere and telomere), as they do no depend of any tissue-specific database.

# How to train a predictive model

In this section, we use the features previously extracted to fit a random forest model. Please note that we optimize the 2 parameters (mtry and ntree) through cross-validation.

```{r}
if(!require(randomForest)) install.packages('randomForest')
if(!require(pROC)) install.packages('pROC')
# load training data
load('training_set.Rdata')

#param grid
NTREE = c(10,100,500,1000)
MTRY = c(1,3,6,9,ncol(feat))
nrep = 10

#store perfs
best_ntree = NULL
best_auc = NULL

for(mtry in MTRY){
  repeats = matrix(0,nrow = nrep, ncol = length(NTREE))
  #fix the random seed
  set.seed(42)
  for(repe in 1:nrep){
    cat("Repeat",repe,"\n")
    nfolds = 5
    folds = sample(1:nfolds,replace=TRUE,length(Y))
    
    #loop over multiple costs
    cpt = 0
    
    auc = rep(0,length(NTREE))
    for(ntree in NTREE){
      #proceed CV
      preds = rep(0,length(Y))
      for(fold in 1:nfolds){
        cat('.')
        idx = which(folds == fold)
        X.train = feat[-idx,]
        Y.train = Y[-idx]
        X.test = feat[idx,]
        Y.test = Y[idx]
        
        #need to scale data
        s=scale(X.train,center=TRUE,scale=TRUE)
        s[which(is.na(s))] = 0
        # Scale the test data
        s2=scale(X.test,attr(s,"scaled:center"),attr(s,"scaled:scale"))
        s2[which(is.na(s2))] = 0
        s2[which(s2 == Inf)] = 0
        #train model
        m <- randomForest(x=s,y=as.factor(Y.train),ntree=ntree,mtry=mtry,strata=as.factor(Y.train),sampsize = rep(min(table(Y.train)),2))
        x=predict(m,s2,type='prob')
        preds[idx] = x[,'1']
        #make preds
      }
      cat('\n')
      #compute acc
      cpt = cpt + 1
      auc[cpt] = auc(Y,preds)
    }
    #remove NA
   if(length(which(is.na(auc))) > 0) auc[which(is.na(auc))] = 0
    repeats[repe,] = auc
  }
  colnames(repeats) = NTREE
  #find the best parameter
  cat('The best model for mtry=',mtry,'is for ntree=',colnames(repeats)[which.max(apply(repeats,2,mean))],'with a mean auc equal to',max(apply(repeats,2,mean)),'\n')
  best_auc = c(best_auc,max(apply(repeats,2,mean)))
  best_ntree = c(best_ntree,colnames(repeats)[which.max(apply(repeats,2,mean))])
}

# define optimal parameters
opti_ntree = as.numeric(best_ntree[which.max(best_auc)])
best_mtry = MTRY[which.max(best_auc)]


########################################
# get final cross-validated predictions

set.seed(42)
nfolds = 5
folds = sample(1:nfolds,replace=TRUE,length(Y))
preds = rep(0,length(Y))
for(fold in 1:nfolds){
  cat('.')
  idx = which(folds == fold)
  X.train = feat[-idx,]
  Y.train = Y[-idx]
  X.test = feat[idx,]
  Y.test = Y[idx]
  #need to scale data
  s=scale(X.train,center=TRUE,scale=TRUE)
  s[which(is.na(s))] = 0
  # Scale the test data
  s2=scale(X.test,attr(s,"scaled:center"),attr(s,"scaled:scale"))
  s2[which(is.na(s2))] = 0
  s2[which(s2 == Inf)] = 0
  #train model
  m <- randomForest(x=s,y=as.factor(Y.train),ntree=opti_ntree,mtry=best_mtry,strata=as.factor(Y.train),sampsize = rep(min(table(Y.train)),2))
  x=predict(m,s2,type='prob')
  preds[idx] = x[,'1']
}
save(preds,file='cv_predictions.Rdata')


############
# Final model + scaling
X.train = feat
Y.train = Y
s=scale(X.train,center=TRUE,scale=TRUE)
s[which(is.na(s))] = 0
m <- randomForest(x=s,y=as.factor(Y.train),ntree=opti_ntree,mtry=best_mtry,strata=as.factor(Y.train),sampsize = rep(min(table(Y.train)),2))  
save(m,s,file='model.Rdata')

```

# How to define the optimal score threshold

In this section, we cover how to define a hard threshold value on the tissue-specific scores, based on a false positive rate of 10%.

```{r}
# load cross-validated performances
load('cv_predictions.Rdata')

# Find thresh value for 1% FPR
FPR = c()

for(THRESH in seq(0,1,by=0.01)){
  
  db = data.frame(Y=Y,proba=preds,TRS=(1-THRESH)/THRESH-(1-preds)/(preds)) # scaled odd-ratio formula
  db$TRS[db$TRS < 0] = 0
  db$TRS = db$TRS/((1-THRESH)/THRESH)

  FPR = c(FPR,length(which(db$TRS[db$Y == 0] > 0 )) / length(which(db$Y == 1))) 

}

THRESH = seq(0,1,by=0.01)[which.min(abs(FPR-0.1))] # here the value is 0.84
save(m,s,THRESH,file='model.Rdata')
```

Now, you could use the saved objects to score any genomic position using the 'get_predictions' vignette.

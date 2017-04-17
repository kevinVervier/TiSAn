#------------------------------------------#
# Tissue genes related features extraction #
#------------------------------------------#

#' Compute proximity with closest tissue and non-tissue genes
#'
#' \code{getClosestGene}
#' @param chr the chromosome number
#' @param position the locus position
#' @return tissue.gene.bool: boolean value equal to TRUE if the position is inside a known tissue gene
#' @return tissue.gene.prox: proximity value to the closest tissue gene
#' @return bgd.gene.bool: boolean value equal to TRUE if the position is inside a non-tissue gene
#' @return bgd.gene.prox: proximity value to the closest non-tissue gene
#' @export
#' @examples #get genes related features for locus chr1:5000000
#' getClosestGene(chr=1,position=5000000)
#'

getClosestGene <- function(chr,position){

  # convert chr and position in string if needed
  if(is.factor(chr)) chr = as.character(chr)
  if(is.factor(position)) position = as.numeric(position)

  db.tissue = read.table('your_tissue_genes.db',header=TRUE) #here is the path to the tissue-related genes database (created by vignette 'create_databases')
  db.bgd = read.table('your_background_genes.db',header=TRUE) #here is the path to the tissue-related genes database (created by vignette 'create_databases')
  #main loop
  if(length(chr) == length(position)){
    ##############
    #tissue genes#
    bool.tissue =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.tissue$chromosome_name == chr[i] & db.tissue$start_position <= position[i] & db.tissue$end_position >= position[i])))))
    tmp = sapply(c(1:length(chr))[!bool.tissue], function(i) min( min(abs(db.tissue$start_position[which(db.tissue$chromosome_name == chr[i])] - position[i])), min(abs(position[i] - db.tissue$end_position[which(db.tissue$chromosome_name == chr[i])]))))
    dist.tissue = rep(0,length(position))
    if(length(tmp) != 0 ) dist.tissue[!bool.tissue] = tmp
    #use weibull estimated parameters
    shape = 0.6944 # here is the tissue-data Weibull shape parameter from 'estimate_distributions' vignette
    scale = 691241 # here is the tissue-data Weibull scale parameter from 'estimate_distributions' vignette
    dist.tissue = (shape/scale) * (dist.tissue/scale)^(shape-1) * exp(-(dist.tissue/scale)^shape) #
    dist.tissue[is.infinite(dist.tissue)] = (shape/scale) * (1/scale)^(shape-1) * exp(-(1/scale)^shape) #compute the max value for dist = 1

    ##################
    #non-tissue genes#
    bool.bgd =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.bgd$chromosome_name == chr[i] & db.bgd$start_position <= position[i] & db.bgd$end_position >= position[i])))))
    tmp = sapply(c(1:length(chr))[!bool.bgd], function(i) min( min(abs(db.bgd$start_position[which(db.bgd$chromosome_name == chr[i])] - position[i])), min(abs(position[i] - db.bgd$end_position[which(db.bgd$chromosome_name == chr[i])]))))
    dist.bgd = rep(0,length(position))
    if(length(tmp) != 0 ) dist.bgd[!bool.bgd] = tmp
    #use weibull estimated parameters
    shape =  0.549 # here is the non-tissue-data Weibull shape parameter from 'estimate_distributions' vignette
    scale= 230731 # here is the non-tissue-data Weibull scale parameter from 'estimate_distributions' vignette
    dist.bgd = (shape/scale) * (dist.bgd/scale)^(shape-1) * exp(-(dist.bgd/scale)^shape) #
    dist.bgd[is.infinite(dist.bgd)] = (shape/scale) * (1/scale)^(shape-1) * exp(-(1/scale)^shape) #compute the max value for dist = 1
    #output
    x = cbind(bool.tissue,dist.tissue,bool.bgd,dist.bgd)
    colnames(x) = c('tissue.gene.bool','tissue.gene.prox','bgd.gene.bool','bgd.gene.prox')
    return(x)
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}

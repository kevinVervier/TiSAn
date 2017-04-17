#---------------------------------------------------------------------------#
# Source for RME methylation level data for tissue and non-tissue databases #
#---------------------------------------------------------------------------#

#' Compute proximity with closest tissue and non-tissue RoadMap Epigenomics(RME) database
#'
#' \code{RME_methyl}
#' @param chr the chromosome number
#' @param position the locus position
#' @return RME.bool: boolean value equal to TRUE if the position is inside a RME methylated region
#' @return RME.prox: proximity value to the closest RME methylated region
#' @return RME.tissue: mean methylation value across all tissues of interest  (only if belongs to a methylated region)
#' @return RME.others: mean methylation value across all other tissues (only if belongs to a methylated region)
#' @export
#' @examples #get tissue-specific methylation related features for locus chr1:5000000
#' RME_methyl(chr=1,position=5000000)
#'

RME_methyl<-function(chr,position,verbose=T){

  # convert chr and position in string if needed
  if(is.factor(chr)) chr = as.character(chr)
  if(is.factor(position)) position = as.numeric(position)

  #load tissue methylation values (from 'create_RME_dataset.R')
  db.tissue = read.table('your_tissue_RME.db',header=TRUE)
  #load non-tissue methylation values
  db.bgd = read.table('your_background_RME.db',header=TRUE)

  #check if the number of provided chr is equal to the number of provided positions
  if(length(chr) == length(position)){
    ############
    # DISTANCE #
    ############

    #for every chr-pos pair, get the ones inside a methylated region (no check on chromosome)
    list.pos = lapply(1:length(chr), function(i) which(db.tissue$start <= position[i] & db.tissue$end >= position[i]))
    #for every chromosome check which methylated regions belong to it
    list.chr = lapply(c(1:24,'X','Y'), function(i) which(db.tissue$chr == i))
    names(list.chr) = c(1:24,'X','Y')

    #find which chr-pos pairs are inside a methylated region
    bool = sapply(1:length(chr), function(i) as.integer(as.logical(length(intersect(db.tissue$chr[list.pos[[i]]],chr[i])))))
    #for chr-pos pairs not in a methylated region, find the closest region
    options(warn=-1)
    tmp = sapply(c(1:length(chr))[!bool], function(i) min(min(abs(position[i] - db.tissue$start[list.chr[[chr[i]]]])), min(abs(position[i] - db.tissue$end[list.chr[[chr[i]]]])))) # Inf value means Chr Y  dist = rep(0,length(position))
    options(warn=0)
    #init distance vector
    dist = rep(0,length(position))
    #bp distance
    if(length(tmp) > 0) dist[!bool] = tmp
    #Weibull estimated parameters (found in 'estimate_distributions' vignette)
    # default values are from Heart model
    shape = 0.4375
    scale = 316
    # use Weibull formula to compute distance
    dist = (shape/scale) * (dist/scale)^(shape-1) * exp(-(dist/scale)^shape)
    #for postions inside region, compute the max value
    dist[is.infinite(dist)] = (shape/scale) * (1/scale)^(shape-1) * exp(-(1/scale)^shape) #compute the max value for dist = 0

    #####################
    # METHYLATION LEVEL #
    #####################

    #output matrix
    features = matrix(NA,nrow=length(chr),ncol=2)
    #get mean methylation values
    if(any(bool) != 0){
      #get methylation level for positions only if they belongs to a methylated region
      idx = sapply(which(as.logical(bool)), function(i) intersect(list.chr[[chr[i]]] , list.pos[[i]]))
      #merge tissue methylation and non-tissue methylation
      RME.list =  cbind(db.tissue[idx,4],db.bgd[idx,4])
      #put corresponding values for positions inside a methylated region
      features[as.logical(bool),] = RME.list
    }
    #merge boolean, distance and methylation level for output
    x = cbind(bool,dist,features)
    colnames(x) = c('RME.bool','RME.prox','RME.tissue','RME.others')
    return(x)
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}

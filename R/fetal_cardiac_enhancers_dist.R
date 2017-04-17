#---------------------------------------------#
# Fetal cardiac enhancers features extraction #
#---------------------------------------------#

#' Compute proximity with closest tissue and non-tissue genes
#'
#' \code{dvpt_enhancer_dist}
#' @param chr the chromosome number
#' @param position the locus position
#' @return dvpt.reg.bool: boolean value equal to TRUE if the position is inside a known fetal heart enhancer
#' @return dvpt.reg.prox: proximity value to the closest fetal heart enhancer
#' @export
#' @examples #get fetal heart enhancer features for locus chr1:5000000
#' dvpt_enhancer_dist(chr=1,position=5000000)
#'

dvpt_enhancer_dist<-function(chr,position){

  #load db
  db.regions = read.csv(file='your_fetal_enhancer.db',header=TRUE) #here is the path to the fetal heart enchancers database (created by vignette 'create_databases')

  if(length(chr) == length(position)){
    #regions related
    bool.reg =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.regions$chromosome == chr[i] & db.regions$start <= position[i] & db.regions$end >= position[i])))))
    options(warn=-1)
    tmp = sapply(c(1:length(chr))[!bool.reg], function(i) min( min(abs(db.regions$start[which(db.regions$chromosome == chr[i])] - position[i])), min(abs(position[i] - db.regions$end[which(db.regions$chromosome == chr[i])]))))
    options(warn=0)
    dist.reg = rep(0,length(position))
    dist.reg[!bool.reg] = tmp
    #use Weibull estimated paramters
    shape = 0.641689 # here is the Weibull shape parameter from 'estimate_distributions' vignette
    scale = 505640 # here is the Weibull scale parameter from 'estimate_distributions' vignette
    dist.reg = (shape/scale) * (dist.reg/scale)^(shape-1) * exp(-(dist.reg/scale)^shape)
    dist.reg[is.infinite(dist.reg)] = (shape/scale) * (1/scale)^(shape-1) * exp(-(1/scale)^shape) #compute the max value for dist = 1
    #output
    x = cbind(bool.reg,dist.reg)
    colnames(x) = c('dvpt.reg.bool','dvpt.reg.prox')
    return(x)
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}

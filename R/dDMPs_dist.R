#----------------------------------------------------------------------#
# Source fetching information on developmental brain methylation spots #
#----------------------------------------------------------------------#

#' Compute proximity with closest developmental brain methylation spots (regions or positions)
#'
#' \code{dDMPs_dist}
#' @param chr the chromosome number
#' @param position the locus position
#' @return dDMP-pos-bool: boolean value equal to TRUE if the position is a known dDMP position
#' @return dDMP-pos-prox: proximity value to the closest dDMP position
#' @return dDMP-reg-bool: boolean value equal to TRUE if the position is a known dDMP region
#' @return dDMP-reg-prox: proximity value to the closest dDMP region
#' @export
#' @examples #get dDMPs related features for locus chr1:5000000
#' dDMPs_dist(chr=1,position=5000000)
#'

#########################################################################
# function for distance/proximity on dDMPs
dDMPs_dist<-function(chr,position,verbose=T){

  #load db
  db.regions = read.csv(file='Bray2015.csv',header=TRUE) # 4,825 developmentally differentially methylated regions, from Spiers et al 2015
  db.pos = read.csv(file='dDMPs.csv',header=TRUE) # 28,718 developmentally differentially methylated positions from Spiers et al., 2015

  if(length(chr) == length(position)){
    ####################
    # position related #
    bool.pos =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.pos$Chr == chr[i] & db.pos$Position == position[i])))))
    options(warn=-1)
    tmp = sapply(c(1:length(chr))[!bool.pos], function(i) min(abs(db.pos$Position[which(db.pos$Chr == chr[i])] - position[i])))
    options(warn=0)
    dist.pos = rep(0,length(position))
    dist.pos[!bool.pos] = tmp

    #use Weibull estimated parameters (from 'estimate_distributions' vignette script)
    shape = 0.7612
    scale = 129555
    dist.pos = (shape/scale) * (dist.pos/scale)^(shape-1) * exp(-(dist.pos/scale)^shape) #
    dist.pos[is.infinite(dist.pos)] = (shape/scale) * (1/scale)^(shape-1) * exp(-(1/scale)^shape) #compute the max value for dist = 1

    ###################
    # regions related #
    bool.reg =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(db.regions$chrom == chr[i] & db.regions$start <= position[i] & db.regions$end >= position[i])))))
    options(warn=-1)
    tmp = sapply(c(1:length(chr))[!bool.reg], function(i) min( min(abs(db.regions$start[which(db.regions$chrom == chr[i])] - position[i])), min(abs(position[i] - db.regions$end[which(db.regions$chrom == chr[i])]))))
    options(warn=0)
    dist.reg = rep(0,length(position))
    dist.reg[!bool.reg] = tmp
    #use weibull estimated parameters (from 'estimate_distributions' vignette)
    shape = 0.726
    scale = 1074209
    dist.reg = (shape/scale) * (dist.reg/scale)^(shape-1) * exp(-(dist.reg/scale)^shape) #
    dist.reg[is.infinite(dist.reg)] = (shape/scale) * (1/scale)^(shape-1) * exp(-(1/scale)^shape) #compute the max value for dist = 1
    #output
    x = cbind(bool.pos,dist.pos,bool.reg,dist.reg)
    colnames(x) = c('dDMP-pos-bool','dDMP-pos-prox','dDMP-reg-bool','dDMP-reg-prox')
    return(x)
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}

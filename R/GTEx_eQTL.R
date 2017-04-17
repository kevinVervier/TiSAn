#---------------------------------------#
# GTEX eQTL related features extraction #
#---------------------------------------#

#' Compute proximity with closest tissue and non-tissue GTEx eQTLs
#'
#' \code{GTEx_eQTL}
#' @param chr the chromosome number
#' @param position the locus position
#' @return bool.tissue: boolean value equal to TRUE if the position is a known tissue GTEx eQTL
#' @return prox.tissue: proximity value to the closest tissue eQTL
#' @return bool.other: boolean value equal to TRUE if the position is a known non-tissue GTEx eQTL
#' @return prox.other: proximity value to the closest non-tissue eQTL
#' @export
#' @examples #get GTEx eQTL related features for locus chr1:5000000
#' GTEx_eQTL(chr=1,position=5000000)
#'


GTEx_eQTL <- function(chr,position){

  #load tissue database
  db.tissue = read.table('your_tissue_eQTL.db',header=TRUE)#here is the path to the tissue-related GTEx eQTL database (created by vignette 'create_databases')
  #load non-tissue database
  db.bgd = read.table('your_background_eQTL.db',header=TRUE)#here is the path to the non-tissue related GTEx eQTL database (created by vignette 'create_databases')
  #infered parameters from Weibull distribution on distance to the closest eQTL
  #here: default values are found for Heart model
  shape.tissue = 3.333712e-01 # here is the tissue-data Weibull shape parameter from 'estimate_distributions' vignette
  scale.tissue = 1.522949e+04 # here is the tissue-data Weibull scale parameter from 'estimate_distributions' vignette
  shape.other = 3.156386e-01 # here is the non-tissue-data Weibull shape parameter from 'estimate_distributions' vignette
  scale.other = 9.153087e+03 # here is the non-tissue-data Weibull scale parameter from 'estimate_distributions' vignette
  # check if the two inputs have the same size
  if(length(chr) == length(position)){
    # get tissue eQTL positions that are also in input positions
    tissue.list.pos = lapply(1:length(chr), function(i) which(db.tissue$snp_pos == position[i]))
    # get tissue eQTL chromosomes that are also in input chromosomes
    tissue.list.chr = lapply(c(1:24,'X','Y'), function(i) which(db.tissue$snp_chr == i))
    names(tissue.list.chr) = c(1:24,'X','Y')
    ##########
    # tissue #
    # boolean variable checking if an input position is a known tissue eQTL
    bool.tissue = sapply(1:length(chr), function(i) as.integer(as.logical(length(intersect(db.tissue$snp_chrom[tissue.list.pos[[i]]],chr[i])))))
    tmp = sapply(c(1:length(chr))[!bool.tissue], function(i) min(abs(position[i] - db.tissue$snp_pos[tissue.list.chr[[chr[i]]]]))) # Inf value means Chr Y
    # compute distance to closest tissue eQTL, using Weibull parameters
    dist.tissue = rep(0,length(position))
    dist.tissue[!bool.tissue] = tmp
    dist.tissue = (shape.tissue/scale.tissue) * (dist.tissue/scale.tissue)^(shape.tissue-1) * exp(-(dist.tissue/scale.tissue)^shape.tissue)
    # replace infinite by maximal proximity value, when position is a known eQTL
    dist.tissue[is.infinite(dist.tissue)] = (shape.tissue/scale.tissue) * (1/scale.tissue)^(shape.tissue-1) * exp(-(1/scale.tissue)^shape.tissue)

    #################
    # other tissues #
    # get tissue eQTL positions that are also in input positions
    other.list.pos = lapply(1:length(chr), function(i) which(db.bgd$snp_pos == position[i]))
    # get tissue eQTL chromosomes that are also in input chromosomes
    other.list.chr = lapply(c(1:24,'X','Y'), function(i) which(db.bgd$snp_chr == i))
    names(other.list.chr) = c(1:24,'X','Y')
    # boolean variable checking if an input position is a known non-tissue eQTL
    bool.other = sapply(1:length(chr), function(i) as.integer(as.logical(length(intersect(db.bgd$snp_chrom[other.list.pos[[i]]],chr[i])))))
    tmp = sapply(c(1:length(chr))[!bool.other], function(i) min(abs(position[i] - db.bgd$snp_pos[other.list.chr[[chr[i]]]]))) # Inf value means Chr Y
    # compute distance to closest non-tissue eQTL, using Weibull parameters
    dist.other = rep(0,length(position))
    dist.other[!bool.other] = tmp
    dist.other = (shape.other/scale.other) * (dist.other/scale.other)^(shape.other-1) * exp(-(dist.other/scale.other)^shape.other)
    # replace infinite by maximal proximity value, when position is a known eQTL
    dist.other[is.infinite(dist.other)] = (shape.other/scale.other) * (1/scale.other)^(shape.other-1) * exp(-(1/scale.other)^shape.other)

    return(data.frame('bool.tissue'=bool.tissue,'prox.tissue'=dist.tissue, 'bool.other'=bool.other,'prox.other'=dist.other))
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}

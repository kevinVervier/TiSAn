#-----------------------------------------------------#
# Centromere and Telomere related features extraction #
#-----------------------------------------------------#

#' Compute proximity with closest telomere and centromere
#' 
#' \code{centroTelo_dist}
#' @param chr the chromosome number
#' @param position the locus position
#' @return telomere.bool: boolean value equal to TRUE if the position is inside a telomere region
#' @return telomere.prox: proximity value to the closest telomere
#' @return centromere.bool: boolean value equal to TRUE if the position is inside a centromere region
#' @return centromere.prox: proximity value to the closest centromere
#' @export
#' @examples #get centro/telo -mere related features for locus chr1:5000000
#' centroTelo_dist(chr=1,position=5000000)

centroTelo_dist <- function(chr,position){
  if(length(chr) == length(position)){
    #load db
    tmp = read.table('gaps.db',header=TRUE) # database can be downloaded at http://genome.ucsc.edu/cgi-bin/hgTables , choose All Tables -> table:gap
    tmp$chrom = substr(tmp$chrom,4,6)
    telo = tmp[which(tmp$type == 'telomere'),]
    centro = tmp[which(tmp$type == 'centromere'),]
    #telomere related
    bool.telo =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(telo$chrom == chr[i] & telo$chromStart <= position[i] & telo$chromEnd >= position[i]))))) 
    tmp = sapply(c(1:length(chr))[!bool.telo], function(i) min( min(abs(telo$chromStart[which(telo$chrom == chr[i])] - position[i])), min(abs(position[i] - telo$chromEnd[which(telo$chrom == chr[i])]))))
    dist.telo = rep(0.0001,length(position))
    dist.telo[!bool.telo] = 1/tmp 
    
    #centromere related
    bool.centro =  sapply(1:length(chr), function(i) as.integer(as.logical(length(which(centro$chrom == chr[i] & centro$chromStart <= position[i] & centro$chromEnd >= position[i]))))) 
    tmp = sapply(c(1:length(chr))[!bool.centro], function(i) min( min(abs(centro$chromStart[which(centro$chrom == chr[i])] - position[i])), min(abs(position[i] - centro$chromEnd[which(centro$chrom == chr[i])]))))
    dist.centro= rep(0.0001,length(position))
    dist.centro[!bool.centro] = 1/tmp 
    #output
    x = cbind(bool.telo,dist.telo,bool.centro,dist.centro)
    colnames(x) = c('telomere.bool','telomere.prox','centromere.bool','centromere.prox')
    return(x)
  }else{
    warning('Chr and Position vectors do not have the same length! \n')
  }
}

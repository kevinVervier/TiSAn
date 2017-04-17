#------------------------------------------------------#
# Main function for tissue-specific feature extraction #
#------------------------------------------------------#

#' Compute the feature representation for one or more genomic positions
#'
#' \code{feature_extraction}
#' @param chr chromosome number(s)
#' @param position genomic position(s)
#' @return a named matrix for feature space representation (~ 360 features), each row corresponds to one input position
#' @export
#' @examples #get features for locus chr1:5000000
#' feature_extraction(chr=c(1,2),position=c(5000000,10000000))
#'

#inputs are chromosome number and positions sequence
feature_extraction <- function(chr=1,position=5000000){
  #n-grams
  db = get_ngrams_count(chr,position)
  #GTEx-eQTLs
  x = GTEx_eQTL(chr,position)
  db = cbind(db,x)
  #RME-methyl
  x = RME_methyl(chr,position)
  db = cbind(db,x)
  # distance to centromere and telomere
  x = centroTelo_dist(chr,position)
  db = cbind(db,x)
  #distance to the closest gene
  x = getClosestGene(chr,position)
  db = cbind(db,x)
  # ... here users are invited to add any additional feature, especially ones coming from tissue-specific genome-wide studies
  # For instance, in online version of TiSAn heart model, we used developmental fetal heart enchancers from Heart Enhancer Compendium
  # Here is an example of the 'fetal_cardiac_enhancer' script call
  x = dvpt_enhancer_dist(chr,position)
  db = cbind(db,x)

  #output
  return(db)
}

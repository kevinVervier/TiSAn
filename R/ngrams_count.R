#-----------------------------------#
# N-nucleotide frequency extraction #
#-----------------------------------#

#' Compute the frequency for each 1,2,3,4-nucleotide in a +/-500bp neighbourhood around a locus
#' 
#' \code{ngrams}
#' @param chr the chromosome number (including 'chr' prefix)
#' @param position the locus position, center of 1kbp window
#' @return a named vector of length 340 (4 + 16 + 64 + 256) with ngram frequencies
#' @export
#' @examples #get ngrams count in chr1-4999500:5000500 interval
#' get_ngrams_count(chr=1,position=5000000)


#take as input a chromosome number and positions sequence
get_ngrams_count <- function(chr='chr1',position=5000000){
  s = BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr]] 
  ### get views on the 100kb sequence for 1000bp sliding window
  v = Views(s,IRanges(start = pmax(1,position - 500),end = position + 500))
  
  ### compute the matrix of [0,1] k-mer frequencies
  fq = cbind(
    oligonucleotideFrequency(v,1,as.prob = TRUE),
    oligonucleotideFrequency(v,2,as.prob = TRUE),
    oligonucleotideFrequency(v,3,as.prob = TRUE),
    oligonucleotideFrequency(v,4,as.prob = TRUE)
  )
  return(fq)
}

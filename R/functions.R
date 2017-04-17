##############################
# Literature mining function # 
##############################

#' Get access to pubmed ID given query on keywords in title/abstract/etc
#' 
#' \code{getIDs}
#' @param term a search term (character vector of length 1)
#' @param position the locus position
#' @return out: a character vector of PMIDs
#' @export
#' @examples #get PubMed IDs for a query
#' getIDs("heart[Title]")


getIDs = function(term){
  ## clean up input
  term = gsub(" ","+",term,ignore.case = FALSE)
  term = gsub("/",".",term)
  term = gsub('"',"%22",term)
  term = gsub("'","%27",term)
  term = gsub("-","%2D",term)
  
  ## construct the URL query
  base.url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&usehistroy=y"
  url = paste(base.url,"&term=",term,"&retmax=1",sep="")
  
  ## determine the total number of papers
  tmp <- RCurl::getURL(url=url)
  doc = XML::xmlTreeParse(tmp, useInternal=TRUE)
  root = XML::xmlRoot(doc)
  retmax = as.numeric(XML::xmlValue(root[[1]]))
  
  ## query again, this time setting retmax to the total
  ## number of articles
  url = paste(base.url,"&term=",term,"&retmax=",retmax,sep="")
  tmp <- RCurl::getURL(url=url)
  doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
  
  ## use XPath to extract the IDs from the XML document
  a = XML::xpathApply(doc,"//IdList/Id",XML::xmlValue)
  out = unlist(a)
  
  return(out)
}

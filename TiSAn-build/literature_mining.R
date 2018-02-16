
#######################################################
# Pubmed eutils can block too large number of queries

repeat_until_it_works <- function(catch, path, query, max_tries = 3, wait_time = 10, 
                                  messages = TRUE, key = NULL, ...) {
  error_handler <- function(e) {
    if (e$message %in% catch) {
      if (messages) warning(paste("Caught error:", e$message))
      return(NA)
    } else {
      stop(e$message)
    }
  }
  for (count in 1:max_tries) {
    tmp <- getURL(url=url)
    doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
    
    output <- tryCatch(res$parse("UTF-8"), error = error_handler)
    if (!is.na(output)) return(output)
    Sys.sleep(wait_time * count)
  }
  return(output)
}

### this is a function that returns a character vector of PMIDs, given a 
### search term
### term=character vector of length 1
### output value: character vector of PMIDs
getIDs = function(term,idset){
  if(!require(XML)){install.packages('XML')}
  if(!require(RCurl)){install.packages('RCurl')}
  
  
  ## clean up input
  term = gsub(" ","+",term,ignore.case = FALSE)
  term = gsub("/",".",term)
  term = gsub('"',"%22",term)
  term = gsub("'","%27",term)
  term = gsub("-","%2D",term)
  
  if(length(idset) > 700){
    
    # loop over all potential papers (limited by number of characters in a query)
    out = NULL
    for(i in 1:(length(idset)%/%700)){
      # store the list of IDs online
      base.url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=",paste(idset[(1+(i-1)*700):(700*i)],collapse=','),sep='')
      tmp <- getURL(url=base.url)
      doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
      webenv = unlist(XML::xpathApply(doc,"//WebEnv",XML::xmlValue))
      query_key = unlist(XML::xpathApply(doc,"//QueryKey",XML::xmlValue))
      
      ## construct the URL query
      base.url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&query_key="
      base.url = paste(base.url,query_key,"&WebEnv=",webenv,sep='') 
      
      # url = paste(base.url,"&term=",term,"&id=",paste(idset[(1+(i-1)*100):(100*i)],collapse=','),sep="")
      url = paste(base.url,"&term=",term,sep="")
      #errors_to_catch <- c("Could not resolve host: eutils.ncbi.nlm.nih.gov")
      #xml_result <- repeat_until_it_works(errors_to_catch,"esearch", query = list(db = "taxonomy", term = term, api_key = key), ...)
      
     tmp <- getURL(url=url)
      doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
      
      ## use XPath to extract the IDs from the XML document
      a = XML::xpathApply(doc,"//IdList/Id",XML::xmlValue)
      Sys.sleep(1)
      #out = unlist(a)
      out = c(out,unlist(a))
      if(i%%100 == 0){
        cat(i,'/',(length(idset)%/%100),'\n')
        Sys.sleep(10)
      }
    }
    
    if(700*i != length(idset)){
      #final loop
      # url = paste(base.url,"&term=",term,"&id=",paste(idset[(1+i*100):length(idset)],collapse=','),sep="")
      # store the list of IDs online
      base.url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=",paste(idset[(1+i*500):length(idset)],collapse=','),sep='')
      tmp <- getURL(url=base.url)
      doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
      webenv = unlist(XML::xpathApply(doc,"//WebEnv",XML::xmlValue))
      query_key = unlist(XML::xpathApply(doc,"//QueryKey",XML::xmlValue))
      
      ## construct the URL query
      base.url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&query_key="
      base.url = paste(base.url,query_key,"&WebEnv=",webenv,sep='') 
      tmp <- getURL(url=url)
      doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
      
      # use XPath to extract the IDs from the XML document
      a = XML::xpathApply(doc,"//IdList/Id",XML::xmlValue)
      out = c(out,unlist(a))
    }
  }else{
    # store the list of IDs online
    base.url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=pubmed&id=",paste(idset,collapse=','),sep='')
    tmp <- getURL(url=base.url)
    doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
    webenv = unlist(XML::xpathApply(doc,"//WebEnv",XML::xmlValue))
    query_key = unlist(XML::xpathApply(doc,"//QueryKey",XML::xmlValue))
    
    ## construct the URL query
    base.url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&query_key="
    base.url = paste(base.url,query_key,"&WebEnv=",webenv,sep='') 
    
    # url = paste(base.url,"&term=",term,"&id=",paste(idset[(1+(i-1)*100):(100*i)],collapse=','),sep="")
    url = paste(base.url,"&term=",term,sep="")
    tmp <- getURL(url=url)
    doc = XML::xmlTreeParse(tmp,useInternalNodes=T)
    
    ## use XPath to extract the IDs from the XML document
    a = XML::xpathApply(doc,"//IdList/Id",XML::xmlValue)
    #out = unlist(a)
    out = unlist(a)
  }
  return(out)
}

# example
#id = getIDs(term = 'brain[title]',idset = pubmed$V2)
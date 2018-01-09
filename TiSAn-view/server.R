#--------------------#
# TiSAn-view  server #
#--------------------#
if(!require(shiny)) install.packages('shiny')
if(!require(DT)) install.packages('DT')
if(!require(genomation)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("genomation")
}
if(!require(GenomicRanges)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

#light TisAn-brain version
load('rounded_brain_score.Rdata')
# DEBUG --> do we want to provide a check box for either heart or brain scores, and not load the whole data in.
load('rounded_heart_score.Rdata')
################ TODO ###################

# genes data (use Granges to be more efficient), with a tag 'brain' as metadata for easy filtering
load('brain_gene_gr.Rdata')
load('non_brain_gene_gr.Rdata')

# GTex eQTL  (granges + brain tag in metadata) --> rsid and tissue ? as metadata
#load('brain_eqtl_gr.Rdata')
load('brain_eqtl_gr.Rdata')
#load('non_brain_eqtl_gr.Rdata')
load('non_brain_eqtl_gr.Rdata')
# RME data (granges) with 2 metadata tags (brain and non-brain values)
load('rme_gr.Rdata')
# fetal ddmR
load('ddmr_gr.Rdata')

# logo

# if a region is provided, return the highest TiSAn location in the table

# plot a Manhattan like figure in case of a region being provided in bed file

shinyServer(function(input, output) {
  
  #read data function
  read_dataset <- reactive({
    # support bed file
    readBed(input$your_data$datapath,zero.based = FALSE)
  })
  
  # selected data set
  output$text1 <- renderText({ 
    paste("You have selected ", input$your_data$datapath,".",sep='')
  })
  # print data dimensions
  output$text_dim <- renderText({ 
    data = read_dataset()
    paste(length(data)[1], "loci.")
  })
#   # get current selected row in brain table
#   output$select_brain_loc <- renderText({
#    paste(input$brain_map_rows_selected[1])
#   })
  
  # load and plot the data
#   output$contents <- renderDataTable({
#     inFile <- input$your_data
#     if(is.null(inFile)) return(NULL)
#     data = read_dataset()#read.csv(inFile$datapath,row.names = 1)
#     # data = cbind(Row.names = rownames(data), data)
#     
#    return(as.data.frame(data))
#     
#     #datatable(as.data.frame(data), selection='single')
#   },options = list(pageLength = 10))

  # simply annotate with brain score and print table
  output$brain_map <- DT::renderDataTable({ 
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "UCSC")
    data <- renameSeqlevels(data, newStyle)
    
    tmp = findOverlaps(data,brain.gr)
    data@elementMetadata$brain.score = 0
    
    if(length(tmp) > 0){
    
    # derive a version in the case a region is provided
   # data@elementMetadata$brain.score[tmp@from] = brain.gr@elementMetadata$score[tmp@to]
      data@elementMetadata$brain.score[unique(tmp@from)] = tapply(INDEX = tmp@from,X = tmp@to,function(i)max(brain.gr@elementMetadata$score[i]))
      # get the corresponding bins to report in the output table
      data@elementMetadata$bins = paste(start(data),end(data),sep='-')
      data@elementMetadata$bins[unique(tmp@from)] = tapply(INDEX = tmp@from,X = tmp@to,function(i) paste(start(brain.gr[i[which.max(brain.gr@elementMetadata$score[i])]]),end(brain.gr[i[which.max(brain.gr@elementMetadata$score[i])]]),sep='-'))
      
    }
    data = as.data.frame(data)
    if(ncol(data) == 7){
      return(data[,c(1,2,3,6,7)])
    }else{
      return(data[,c(1,2,3,6)])
    }
  })
  
  # simply annotate with heart score and print table
  output$heart_map <- DT::renderDataTable({ 
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "UCSC")
    data <- renameSeqlevels(data, newStyle)
    
    tmp = findOverlaps(data,heart.gr)
    data@elementMetadata$brain.score = 0
    
    if(length(tmp) > 0){
      
      # derive a version in the case a region is provided
      # data@elementMetadata$brain.score[tmp@from] = brain.gr@elementMetadata$score[tmp@to]
      data@elementMetadata$brain.score[unique(tmp@from)] = tapply(INDEX = tmp@from,X = tmp@to,function(i)max(heart.gr@elementMetadata$score[i]))
      # get the corresponding bins to report in the output table
      data@elementMetadata$bins = paste(start(data),end(data),sep='-')
      data@elementMetadata$bins[unique(tmp@from)] = tapply(INDEX = tmp@from,X = tmp@to,function(i) paste(start(heart.gr[i[which.max(heart.gr@elementMetadata$score[i])]]),
                                                                                                                            end(heart.gr[i[which.max(heart.gr@elementMetadata$score[i])]]),sep='-'))
      
    }
    data = as.data.frame(data)
    if(ncol(data) == 7){
      return(data[,c(1,2,3,6,7)])
    }else{
      return(data[,c(1,2,3,6)])
    }
  })
  
  # brain eQTL
  output$brain_eqtl <- renderDataTable({
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    #data = data[input$brain_map_rows_selected[1]]
    data = data[input$brain_map_row_last_clicked]
    #newStyle <- mapSeqlevels(seqlevels(data), "UCSC")
    #data <- renameSeqlevels(data, newStyle)
    tmp = distanceToNearest(data,brain.eqtl)
    tmp = cbind(as.data.frame(brain.eqtl[tmp@to])[,c(1:2,6:7)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','pos','rsid','symbol','distance')
    row.names(tmp)=NULL
    return(tmp)
    
  })
  # non-brain eQTL
  output$nonbrain_eqtl <- renderDataTable({
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    data = data[input$brain_map_row_last_clicked]
    tmp = distanceToNearest(data,other.eqtl)
    tmp = cbind(as.data.frame(other.eqtl[tmp@to])[,c(1:2,6:7)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','pos','rsid','symbol','distance')
    row.names(tmp)=NULL
    return(tmp)
    
  })
  # brain gene
  output$brain_gene <- renderDataTable({
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    data = data[input$brain_map_row_last_clicked]
    tmp = distanceToNearest(data,brain.gene)
    tmp = cbind(as.data.frame(brain.gene[tmp@to])[,c(1:3,6)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','symbol','distance')
    row.names(tmp)=NULL
    return(tmp)
    
  })
  # non-brain gene
  output$nonbrain_gene <- renderDataTable({
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    data = data[input$brain_map_row_last_clicked]
    tmp = distanceToNearest(data,other.gene)
    tmp = cbind(as.data.frame(other.gene[tmp@to])[,c(1:3,6)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','symbol','distance')
    row.names(tmp)=NULL
    return(tmp)
  })
  # RME
  output$rme <- renderDataTable({
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    data = data[input$brain_map_row_last_clicked]
    tmp = distanceToNearest(data,rme)
    tmp = cbind(as.data.frame(rme[tmp@to])[,c(1:3,6:7)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','brain.methyl','non.brain.methyl','distance')
    row.names(tmp)=NULL
    return(tmp)
  })
  # ddmr
  output$ddmr <- renderDataTable({
    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    data = data[input$brain_map_row_last_clicked]
    tmp = distanceToNearest(data,ddmr)
    tmp = cbind(as.data.frame(ddmr[tmp@to])[,c(1:3)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','distance')
    row.names(tmp)=NULL
    return(tmp)
  })
  
  # renderPlotly() also understands ggplot2 objects!
  output$plot <- renderPlotly({

    data = read_dataset()
    mcols(data) = NULL
    newStyle <- mapSeqlevels(seqlevels(data), "UCSC")
    data <- renameSeqlevels(data, newStyle)
    
    tmp = findOverlaps(data,brain.gr)
    data@elementMetadata$brain.score = 0
    
    if(length(tmp) > 0){
      
      # derive a version in the case a region is provided
      # data@elementMetadata$brain.score[tmp@from] = brain.gr@elementMetadata$score[tmp@to]
      data@elementMetadata$brain.score[unique(tmp@from)] = tapply(INDEX = tmp@from,X = tmp@to,function(i)max(brain.gr@elementMetadata$score[i]))
      # get the corresponding bins to report in the output table
      data@elementMetadata$bins = paste(start(data),end(data),sep='-')
      data@elementMetadata$bins[unique(tmp@from)] = tapply(INDEX = tmp@from,X = tmp@to,function(i) paste(start(brain.gr[i[which.max(brain.gr@elementMetadata$score[i])]]),end(brain.gr[i[which.max(brain.gr@elementMetadata$score[i])]]),sep='-'))
      
    }
    data = sort(data)
    # read chr_size
    chr_size = matrix(c('chr1',	249250621,
                        'chr2',	243199373	,
                        'chr3',	198022430,
                        'chr4',	191154276,
                        'chr5',	180915260,
                        'chr6',	171115067,
                        'chr7',159138663,
                        'chr8',146364022,
                        'chr9',141213431,
                        'chr10',135534747,
                        'chr11',135006516,
                        'chr12',133851895,
                        'chr13',115169878,
                        'chr14',107349540,
                        'chr15',102531392,
                        'chr16',90354753,
                        'chr17',81195210,
                        'chr18',78077248,
                        'chr19',59128983,
                        'chr20',63025520,
                        'chr21',48129895,
                        'chr22',51304566,
                        'chrX',155270560,
                        'chrY',59373566),ncol = 2,byrow = TRUE)
    
    d = as.data.frame(data)
    d$chr = match(d$seqnames,chr_size[,1])
    d$pos=NA
    ticks=NULL
    lastbase=0
    #colors <- rep(colors,max(d$chr))[1:max(d$chr)]
    
    for (i in 1:nrow(chr_size)) {
      if (i==1) {
        d[d$chr==i, ]$pos=d[d$chr==i, ]$start
      } else {
        lastbase=lastbase+as.numeric(chr_size[i-1,2])#tail(subset(d,chr==i-1)$start, 1)
        d[d$chr==i, ]$pos=d[d$chr==i, ]$start+lastbase
      }
      ticks=c(ticks, lastbase+as.numeric(chr_size[i,2])/2)#d[d$chr==i, ]$pos[floor(length(d[d$chr==i, ]$pos)/2)+1])
    }
    
    a <- list(
      title = 'Chromosome',
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 0.25,
      ticklen = 5,
      tickwidth = 2,
      tickcolor = toRGB("blue"),
      tickmode = 'array',
      tickvals = ticks,
      ticktext = chr_size[,1]#unique(d$seqnames)
    )
    
    b <- list(
      title = 'TiSAn score',
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 0.25,
      ticklen = 5,
      tickwidth = 2,
      tickcolor = toRGB("blue"),
      range = c(0,1.2)
    )
    
    plot_ly(d, x=~pos, y=~brain.score , color = ~brain.score, colors = 'Purples',size = ~brain.score, type = 'scatter', mode = 'markers',source="source",
            text = ~paste('pos: ', seqnames,':',start,sep=''))%>%
      layout(xaxis = a, yaxis = b)
    #axis(1, at=ticks, lab=unique(d$seqnames))
   # icol=1
  #  for (i in unique(d$chr)) {
  #    with(d[d$chr==i, ],points(pos, brain.score, col=colors[icol],pch=16))
  #    icol=icol+1
  #  }
    
  })
  
  output$event <- renderPrint({
    
    out.list = list()
    eventdata <- event_data("plotly_click", source = "source")
    validate(need(!is.null(eventdata), "Click on a specific locus to get annotations"))
    datapoint <- as.numeric(eventdata$pointNumber)[1]
    
    data = read_dataset()
    mcols(data) = NULL
    
    newStyle <- mapSeqlevels(seqlevels(data), "NCBI")
    data <- renameSeqlevels(data, newStyle)
    data = sort(data)[datapoint+1]
    #data = data[datapoint+1]

    # brain gene
    tmp = distanceToNearest(data,brain.gene)
    tmp = cbind(as.data.frame(brain.gene[tmp@to])[,c(1:3,6)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','symbol','distance')
    row.names(tmp)=NULL
    out.list[['brain.gene']] = tmp
    # non brain gene
    tmp = distanceToNearest(data,other.gene)
    tmp = cbind(as.data.frame(other.gene[tmp@to])[,c(1:3,6)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','symbol','distance')
    row.names(tmp)=NULL
    out.list[['other.gene']] = tmp
    # brain eQTL
    tmp = distanceToNearest(data,brain.eqtl)
    tmp = cbind(as.data.frame(brain.eqtl[tmp@to])[,c(1:2,6:7)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','pos','rsid','symbol','distance')
    row.names(tmp)=NULL
    out.list[['brain eQTL']] = tmp
    # non brain eqtl
    tmp = distanceToNearest(data,other.eqtl)
    tmp = cbind(as.data.frame(other.eqtl[tmp@to])[,c(1:2,6:7)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','pos','rsid','symbol','distance')
    row.names(tmp)=NULL
    out.list[['other eQTL']] = tmp
    
    # RME
    tmp = distanceToNearest(data,rme)
    tmp = cbind(as.data.frame(rme[tmp@to])[,c(1:3,6:7)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','brain.methyl','non.brain.methyl','distance')
    row.names(tmp)=NULL
    out.list[['RME']] = tmp
    #ddmr
    tmp = distanceToNearest(data,ddmr)
    tmp = cbind(as.data.frame(ddmr[tmp@to])[,c(1:3)],tmp@elementMetadata$distance)
    colnames(tmp) = c('chr','start','end','distance')
    row.names(tmp)=NULL
    out.list[['DDMR']] = tmp

    return(out.list)
    
  })
  
})

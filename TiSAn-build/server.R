#----------------------------#
# server part of TiSAn-build #
#----------------------------#

# install required packages if needed
if(!require(shiny)) install.packages('shiny')
if(!require(shinyFiles)) install.packages('shinyFiles')
if(!require(DT)) install.packages('DT')
if(!require(biomaRt)) install.packages('biomaRt')
if(!require(genomation)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("genomation")
}
if(!require(GenomicRanges)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}


# literature mining (in-house) functions
source('literature_mining.R')

shinyServer(function(input, output,session) {

  # GTEX eQTL database location
  shinyDirChoose(input, 'eqtl_directory', session=session, roots=c(wd='.'))
  path_gtex <- reactive({
    return(print(parseDirPath(roots=c(wd='.'), input$eqtl_directory)))
  })
  # get all tissue names in GTEX
  output$gtex_tissue <- renderUI({
    if(is.null(input$eqtl_directory))
      return()

    gtex.files <- list.files(path = path_gtex(), pattern = '*.gene_pairs.txt.gz')
    # extract tissue names
    gtex.files <- gsub(pattern = '.v7.*',replacement = '',x = gtex.files)
    checkboxGroupInput("eQTL_columns", "Choose columns",
                       choices  =  gtex.files,
                       selected =  gtex.files[grep(pattern = 'Brain',x = gtex.files)])
  })


  # split all GTEx tissues in two files (tissue and non-tissue) only after 'Done' button is hit
  observeEvent(input$gtex_button, {
    if(is.null(input$eqtl_directory))
      return()
    # check if the databases already exist
    if(file.exists('tissue_eqtl_gr.Rdata'))
      return('It seems that the database already exists.') # TODO: print this message

    gtex.files <- list.files(path = path_gtex(), pattern = '*.gene_pairs.txt.gz',full.names = TRUE)
    # extract tissue names
    tissue.files <- gtex.files[grep(pattern = paste(input$eQTL_columns,collapse="|") ,x = gtex.files)]
    tmp = NULL
    withProgress(message = 'eQTL', value = 0, {
      n <- length(gtex.files)
      for(f in tissue.files){
        tmp = rbind(tmp,read.table(gzfile(f), header=T)[,1:2])
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Processing ", f))
      }
      # filter duplicated variants
      tmp = tmp[-which(duplicated(tmp[,1])),]
      tmp2 = do.call('rbind',strsplit(as.character(tmp[,1]),split='_'))

      tissue.eqtl = GRanges(seqnames = Rle(tmp2[,1]),
                            ranges = IRanges(start=as.numeric(tmp2[,2]), end=as.numeric(tmp2[,2])),
                            gene = Rle(tmp[,2])
      )
      save(tissue.eqtl,file='tissue_eqtl_gr.Rdata')

      # check if the databases already exist
      if(file.exists('non_tissue_eqtl_gr.Rdata'))
        return()
      tmp = NULL
      for(f in gtex.files[-which(gtex.files %in% tissue.files)]){
        tmp = rbind(tmp,read.table(gzfile(f), header=T)[,1:2])
        incProgress(1/n, detail = paste("Processing ", f))
      }
      # filter duplicated variants
      tmp = tmp[-which(duplicated(tmp[,1])),]
      tmp2 = do.call('rbind',strsplit(as.character(tmp[,1]),split='_'))

      other.eqtl = GRanges(seqnames = Rle(tmp2[,1]),
                           ranges = IRanges(start=as.numeric(tmp2[,2]), end=as.numeric(tmp2[,2])),
                           gene = Rle(tmp[,2])
      )
      save(other.eqtl,file='non_tissue_eqtl_gr.Rdata')
    })
  })

  #-----#
  # RME #

  # get RME file
  shinyFileChoose(input, 'rme_file', session=session, roots=c(wd='.'))
  path_rme <- reactive({
    return(print(parseFilePaths(roots=c(wd='.'), input$rme_file)))
  })
  # get RME metadata file
  shinyFileChoose(input, 'rme_metadata_file', session=session, roots=c(wd='.'))
  path_rme_meta <- reactive({
    return(print(parseFilePaths(roots=c(wd='.'), input$rme_metadata_file)))
  })

  # get all tissue names in RME
  output$rme_tissue <- renderUI({
    if(is.null(input$rme_file) | is.null(input$rme_metadata_file))
      return()

    # onlt need the header to get sample names
    #  methyl <- read.delim(as.character(rme_file_selected$datapath),header=T,nrows = 5)
    methyl <- read.delim(as.character(path_rme()$datapath),header=T,nrows = 5)
    methyl = methyl[,-24] #E064 is  missing, not part of the study anymore ?
    samples = unlist(lapply(strsplit(colnames(methyl)[-c(1:3)],split='_'),function(y)y[3]))

    # metadata <- read.delim(as.character(rme_meta_file_selected$datapath),header=T,as.is=T)
    metadata <- read.delim(as.character(path_rme_meta()$datapath),header=T,as.is=T)
    # keep only tissues found in samples reported in the methylation file
    idx = match(samples,metadata$Epigenome.ID..EID.)
    tissues = cbind(samples,metadata$ANATOMY[idx])

    # extract tissue names
    checkboxGroupInput("rme_columns", "Choose columns",
                       choices  =  sort(unique(tissues[,2]),decreasing=FALSE),
                       selected =  "BRAIN")
  })


  # split all RME tissues in two files (tissue and non-tissue) only after 'Done' button is hit
  observeEvent(input$rme_button, {
    if(is.null(input$rme_file) | is.null(input$rme_metadata_file))
      return()
    # check if the database already exists
    if(file.exists('rme_gr.Rdata'))
      return('It seems that the database already exists.') # todo: print this warning

    tmp = NULL
    withProgress(message = 'Methylation', value = 0, {
      n=5
      incProgress(1/n, detail = paste("Processing ", "methylation level"))
      methyl <- read.delim(as.character(path_rme()$datapath),header=T)
      methyl = methyl[,-24] #E064 is  missing, not part of the study anymore ?
      samples = unlist(lapply(strsplit(colnames(methyl)[-c(1:3)],split='_'),function(y)y[3]))
      incProgress(1/n, detail = paste("Processing ", "metadata"))
      metadata <- read.delim(as.character(path_rme_meta()$datapath),header=T,as.is=T)

      idx = match(samples,metadata$Epigenome.ID..EID.)
      tissues = cbind(samples,metadata$ANATOMY[idx])

      idx.tissue = which(tissues[,2] %in% input$rme_columns) +3
      idx.other =   which(!(tissues[,2] %in% input$rme_columns)) +3

      incProgress(1/n, detail = paste("Splitting ", "data"))

      tissue.db = methyl[,c(1:3,idx.tissue)]
      other.db = methyl[,c(1:3,idx.other)]

      # get average value for tissue and non-tissue samples
      av.tissue = apply(tissue.db[,-c(1:3)],1,function(x) mean(x,na.rm=TRUE))
      av.other = apply(other.db[,-c(1:3)],1,function(x) mean(x,na.rm=TRUE))

      rme = GRanges(seqnames = Rle(methyl[,1]),
                    ranges = IRanges(start=as.numeric(methyl[,2]), end=as.numeric(methyl[,3])),
                    tissue.methyl = Rle(av.tissue),
                    other.methyl = Rle(av.other)
      )
      incProgress(1/n, detail = paste("Saving ", "data"))
      save(rme,file='rme_gr.Rdata')
    })
  })
  ####################
  # literature mining

  # get Pubmed file
  shinyFileChoose(input, 'pubmed_file', session=session, roots=c(wd='.'))


  path_pubmed <- reactive({
    return(print(parseFilePaths(roots=c(wd='.'), input$pubmed_file)))
  })


  # split all genes in two files (tissue and non-tissue) only after 'Done' button is hit
  observeEvent(input$pubmed_button, {
    if(is.null(input$pubmed_file))
      return()
    # check if the database already exists
    if(file.exists('tissue_gene_gr.Rdata'))
      return('It seems that the databases already exist.') # todo: print this warning

    withProgress(message = 'Literature mining', value = 0, {
      n <- 5
      incProgress(1/n, detail = paste("Reading ", "gene2pubmed"))
      pubmed <- read.delim(as.character(path_pubmed()$datapath),header=T)
      # keep human only
      pubmed = pubmed[pubmed[,1] == '9606',]
      incProgress(1/n, detail = paste("Querying ", "Pubmed"))
      id = getIDs(term = input$pubmed_keyword,idset = pubmed[,3])
      # get papers found in gene2pubmed
      idx = which(pubmed[,3] %in% id)
      # get 3 citations or more
      idx.support = which(pubmed[,2] %in% names(which(table(pubmed[idx,2]) > 2)))

      tissue.genes = unique(pubmed[idx.support,2])
      other.genes = unique(pubmed[-idx.support,2])

      incProgress(1/n, detail = "Get tissue genes location")
      # create a bioMart to find gene location based on Entrez IDs
      ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
      ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

      tissue.anno <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position','hgnc_symbol'),
                           filters='entrezgene',
                           values=tissue.genes,
                           mart=ensembl)
      # remove non chromosomal locations
      tissue.anno <- tissue.anno[which(tissue.anno$chromosome_name %in% c(1:22,'X','Y')),]


      incProgress(1/n, detail = "Get non-tissue genes location")
      other.anno <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position','hgnc_symbol'),
                          filters='entrezgene',
                          values=other.genes,
                          mart=ensembl)
      # remove non chromosomal locations
      other.anno <- other.anno[which(other.anno$chromosome_name %in% c(1:22,'X','Y')),]

      incProgress(1/n, detail = paste("Saving ", "data"))
      # tissue genes
      tissue.gene = GRanges(seqnames = Rle(tissue.anno$chromosome_name),
                            ranges = IRanges(start= tissue.anno$start_position, end=tissue.anno$end_position),
                            symbol = tissue.anno$hgnc_symbol
      )

      save(tissue.gene,file='tissue_gene_gr.Rdata')
      # non-tissue genes
      other.gene = GRanges(seqnames = Rle(other.anno$chromosome_name),
                           ranges = IRanges(start= other.anno$start_position, end=other.anno$end_position),
                           symbol = other.anno$hgnc_symbol
      )
      save(other.gene,file='other_gene_gr.Rdata')

    })

  })

  #########################
  # training set - lincSNP

  # lincSNP database location
  shinyDirChoose(input, 'lincsnp_directory', session=session, roots=c(wd='.'))
  path_lincsnp <- reactive({
    return(print(parseDirPath(roots=c(wd='.'), input$lincsnp_directory)))
  })
  # get all disease names in lincSNP (at least frequent ones)
  output$linc_disease <- renderUI({
    if(is.null(input$lincsnp_directory))
      return()
    withProgress(message = 'Loading LincSNP database', value = 0, {
      linc.files <- list.files(path = path_lincsnp(), pattern = '.txt',full.names = TRUE)
      # read the files and find disease names
      dis.names <- NULL
      n <- length(linc.files)

      for(f in linc.files){
        incProgress(1/n, detail = paste("Reading", f))

        tmp <- read.delim(f,header=FALSE, colClasses=c(rep("NULL",21), NA, "NULL","NULL"),as.is = TRUE)[,1]
        dis.names = c(dis.names,tmp)
      }
      #dis.names <- as.character(dis.names)
      # can change this filter to see more diseases
      dis.names <- names(which(table(dis.names) > quantile(table(dis.names),p=0.999)))

      checkboxGroupInput("dis_names", "Choose columns",
                         choices  =  dis.names,
                         selected =   dis.names[grep(pattern = paste(c('Suicide attempts in bipolar disorder','Alzheimer`s disease','Parkinson`s disease','Bipolar disorder','Depression--quantitative trait','Cognitive performance','Major depressive disorder',
                                                                       'Schizophrenia','Alcohol','Personality dimensions','Narcolepsy','Information processing speed','Smoking behavior','White matter hypersensitivity burden','Attention deficit hyperactivity disorder',
                                                                       'Nicotine dependence','Alcoholism','Intelligence','Coffee consumption','Panic disorder','Neuroblastoma','Cognition phenotypes','Asperger disorder','Brain lesion load','Caffeine consumption',
                                                                       'Autism','Drinking behavior','Hippocampal atrophy','Response to lithium treatment in bipolar disorder','Conduct disorder','Suicidal ideation','Hoarding','Anorexia nervosa','Response to antipsychotic treatment'),collapse="|") ,x = dis.names)])
    })
  })


  # split all lincSNP in two files (tissue and non-tissue) only after 'Done' button is hit
  observeEvent(input$lincsnp_button, {
    if(is.null(input$lincsnp_directory))
      return()
    # check if the databases already exist
    if(file.exists('train_lincsnp_pos.txt'))
      return('It seems that the database already exists.') # TODO: print this message

    linc.files <- list.files(path = path_lincsnp(), pattern = '.txt',full.names = TRUE)

    withProgress(message = 'Screening LincSNP database', value = 0, {
      linc.files <- list.files(path = path_lincsnp(), pattern = '.txt',full.names = TRUE)
      # read the files and find disease names
      pos <- NULL
      neg <- NULL

      n <- length(linc.files) + 1

      for(f in linc.files){
        incProgress(1/n, detail = paste("Reading", f))
        tmp <- read.delim(f,header=FALSE,as.is = TRUE)[,c(1,16,18,22)]
        idx = which(tmp$V22 %in% input$dis_names)

        pos = rbind(pos,tmp[idx,])
        neg = rbind(neg,tmp[-idx,])
      }
      # filter missing rsIDs
      pos <- pos[grep(pattern = 'rs',pos[,3]),]
      neg <- neg[grep(pattern = 'rs',neg[,3]),]
      # filter duplicated IDs
      pos <- pos[!duplicated(pos[,3]),]
      neg <- neg[!duplicated(neg[,3]),]
      incProgress(1/n, detail = paste("Writing training sets"))
      write.table(pos,file='train_lincsnp_pos.txt',quote = FALSE,row.names = FALSE,col.names = c("chr","pos","rsid","disease"),sep = '\t')
      write.table(neg,file='train_lincsnp_neg.txt',quote = FALSE,row.names = FALSE,col.names = c("chr","pos","rsid","disease"),sep = '\t')

    })

  })




})

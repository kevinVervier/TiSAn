#------------------------#
# UI part of TiSAn-build #
#------------------------#

# install required packages if needed
if(!require(plotly)) install.packages('plotly')
if(!require(shinyFiles)) install.packages('shinyFiles')

shinyUI(fluidPage(#theme = "bootstrap.css",
  
  titlePanel("TiSAn-build: train your Tissue Specific Annotation!"),
  
  navlistPanel(
    "Introduction",
    tabPanel("Before starting...",mainPanel(
      h1('Few words before starting:'),
      p("Welcome to the TiSAn-build application."),
      p("This application aims at providing an interactive way for extracting tissue-specific features, and train your own TiSAn model."),
      p("Default settings allow to train the TiSAn-brain model (publically available).")
    )),
    #eQTL panel
    tabPanel("expression Quantitative Trait Loci (eQTL) features",mainPanel(
      p('First, download from GTEx portal (http://gtexportal.org/home/datasets): GTEx_Analysis_v7_eQTL.tar.gz and untar it in the TISAn-build folder.'),
      p('Once it is done, provide its directory location:'),
      tags$p(shinyDirButton('eqtl_directory', 'Push to select a directory', 'Please select a folder'), align= "center"),#, value = 'GTEx_Analysis_v7_eQTL'),
      p('NB: '),
      tags$ul(
        tags$li('This folder contains one file for each tissue, with all the significant variant-gene associations.'),
        tags$li('It is technically possible to use a different database, only if the file format is similar to GTEx one.')
      ),
     # if(!is.null(input$eqtl_directory)){ 
        p('Once the file selection is completed, hit the button to create tissue and non-tissue databases'),
        actionButton("gtex_button", "eQTL tissue selection done"),
      #  },

      uiOutput("gtex_tissue")
    )),
    #Methylation panel
    tabPanel("RoadMap Epigenomics (RME) features",mainPanel(
      p('First, download http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/DMRs/WGBS_DMRs_v2.tsv.gz and untar it in the TISAn-build folder.'),
      p('Then, we also need the metadata file mapping sample ID to tissue (sheet 1 in http://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM).'),
    
      p('Once it is done, provide the location of the methylation file:'),
      tags$p(shinyFilesButton('rme_file', 'Push to select a file', 'Please select a file', multiple = FALSE), align= "center"),

      p('and also provide the location of the metadata file:'),
      tags$p(shinyFilesButton('rme_metadata_file', 'Push to select a file', 'Please select a file', multiple = FALSE), align= "center"),
      
      p('NB: '),
      tags$ul(
        tags$li('The methylation file contains one column per cell line with observed DNA methylation.'),
        tags$li('It is technically possible to use a different database, only if the file format is similar to RME one.')
      ),
      #if(!is.null(input$rme_file)){ 
        p('Once the file selection is completed, hit the button to create tissue and non-tissue databases'),
        actionButton("rme_button", "Methylation tissue selection done"),
      #},
      
      uiOutput("rme_tissue")
    )),
    #Literature mining for tissue-specific genes
    tabPanel("Literature Mining Genes features",mainPanel(
      p('First, download ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz and untar it in the TISAn-build folder.'),
      
      p('Once it is done, provide the location of the gene2pubmed file:'),
      tags$p(shinyFilesButton('pubmed_file', 'Push to select a file', 'Please select a file', multiple = FALSE), align= "center"),
      p('Then, provide a list of keywords related to the tissue of interest:'),
      tags$p(textInput("pubmed_keyword", "tissue-specific keywords", "brain[title]"), align= "center"),
      p('NB: '),
      tags$ul(
        tags$li('This file contains a column for gene ENTREZ ID and one column for PubMed ID.'),
        tags$li('Keywords will be queried in publications title.'),
        tags$li('Multiple keywords can be provided, as long as they are separated by "|" (example: "brain|neuron").'),
        tags$li('It is technically possible to use a different database, only if the file format is similar to gene2pubmed one.')
      ),
      p('Once the file selection is completed, hit the button to create tissue and non-tissue databases'),
      actionButton("pubmed_button", "tissue-specific keywords provided")
    )),
    #Custom database
    tabPanel("Custom database",mainPanel(
      p('This part is left for users to include more features with domain-specific data.'),
      
      p('One can adapt what is done in the previous three panels depending on the data type.')
    )),
    # Panel where training position are extracted
    tabPanel("Training set composition",mainPanel(
      p('In this section, user provides sets of both tissue and non-tissue related loci.'),
      p('We recommend those loci to be associated with diseases in both sets.'),
      
      #tags$h1("Non-coding loci"),

      p('To train our models, we relied on LinCSNP database (http://210.46.80.146/lincsnp/LncRNA-ldSNP.zip) and unzip it in the TISAn-build folder.'),
      p('Once it is done, provide its directory location:'),
      tags$p(shinyDirButton('lincsnp_directory', 'Push to select a directory', 'Please select a folder'), align= "center"),
      p('Once the file selection is completed, hit the button to create tissue and non-tissue databases'),
      actionButton("lincsnp_button", "lincSNP disease selection done"),
      uiOutput("linc_disease"),

      p('NB: '),
      tags$ul(
        tags$li('This directory contains a file for each chromosome containing disease-associated SNPs.'),
        tags$li('Multiple diseases can be provided.'),
        tags$li('It is technically possible to use a different database, only if the file format is similar to lincSNP one.')
      )
      
      
    ))
  )
))
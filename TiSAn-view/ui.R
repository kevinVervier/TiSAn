#----------------#
# TiSAn-view UI  #
#----------------#

if(!require(plotly)) install.packages('plotly')

shinyUI(fluidPage(#theme = "bootstrap.css",
  
  titlePanel("TiSAn: Tissue Specific Annotation for Genetic Variations"),
  
  navlistPanel(
    "Introduction",
    # warm-up
    tabPanel("Before starting...",mainPanel(
      h1('Few words before starting:'),
      p("Welcome to the TiSAn application."),
      p("This application aims at providing an interactive way for annotating genetic loci with tissue-specific features."),
      p("First, we propose you to load the loci of interest.")
    )),
    #choose data to work on
    tabPanel("Load data",mainPanel(
      p('Please select your working data (supported bed format):'),
      fileInput('your_data', label='Select your loci file', multiple = FALSE,accept = '*.bed'),
      p('NB: '),
      tags$ul(
        tags$li('current annotation databases are built on GrCh37/hg19 coordinates.'),
        tags$li('different output is provided depending of loci nature.'),
        tags$li('for single nucleotide loci: position-specific annotation.'),
        tags$li('for genomic regions: Mahnattan plot for TiSAn distribution, and highlight of top hits.'),
        tags$li('no chromosome X and Y annotation due to lack of available data.'),
        tags$li('To fasten the annotation process, TiSAn scores were rounded.'),
        tags$li('Link to the Complete Tisan Scores can be found in the Github repository.')
      ),
      textOutput("text1")
    )),
    # descriptive analysis
    tabPanel("Data viewer",mainPanel(
      #p('The current data you are using contains:'),
      #textOutput("text_dim"),
      #dataTableOutput('contents')
      DT::dataTableOutput("brain_map")
    )),
   # tags$style(HTML('table.dataTable tr.selected td, table.dataTable td.selected {background-color: white !important;}')),
    tabPanel("Brain annotation",mainPanel(
      #p('Apply TiSAn-Brain score to the loci of interest:'),
      #DT::dataTableOutput("brain_map"),
      # tags$head(tags$style("#brain_map table {background-color: white; }", media="screen", type="text/css")),
      # p('Current selected locus:'),
      # textOutput("select_brain_loc"),
#       checkboxInput("close_gene_brain", "closest brain gene"),
#       conditionalPanel(
#         condition = "input.close_gene_brain == true",
#         dataTableOutput("brain_gene")
#       ),
#       checkboxInput("close_eqtl_brain", "closest eQTL from GTEx brain tissues"),
#       conditionalPanel(
#         condition = "input.close_eqtl_brain == true",
#         dataTableOutput("brain_eqtl")
#       ),
#       checkboxInput("close_gene_nonbrain", "closest non-brain gene"),
#       conditionalPanel(
#         condition = "input.close_gene_nonbrain == true",
#         dataTableOutput("nonbrain_gene")
#       ),
#       checkboxInput("close_eqtl_nonbrain", "closest eQTL from GTEx non-brain tissues"),
#       conditionalPanel(
#         condition = "input.close_eqtl_nonbrain == true",
#         dataTableOutput("nonbrain_eqtl")
#       ),
#       checkboxInput("close_methyl", "closest methylated region (RME)"),
#       conditionalPanel(
#         condition = "input.close_methyl == true",
#         dataTableOutput("rme")
#       ),
#       checkboxInput("ddmr", "differentially methylated region in fetal brain"),
#       conditionalPanel(
#         condition = "input.ddmr == true",
#         dataTableOutput("ddmr")
#       ),
      plotlyOutput("plot"),
      verbatimTextOutput("event")
    )),
    tabPanel("Heart annotation",mainPanel(
      p('Apply TiSAn-Heart score to the loci of interest:'),
      p('You can also click on a given row to get details on how the score was estimated for the correspoding locus.'),
      p('It is possible to export the whole table.'),
      DT::dataTableOutput("heart_map")
    ))
  )
))


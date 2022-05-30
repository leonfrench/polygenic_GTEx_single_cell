library(readxl)
library(here)
library(shiny)
library(readr)
library(magrittr)
library(tibble)
library(tidyr)
library(pheatmap)
library(dplyr)
library(purrr)
#library(shinyjs)
source("./AUCFunction.R")
source("./string_processing.R")

#deploy with:
#rsconnect::deployApp(appFileManifest = "appFileManifest.txt", account = 'polygenic')


apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  #tags$head(includeHTML("google-analytics.html")),
  # App title ----
  titlePanel("Polygenic tester for the GTEx single cell snRNA-seq data from Eraslan et al."),
  
  # Sidebar panel for inputs ----
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = 'dataset',
        label = 'Resolution:',
        choices = c(
          'Granular (170)',
          'Broad (113)'
        )
      ),
      
      textAreaInput(
        inputId = "genelist",
        label = "Input your gene list:",
        value = 'IDH1\nLIPE\nACSL1\nADH1B\nPPARG\nGPAM\nADIPOQ\nEBF1\nACACB\n',
        rows = 10
      ),
      textAreaInput(
        inputId = "background_genelist",
        label = "Background gene list (optional):",
        value = '',
        rows = 2
      ),
      selectInput(
        inputId = 'species',
        label = 'Species of input genes (data file used appears to be human gene symbols):',
        choices = c('Human', 'Mouse', 'Monkey')
      ),
      actionButton(inputId = "submit",
                   label = "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv"),
      downloadButton(outputId = "download_heatmap", label = "Download heatmap"),
      hr(),
      
      
      tags$b("Thank you to Gökcen Eraslan for sharing the full differential expression matrices used here."),
      tags$b("Full data from the Atlas are described in and available at:"),
      br(),
      tags$a(href="https://www.science.org/doi/10.1126/science.abl4290", "Gökcen Eraslan et al. paper"),
      br(),
      tags$a(href="www.gtexportal.org", "GTEx Portal"),
      br(),
      tags$a(href="https://singlecell.broadinstitute.org/single_cell/study/SCP1479", "Broad Institute Single Cell Portal"),
      br(),
      hr(),
      tags$a(href="https://github.com/leonfrench/polygenic_GTEx_single_cell", "Source code"),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(
        id = "main",
        # Output: Verbatim text for data summary ----
        verbatimTextOutput("summary"),
        br(),
        dataTableOutput("view"),
        br(),
        verbatimTextOutput("info")
      )
    )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  
  
  output$summary <- renderPrint({
    cat("\nTissues used: heart, prostate, skin, lung, breast, esophagus muscularis, esophagus mucosa and skeletal muscle")
    cat("\nResults will load here when complete")
    cat("\n")
    #print(gc())
    #print(Sys.info()['nodename'])
  })
  
  observeEvent(input$submit, {
    if (input$dataset == 'Broad (113)') {
      cell_expression_ranks <- read_csv(here("data", "de_broad.pkl.slim.csv.zip.ranks.csv"), show_col_types = FALSE)
    } else {
      cell_expression_ranks <- read_csv(here("data", "de_granular.pkl.slim.csv.zip.ranks.csv"), show_col_types = FALSE)
    }
    
    
    start <- Sys.time()
    
    cleaned_gene_list <-
      isolate(process_input_genes(input$genelist))
    
    #gene_list <- 'IDH1\nLIPE\nACSL1\nADH1B\nPPARG\nGPAM\nADIPOQ\nEBF1\nACACB\n'
    #cleaned_gene_list <- process_input_genes(gene_list)
    #cleaned_gene_list <- convert2mac(input_genes = cleaned_gene_list, in_species = "Human")
    # load reference data
    
    cleaned_gene_list <- convert2human(input_genes = cleaned_gene_list, in_species = input$species)
    first_cleaned_gene_list <- cleaned_gene_list
    background_cleaned_gene_list <- isolate(process_input_genes(input$background_genelist))
    background_cleaned_gene_list <- convert2human(input_genes = background_cleaned_gene_list, in_species = input$species)
    
    if (length(background_cleaned_gene_list) < 2 || background_cleaned_gene_list == "") {
      background_cleaned_gene_list <- cell_expression_ranks$gene_symbol
      
    } else { #if given a background
      print(length(background_cleaned_gene_list))
      cell_expression_ranks %<>% filter(gene_symbol %in% background_cleaned_gene_list)
      print(dim(cell_expression_ranks))
    }
    cleaned_gene_list <- intersect(cleaned_gene_list, background_cleaned_gene_list)
    #re rank based on new background
    cell_expression_ranks %<>% mutate_if(is.numeric, rank)
    
    
    print(paste0("Before time taken:", Sys.time() - start))
    #for indices - use dplyr for ease
    forIndices <- as_tibble(cell_expression_ranks$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- cell_expression_ranks %>%
      select(-gene_symbol)
    
    AUROC <- map_df(df, auroc_analytic, targetIndices)
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = cell_type, value = AUROC), 
                       gather(wilcox_tests, value = pValue)) %>% select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    # these are the values for the results table
    table %<>% mutate(pValue = signif(pValue, digits = 3), 
                      AUROC = signif(AUROC, digits = 3),
                      adjusted_P = signif(p.adjust(pValue, method = "bonferroni"), digits = 3))
    
    table %<>% mutate(tissue = gsub(".*[|]", "", cell_type))
    table %<>% mutate(cell_type = gsub("[|].*", "", cell_type))
    table %<>% select(cell_type, tissue, everything())
    
    table %<>% mutate(cell_type = gsub("DC", "dendritic cell", cell_type))
    table %<>% mutate(cell_type = gsub("LAM", "lipid-associated macrophage", cell_type))
    #
    table %<>% mutate(cell_type = gsub("SMC", "smooth muscle cell", cell_type))
    
    table %<>% arrange(-AUROC)
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste(
        "\nGenes found in data:",
        length(cleaned_gene_list),
        "of",
        length(first_cleaned_gene_list)
      ))
      cat(paste(
        "\nBackground genes:", length(background_cleaned_gene_list)
      ))
    })
    
    output$view <- renderDataTable(    table, escape = FALSE, options = list(
      paging =FALSE,
      pageLength =  51 
    ))
    #code duplication with HPA 
    output$download_heatmap <-
      downloadHandler(
        filename = "polygenic_GTEx_single_heatmap.pdf",
        content = function(file) {
          df <- as.data.frame(cell_expression_ranks %>% filter(gene_symbol %in% cleaned_gene_list))
          rownames(df) <- df$gene_symbol
          df <- df[ , (names(df) != "gene_symbol"), drop = T]
          
          colnames(df) <- gsub("[|]", " ", colnames(df))
          
          plot_out <- pheatmap(df, main = paste0("GTEx single cell heatmap for ", length(cleaned_gene_list), " genes\ncolor scale represents specific expression rank (higher is more specific)"))
          
          pdf(file, width=(23/170) * ncol(df), height=7.5+nrow(df) *.2)
          grid::grid.newpage()
          grid::grid.draw(plot_out$gtable)
          dev.off()
        }
      )
    
    output$download_data <-
      downloadHandler(
        filename = "polygenic_GTEx_single_cell_results.csv",
        content = function(file) {
          write_csv(table, file)
        }
      )
  })
}

shinyApp(ui, server)
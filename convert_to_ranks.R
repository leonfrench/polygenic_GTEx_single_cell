library(here)
library(shiny)
library(ggplot2)
library(readr)
library(magrittr)
library(plotly)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(reshape2)

target_file <- here("data", "de_granular.pkl.slim.csv.zip")
#target_file <- here("data", "de_broad.pkl.slim.csv.zip")
GTEx_single_deg <- read_csv(target_file )
GTEx_single_deg %>% pull(gene) %>% unique() %>% length()
GTEx_single_deg %<>% rename(gene_symbol = gene) 


#create single col for description for convienence
GTEx_single_deg %<>% mutate(cell_by_tissue = paste0(celltype, "|", tissue))
GTEx_single_deg %<>% select(-celltype, -tissue)
GTEx_single_deg
number_of_cellbytissues <- GTEx_single_deg %>% pull(cell_by_tissue) %>% unique() %>% length()
genes_with_all_tests <- GTEx_single_deg %>% group_by(gene_symbol) %>% 
  select(-tstat) %>% distinct() %>% count() %>% filter(n==number_of_cellbytissues) %>% pull(gene_symbol) 

GTEx_single_deg %<>% filter(gene_symbol %in% genes_with_all_tests)
GTEx_single_deg %<>% group_by(cell_by_tissue) %>% mutate(rank = rank(tstat))

GTEx_single_deg %<>% select(-tstat) %>% pivot_wider(names_from= "cell_by_tissue", values_from = "rank") 

dim(GTEx_single_deg)
GTEx_single_deg %>% write_csv(paste0(target_file, ".ranks.csv"))



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

#target_file <- here("data", "de_granular.pkl.slim.csv.zip")
target_file <- here("data", "de_broad.pkl.slim.csv.zip")
GTEx_single_deg <- read_csv(target_file )
GTEx_single_deg %>% pull(gene) %>% unique() %>% length()
GTEx_single_deg %<>% rename(gene_symbol = gene) 
GTEx_single_deg

#create single col for description for convienence
GTEx_single_deg %<>% mutate(cell_by_tissue = paste0(celltype, "|", tissue))



GTEx_single_deg %<>% select(-celltype, -tissue, -pvals_fdr)
GTEx_single_deg

GTEx_single_deg %>% group_by(cell_by_tissue) %>% filter(pvals_fdr < 0.05) %>% count() %>% arrange(n)

number_of_cellbytissues <- GTEx_single_deg %>% pull(cell_by_tissue) %>% unique() %>% length()
genes_with_all_tests <- GTEx_single_deg %>% group_by(gene_symbol) %>% 
  select(-tstat) %>% distinct() %>% count() %>% filter(n==number_of_cellbytissues) %>% pull(gene_symbol) 

GTEx_single_deg %<>% filter(gene_symbol %in% genes_with_all_tests)
GTEx_single_deg %<>% group_by(cell_by_tissue) %>% mutate(rank = rank(tstat))

dir.create(here("results", "markers", "broad"), recursive = T)
GTEx_single_deg %>% group_by(cell_by_tissue) %>% arrange(rank) %>% top_n(n = 100) %>% arrange(cell_by_tissue) %>% 
  group_walk(~ write_csv(.x %>% select(gene_symbol) %>% distinct(), here("results", "markers", "broad", paste0(gsub("[/]", " or ", gsub("[|]", "_", .y$cell_by_tissue)), ".txt")), col_names = F, quote = "none"))

#write out markers
GTEx_single_deg %<>% select(-tstat) %>% pivot_wider(names_from= "cell_by_tissue", values_from = "rank") 

dim(GTEx_single_deg)
GTEx_single_deg %>% write_csv(paste0(target_file, ".ranks.csv"))




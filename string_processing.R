library(stringr)
library(homologene)

process_input_genes <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}

convert2human <- function(input_genes, in_species) {
  if (in_species == 'Mouse') {
    human_genes <- mouse2human(input_genes)
    return(unique(human_genes$humanGene))
  } else if (in_species == 'Rhesus Macaque') {
    human_genes <- homologene(input_genes, inTax = 9544, outTax = 9606)
    return(unique(human_genes$humanGene))
  } else {
    return(unique(input_genes))
  }
}
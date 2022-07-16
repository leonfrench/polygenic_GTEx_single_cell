# from https://github.com/sarbal/EGAD/blob/master/EGAD/R/auroc_analytic.R
# by Sara Ballouz
auroc_analytic <- function(scores, labels) {
  
  negatives <- which(labels == 0, arr.ind = TRUE)
  scores[negatives] <- 0
  
  p <- sum(scores, na.rm = TRUE)
  nL <- length(labels)
  np <- sum(labels, na.rm = TRUE)
  nn <- nL - np
  
  auroc <- (p/np - (np + 1)/2)/nn
  
  return(auroc)
} 

#for ranked_profiles, the columns are the samples or experiments and the rows are what is being marked
#assumes the indices are lined up
auroc_analytic_ranked_profiles <- function(ranked_profiles, targetIndices) {
  #check the length of targetIndices is the row length of ranked_profiles
  if (length(targetIndices ) != nrow(ranked_profiles)) {
    stop("Index length does not equal input table row count!")
  }
  summed_ranks <- targetIndices %*% as.matrix(ranked_profiles)
  
  nL = nrow(ranked_profiles)  #number of labels/ total length
  np <- sum(targetIndices)     #number of positives
  nn <- nL - np #number of negatives
  
  AUCs <- (summed_ranks/np - (np + 1)/2)/nn
  return(as_tibble(AUCs))
}

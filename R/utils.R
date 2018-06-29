#' Define the composition of all samples in a metatranscriptomics dataset
#'
#' \code{loadPasilla} imports the Pasilla dataset from the corresponding package.
#'
#' This is an internal function and should not be called directly by the user, only internally by the \code{simulateSingleGenome} function
#'
#' @import Pasilla
#' @value Matrix of count values for the untreated samples of the Pasilla dataset, each column representing one sample and each row one gene
#' @example
#' counts_pasilla <- loadPasilla()
## TODO: find another public dataset that is better suitable (i.e. prokaryotic?)
loadPasilla=function(){
  pasCts = system.file("extdata","pasilla_gene_counts.tsv",package="pasilla", mustWork=TRUE)
  countsPasilla = read.table(pasCts, header = T, stringsAsFactors=F, sep="\t")
  countsPasilla=countsPasilla[,grep(colnames(countsPasilla),patter="untreated")]
  countsPasilla=cbind(countsPasilla, apply(countsPasilla, MARGIN=1, FUN=function(X){length(which(X>0))}))
  countsPasilla=countsPasilla[countsPasilla[,ncol(countsPasilla)]>1,]
  countsPasilla=unique(countsPasilla)
  countsPasilla=countsPasilla[,grep(colnames(countsPasilla),pattern="untreated")]
  return(countsPasilla)
}



#' Define the composition of all samples in a metatranscriptomics dataset
#'
#' \code{empiricalComposition} returns a composition matrix with rows as species/genomes and columns as samples (cases or controls). This composition matrix follows an empirical distribution which has been fitted from real metatranscriptomics datasets (Metatrans paper).
#' This is an internal function and should not be called directly by the user, only internally by the \code{compositionGenomesMetaT} function
#'
#' @param empiricalSeed A single number or a numeric vector of length equal to nReplicates. Indicates the random seed to assign the reads to different species in each sample. If NULL or a single number, the same seed will be applied to all samples so that they will have the same composition, with differences only due to random sampling.
#' @param genomes Character vector of genome names or genome IDs for the genomes to include in the simulation.
#' @param nReplicates A number, indicating the total samples (cases and controls).
#' @return A microbial composition matrix of \code{nReplicates} columns and \code{nrow(genomes)} rows.
#' @examples
#' # define a list of genomes to simulate
#' genomesToSimulate <- c("F. prausnitzii", "R. intestialis", "L. lactis", "E. faecalis", "R. obeum")
#' # obtain the empirical composition matrix for this 5 species
#' empiricalComposition(empiricalSeed=1, genomes=genomesToSimulate, nReplicates=10)
#' @seealso compositionGenomesMetaT
empiricalComposition <- function(empiricalSeed = NULL, genomes, nReplicates){
  compositionMatrix <- data.frame(row.names = genomes)
  if (is.null(empiricalSeed)){
    empiricalSeed <- rep(1, times = nReplicates)
  } else if(length(empiricalSeed) == 1){
    empiricalSeed <- rep(empiricalSeed, times = nReplicates)
  } else if(length(empiricalSeed) != nReplicates){
    print("ERROR! Please provide an empirical seed of length equal to the number of replicates")
    break
  }

  # 1. load data -----------------------------------------------------
  load("data/speciesDist.RData")

  # 2. remove those species present in less than one sample ----------
  zerovector<- c()
  for (i in 1:nrow(ssmatrix_real)){
    total_zeros <- sum(ssmatrix_real[i, ] == 0)
    zerovector <- c(zerovector, total_zeros)
  }
  to_remove <- which(zerovector == ncol(ssmatrix_real) | zerovector == ncol(ssmatrix_real) -1)
  ssmatrix_real <- ssmatrix_real[-to_remove, ]

  # 3. check the parameters of the negative binomial distribution ----
  r_vector <- c()
  mu_vector <- c()
  p_vector <- c()

  for (i in 1:nrow(ssmatrix_real)){
    fitted <- fitdist(as.vector(c(unlist(ssmatrix_real[i,]))), "nbinom")
    r <- as.vector(fitted[[1]][1])
    r_vector <- c(r_vector, r)
    mu <- as.vector(fitted[[1]][2])
    mu_vector <- c(mu_vector, mu)
    p <- r/(r + mu)
    p_vector <- c(p_vector, p)
  }

  # 4. write the composition matrix ----------------------------------
  for (i in 1:length(empiricalSeed)){
    set.seed(empiricalSeed[i])
    for (j in 1:length(genomes)){
      index = sample(1:length(p_vector),1)
      compositionMatrix[j,i] <- rnbinom(1, size = r_vector[index], prob = p_vector[index])
    }
  }
  compositionMatrix[is.na(compositionMatrix)] <- 0
  return(compositionMatrix)
}




#' Define the species composition of all samples in a metatranscriptomics dataset
#'
#' \code{compositionGenomesMetaT} returns a composition matrix with rows as species/
#'    genomes and columns as samples (cases or controls)
#'
#' @param composition A character string indicating which method to use to determine
#'    the composition of the matrix. Any of "custom", "empirical" or "even" can be used.
#'    See "details"
#' @param empiricalSeed A single number or a numeric vector of length equal to
#'    nReplicates. Indicates the random seed to assign the reads to different species
#'    in each sample. If NULL or a single number, the same seed will be applied to all
#'    samples so that they will have the same composition, with differences only due to
#'    random sampling
#' @param genomes Character vector of genome names or genome IDs for the genomes to
#'    include in the simulation
#' @param compositionMatrix A composition matrix predefined by the user, in which rows
#'    represent species and columns represent samples. Id user provides a composition
#'    matrix, this function only performs the random sampling to scale to the number of
#'    reads desired. Defaults to NULL
#' @param nReads An integer, indicating the number of reads to simulate per sample.
#' @param nReplicates An integer, indicating the total samples (cases and controls).
#' @param seed An integer, sets the random seed for the read distribution.
#' @details There are three options to distribute the total reads among the different
#'    species to simulate:
#'    - custom: requires that the user provides a custom composition matrix and scales it
#'    via random sampling (with replacement) to the specified number of reads in
#'    \code{nReads}
#'    - even: reads are distributed evenly among all the species listed in the
#'    \code{genomes} argument.
#'    - empirical: fits a negative binomial distribution (Vandeputte D. et al, 2017) to
#'    a real dataset (Martinez X. et al, 2016) of metatranscriptomics data. Uses this
#'    distribution to generate a composition matrix from the species listed in the
#'    \code{genomes} argument.
#' @references
#'   Martinez X. et al (2016): MetaTrans: an open source pipeline for metatranscriptomics
#'   Scientific Reports 6, Article number: 26447
#'   Vandeputte D. et al (2017): Quantitative microbiome profiling links gut community
#'   variation to microbial load. Nature 551:507–511
#' @import stats
#' @import utils
#' @export
#' @return A microbial composition matrix of \code{nReplicates} columns and
#'     \code{nrow(genomes)} rows.
#' @examples
#' # define a list of genomes to simulate
#' genomesToSimulate <- c("F. prausnitzii", "R. intestialis", "L. lactis", "E. faecalis",
#'                        "R. obeum")
#' # obtain the empirical composition matrix for this 5 species
#' compositionGenomesMetaT(composition="empirical", empiricalSeed=1,
#'                         genomes=genomesToSimulate, nReads=500000,nReplicates=10)
# TODO: important - make nReads possibly a vector, so that experiment can have a different sequencing depth
compositionGenomesMetaT <- function(composition=c("custom", "empirical", "even"), empiricalSeed=NULL, genomes, compositionMatrix=NULL, nReads=1000000, nReplicates=10, seed=42){
  genomeReadMatrix <- data.frame(row.names = genomes)

  # generate or check that the user has provided a composition matrix ----
  if(composition == "custom" & is.null(compositionMatrix)){
    print("Please provide a composition matrix, with rows for each genome and columns for each sample")
  } else if(composition == "empirical"){
    if(is.null(empiricalSeed)){
      empiricalSeed <- rep(1, times=nReplicates)
    } else if(length(empiricalSeed) == 1){
      empiricalSeed <- rep(empiricalSeed, times = nReplicates)
    } else if(length(empiricalSeed) != nReplicates){
      print("ERROR! Please provide an empirical seed of length equal to the number of replicates")
      break
    }
    compositionMatrix <- suppressWarnings(empiricalComposition(empiricalSeed, genomes, nReplicates))
  } else if(composition == "even"){
    compositionMatrix <- as.data.frame(matrix(1, nrow=length(genomes), ncol=nReplicates))
    rownames(compositionMatrix) <- genomes
  }

  # distribute reads per genome ----
  set.seed(seed)
  for(col in 1:ncol(compositionMatrix)){
    temp <- sample(rownames(compositionMatrix), size = nReads ,prob=compositionMatrix[,col], replace = T)
    genomeReadMatrix <- cbind(genomeReadMatrix,table(temp)[rownames(genomeReadMatrix)])
  }
  genomeReadMatrix <- genomeReadMatrix[,c(F,T)]
  genomeReadMatrix[is.na(genomeReadMatrix)] <- 0
  colnames(genomeReadMatrix) <- c(paste("sample", 1:ncol(genomeReadMatrix), sep="_"))
  return(genomeReadMatrix)
}





#' Calculate reads for one genome, for all samples
#'
#' \code{compositionGenomesMetaT} returns a composition matrix with rows as species/
#'    genomes and columns as samples (cases or controls)
#'
#' @param genomename A character string indicating which method to use to determine
#'    the composition of the matrix. Any of "custom", "empirical" or "even" can be used.
#'    See "details"
#' @param fasta A single number or a numeric vector of length equal to
#'    nReplicates. Indicates the random seed to assign the reads to different species
#'    in each sample. If NULL or a single number, the same seed will be applied to all
#'    samples so that they will have the same composition, with differences only due to
#'    random sampling
#' @param genomeReadMatrix Composition matrix, as calculated by the function
#'    \code{compositionGenomesMetaT}.
#' @param modelMatrix A composition matrix predefined by the user, in which rows
#'    represent species and columns represent samples. Id user provides a composition
#'    matrix, this function only performs the random sampling to scale to the number of
#'    reads desired. Defaults to NULL
#' @param DE Logical, whether or not to simulate differential expression between cases
#'    and controls (defaults to FALSE)
#' @param foldChanges Numeric vector, containing the fold changes to simulate
#' @param foldProbs Numeric vector, containing the probabilities for each of the fold-
#'    changes specified in the parameter \code{foldChanges}. See examples.
#' @param nSamples An integer, must be specified if DE is set to TRUE. Number of cases
#'    in the simulated experiment. nSamples + nControls must be equal to the number of
#'    columns in the composition matrix
#' @param nControls An integer, must be specified if DE is set to TRUE. Number of controls
#'    in the simulated experiment. nSamples + nControls must be equal to the number of
#'    columns in the composition matrix
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
# function to calculate reads for one genome, for all samples
simulateSingleGenome=function(genomeName, fasta, genomeReadMatrix, modelMatrix=NULL,
                              DE=F, foldChanges=NULL, foldProbs=NULL, nSamples=NULL,
                              nControls=NULL, seed=42){
  # this first part loads the count matrix from the pasilla dataset,
  # if the user didn't provide a custom one ----
  if(is.null(modelMatrix)){
    modelMatrix <- loadPasilla()
  }

  # now, randomly assign expression from some genes of the pasilla dataset
  # to the genes of the simulated bacteria
  # set.seed(round(as.numeric(genomeID)/1000)) ----
  fastaGenes <- read.fasta(fasta,as.string = T)
  numGenesToSimulate <- length(fastaGenes)
  countsBacteria <- modelMatrix[sample(numGenesToSimulate, replace = TRUE),]

  # get parameters from this sampling
  par <- get_params(countsBacteria)

  # generate new dataset, separate "controls" and "treated" samples and
  # calculate differential expression with the foldChanges and foldProbs provided
  simData <- create_read_numbers(par$mu, par$fit, par$p0, m = numGenesToSimulate, n = ncol(genomeReadMatrix)) ## normally this creates ~5M reads per sample
  colnames(simData) <- c(paste("sample", 1:ncol(simData), sep = "_"))
  rownames(simData) <- names(fastaGenes)

  if(DE == T){
    if(ncol(genomeReadMatrix) != nSamples+nControls){
      print("Please provide correct number of samples and controls so that nSamples+nControls=total number of simulated experiments")
    } else{
      controlSimData <- simData[,1:nControls]
      treatedSimData <- simData[,(nControls+1):(nSamples+nControls)]

      DEgenes <- sample(foldChanges, size = numGenesToSimulate, prob = foldProbs, replace=T) ## sample number of genes up/down-regulated or staying the same (thus with FC=1)

      treatedSimData <- round(treatedSimData*DEgenes)

      simData <- cbind(controlSimData, treatedSimData)
      colnames(simData) <- c(paste("control", 1:nControls, sep="_"), paste("treated", 1:nSamples, sep="_"))
      rownames(simData) <- names(fastaGenes)
      DEgenesNames <- cbind(names(fastaGenes), DEgenes)
    }
  } else{
    DEgenes <- NULL
    DEgenesNames <- NULL
  }

  ##now, downsize each of the samples to the corresponding number of reads
  reducedSimData <- sampleReads(countMatrix=simData,readNumber=c(unlist(genomeReadMatrix[genomeName,]))) ## we downsize the number of reads for this genome, according to the reads of each sample
  reducedSimData[is.na(reducedSimData)] <- 0
  ## make final object collecting all the info from the simulation of this genome
  simulatedGenome=list(simulationData=reducedSimData, DEgenes=DEgenesNames, numSamples=nSamples, numControls=nControls) ##TODO: check seeds and include in the info for reproducibility
  return(simulatedGenome)
}






sampleReads <- function(countMatrix, readNumber, seed=42){ # downsize the number of simulated reads to the actual number per genome
  reducedReadMatrix <- data.frame(row.names = row.names(countMatrix))
  if (length(readNumber) != 1 & length(readNumber) != ncol(countMatrix)){
    print("Error! You must provide an integer or a vector of length equal to the number of columns in countMatrix")
  }
  if (length(readNumber) == 1){
    readNumber <- rep(readNumber, times=ncol(countMatrix))
  }
  for (col in 1:ncol(countMatrix)){
    temp <- sample(rownames(countMatrix), size = readNumber[col],prob=countMatrix[,col], replace = T)
    reducedReadMatrix <- cbind(reducedReadMatrix,table(temp)[rownames(reducedReadMatrix)])
  }
  reducedReadMatrix <- reducedReadMatrix[,c(F,T)]
  colnames(reducedReadMatrix) <- colnames(countMatrix)
  return(reducedReadMatrix)
}




## function to iterate over all genomes to produce final matrix
simulateMetaTranscriptome <- function(genomeFileDir, genomeReadMatrix, modelMatrix=NULL, DE=F, foldChanges=NULL, foldProbs=NULL, nSamples=NULL, nControls=NULL, seed=42){
  simulatedDataSet <- list(simulationData=NULL, DEgenes=NULL, nSamples=nSamples, nControls=nControls)
  for(i in 1:nrow(genomeReadMatrix)){
    genomeName <- rownames(genomeReadMatrix)[i]
    fastaFile <- c(file.path(genomeFileDir, paste(rownames(genomeReadMatrix)[i],".fa", sep="")),
                file.path(genomeFileDir, paste(rownames(genomeReadMatrix)[i],".fasta", sep="")),
                file.path(genomeFileDir, paste(rownames(genomeReadMatrix)[i],".fna", sep="")),
                file.path(genomeFileDir, paste(rownames(genomeReadMatrix)[i],".genes.fa", sep="")),
                file.path(genomeFileDir, paste(rownames(genomeReadMatrix)[i],".genes.fasta", sep="")),
                file.path(genomeFileDir, paste(rownames(genomeReadMatrix)[i],".genes.fna", sep="")))  ## check all possible fasta file extensions
    fastaFile <- fastaFile[file.exists(fastaFile)] ## select only the one existing
    singleGenome <- simulateSingleGenome(genomeName, fasta=fastaFile, genomeReadMatrix, modelMatrix, DE, foldChanges, foldProbs, nSamples, nControls, seed)
    simulatedDataSet$simulationData <- rbind(simulatedDataSet$simulationData, singleGenome$simulationData)
    simulatedDataSet$DEgenes <- bind(simulatedDataSet$DEgenes, singleGenome$DEgenes)
  }
  return(simulatedDataSet)
}






## TODO: set class "simulatedDataSet" to store the info of the simulation and replace the list structure // check how to initialize an empty instance of this class
##  write wrapper for the simulate_experiment_countmat function in the polyester package
simulateFastaReads <- function(genomeFileDir, simulatedDataSet, outdir=".", paired=T, seed=42, distr="empirical", error_model="illumina5", bias="rnaf", strand_specific=T){
  ## if pre-existing, remove fasta file with all genomes
  fastaFullFile <- file.path(genomeFileDir, ".tmpFile_allFastas.fa")
  if(file.exists(fastaFullFile)){
    system(paste("rm", fastaFullFile, sep=" "))
  }
  ## first, generate a fasta from all the fasta files in genomeFileDir
  fastaFiles <- unique(c(list.files(genomeFileDir, pattern = ".fasta", full.names = T),
                      list.files(genomeFileDir, pattern = ".fa", full.names = T),
                      list.files(genomeFileDir, pattern = ".fna", full.names = T),
                      list.files(genomeFileDir, pattern = ".genes.fasta", full.names = T),
                      list.files(genomeFileDir, pattern = ".genes.fa", full.names = T),
                      list.files(genomeFileDir, pattern = ".genes.fna", full.names = T)))
  for(file in fastaFiles){
    system(paste("cat", file, ">>", fastaFullFile, sep=" "))
  }
  countMatrix <- as.matrix(simulatedDataSet$simulationData)
  simulate_experiment_countmat(fasta=fastaFullFile, readmat=countMatrix, outdir=outdir, paired=paired, seed=seed, distr=distr, error_model=error_model, bias=bias, strand_specific=T)
  system(paste("rm", fastaFullFile, sep=" "))
  out_countMatrix <- file.path(outdir, "countMatrix.txt")
  out_DEgenes <- file.path(outdir, "DEgenes.txt")
  write.table(simulatedDataSet$simulationData, file=out_countMatrix, col.names=T, row.names=T, quote=F, sep="\t")
  write.table(simulatedDataSet$DEgenes, file=out_DEgenes, col.names=F, row.names=F, quote=F, sep="\t")
}



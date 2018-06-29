

# function to calculate reads for one genome, for all samples
simulateSingleGenome=function(genomeName, fasta, genomeReadMatrix, modelMatrix=NULL, DE=F, foldChanges=NULL, foldProbs=NULL, nSamples=NULL, nControls=NULL, seed=42){
  # this first part loads the count matrix from the pasilla dataset,
  # if the user didn't provide a custom one ----
  if(is.null(modelMatrix)){
    modelMatrix <- loadPasilla()
  }

  # now, randomly assign expression from some genes of the pasilla dataset to the genes of the simulated bacteria
  # set.seed(round(as.numeric(genomeID)/1000)) ----
  fastaGenes <- read.fasta(fasta,as.string = T)
  numGenesToSimulate <- length(fastaGenes)
  countsBacteria <- modelMatrix[sample(numGenesToSimulate, replace = TRUE),]

  ## get parameters from this sampling
  par <- get_params(countsBacteria)

  ## generate new dataset, separate "controls" and "treated" samples and calculate differential expression with the foldChanges and foldProbs provided
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



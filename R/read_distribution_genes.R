
#' Calculate reads for one genome, for all samples
#'
#' \code{simulateSingleGenome} simulates a gene count matrix for a single
#'    genome. This is a wrapper for the function \code{create_read_numbers} in
#'    the polyester package. Users should not call it directly but should only
#'    be called within the \code{simulateMetaTranscriptome} function
#' @param genomeName Name of the genome to be simulated. It should match the
#'    basename of the fasta file for this genome
#' @param fasta Either the path of the multi-fasta file containing the sequences
#'    of all genes of the selected genome, or the fasta object read with the
#'    function \code{read.fasta} from the \code{seqinr} package.
#'    The names of the genes will be taken from the fasta headers
#' @param genomeReadMatrix Microbial composition matrix containing the number of
#'    reads per genome and per sample. Can be obtained using the function
#'    \code{compositionGenomesMetaT}
#' @param modelMatrix A composition matrix of gene expression, in which rows
#'    represent genes and columns represent replicates. User can provide one of
#'    their own, otherwise the matrix from the Pasilla dataset will be used.
#'    It's usedto fit a zero-inflated negative binomial and set the parameters
#'    to randomly assign gene expression to the genes from the microbial genome.
#' @param DE Logical, whether or not to simulate differential expression between
#'    cases and controls (defaults to FALSE)
#' @param foldChanges Numeric vector, containing the fold changes to simulate.
#'    It should contain the value 1, for genes which are not differentially
#'    expressed. Required if DE set to TRUE
#' @param foldProbs Numeric vector, containing the probabilities for each of the
#'    fold-changes specified in the parameter \code{foldChanges}. Required if DE
#'    is set to TRUE. See examples
#' @param nSamples An integer, must be specified if DE is set to TRUE. Number of
#'    cases in the simulated experiment. nSamples + nControls must be equal to
#'    the number of columns in the composition matrix \code{genomeReadMatrix}
#' @param nControls An integer, must be specified if DE is set to TRUE. Number
#'    of controls in the simulated experiment. nSamples + nControls must be
#'    equal to the number of columns in the composition matrix
#'    \code{genomeReadMatrix}
#' @param seed An integer, sets the random seed for the read distribution.
#' @details This function is a wrapper for the \code{create_read_numbers}
#'    function from the polyester package. Takes a gene count matrix to fit a
#'    zero-inflated negative binomial distribution, then uses this as a model to
#'    randomly assign gene expression to each of the genes from the genome, as
#'    specified in a fasta file
#' @references
#'    - Huber W, Reyes A (2018). pasilla: Data package with per-exon and
#'    per-gene read counts of RNA-seq samples of Pasilla knock-down by Brooks et
#'    al. R package version 1.8.0
#'    - Alyssa C. Frazee, Andrew E. Jaffe, Rory Kirchner and Jeffrey T. Leek
#'    (2018). polyester: Simulate RNA-seq reads. R package version 1.16.0.
#' @import stats
#' @import utils
#' @import polyester
#' @import seqinr
#' @export
#' @return A list, containing the following elements:
#'    - simulationData: a data.frame with the read counts for each gene and each sample.
#'      Each row represents a gene and each column a sample. If there is differential
#'      expression, column names indicate whether each sample is a case or a control
#'    - numSamples: if DE is set to TRUE, the number of cases specified, otherwise NULL
#'    - numControls: if DE is set to TRUE, the number of controls, otherwise NULL
#'    - DEgenes: if DE is set to TRUE, a two-column data.frame, the first column
#'      indicating gene names and the second column the fold change applied to each gene
#' @examples
#' # First, define a list of genomes to simulate. The names of these genomes need to match
#' # the names of the fasta files (without the extension). The genomes used are:
#' # - F. prausnitzii
#' # - R. intestinalis
#' # - L. johnsonii
#' # - E. faecalis
#' # - B. obeum
#' genomesToSimulate <- c("fprausnitzii", "rintestinalis", "ljohnsonii", "efaecalis",
#'                        "bobeum")
#'
#' # Then, obtain the empirical composition matrix for this 5 species
#' compMatrix <- compositionGenomesMetaT(composition="empirical", empiricalSeed=1,
#'                                    genomes=genomesToSimulate, nReads=500000,
#'                                    nReplicates=10)
#'
#' # Obtain the gene expression matrix for E. faecalis, assuming no differential
#' # expression. No composition matrix is provided, so the one from the pasilla dataset
#' # will be used
#' data(efaecalis)
#' countMatrix_efaecalis <- simulateSingleGenome(genomeName="efaecalis", fasta=efaecalis,
#'                                            genomeReadMatrix=compMatrix)
#'
#' # Obtain the gene expression matrix for B. obeum, incorporating differential
#' # expression: 10% genes have a 2-fold overexpression and 10% have a 0.5-fold depletion.
#' # No composition matrix is provided, so the one from the pasilla dataset will be used.
#' # As there are 10 samples in the count matrix, we assign 5 cases and 5 controls.
#' data(bobeum)
#' countMatrix_bobeum <- simulateSingleGenome(genomeName="bobeum", fasta=bobeum,
#'                                         genomeReadMatrix=compMatrix, DE=TRUE,
#'                                         foldChanges=c(0.5,1,2), foldProbs=c(10,80,10),
#'                                         nSamples=5, nControls=5)
# function to calculate reads for one genome, for all samples
simulateSingleGenome=function(genomeName, fasta, genomeReadMatrix, modelMatrix=NULL,
                              DE=F, foldChanges=NULL, foldProbs=NULL, nSamples=NULL,
                              nControls=NULL, seed=42){
  # this first part loads the count matrix from the pasilla dataset,
  # if the user didn't provide a custom one ----
  if(is.null(modelMatrix)){
    modelMatrix <- loadPasilla()
  }
  # if fasta is a string, it should be the path to the fasta file
  # if it is a list, it should be
  if(is.character(fasta)){
    fastaGenes <- read.fasta(fasta,as.string = TRUE)
  } else if(is.list(fasta)){
    fastaGenes <- fasta
  } else{
    print("Please provide a valid fasta argument. Check help for ?simulateSingleGenome")
  }
  # now, randomly assign expression from some genes of the pasilla dataset
  # to the genes of the simulated bacteria
  numGenesToSimulate <- length(fastaGenes)
  countsBacteria <- modelMatrix[sample(1:nrow(modelMatrix), size=numGenesToSimulate, replace = TRUE),]
  ## While there are less than four unique 'x' values in countsBacteria, repeat sampling
  while( length(unique(countsBacteria$untreated1))<4 || length(unique(countsBacteria$untreated2))<4 || length(unique(countsBacteria$untreated3))<4 || length(unique(countsBacteria$untreated4))<4){
	countsBacteria <- modelMatrix[sample(1:nrow(modelMatrix), size=numGenesToSimulate, replace = TRUE),]
	}

  # get parameters from this sampling
  par <- get_params(countsBacteria)

  # generate new dataset, separate "controls" and "treated" samples and
  # calculate differential expression with the foldChanges and foldProbs provided
  # normally this creates ~5M reads per sample
  simData <- create_read_numbers(par$mu, par$fit, par$p0, m = numGenesToSimulate,
                                 n = ncol(genomeReadMatrix))
  colnames(simData) <- c(paste("sample", 1:ncol(simData), sep = "_"))
  rownames(simData) <- names(fastaGenes)

  if(DE == T){
    if(ncol(genomeReadMatrix) != nSamples+nControls){
      print("Please provide correct number of samples and controls so that nSamples+nControls=total number of simulated experiments")
    } else{
      controlSimData <- simData[,1:nControls]
      treatedSimData <- simData[,(nControls+1):(nSamples+nControls)]

      # sample number of genes up/down-regulated or staying the same (thus with FC=1)
      DEgenes <- sample(foldChanges, size = numGenesToSimulate, prob = foldProbs,
                        replace=T)

      treatedSimData <- round(treatedSimData*DEgenes)

      simData <- cbind(controlSimData, treatedSimData)
      colnames(simData) <- c(paste("control", 1:nControls, sep="_"),
                             paste("treated", 1:nSamples, sep="_"))
      rownames(simData) <- names(fastaGenes)
      DEgenesNames <- as.data.frame(cbind(names(fastaGenes), DEgenes), stringsAsFactors=F)
    }
  } else{
    DEgenes <- NULL
    DEgenesNames <- NULL
  }

  # now, downsize each of the samples to the corresponding number of reads
  # we downsize the number of reads for this genome, according to the reads of each sample
  reducedSimData <- sampleReads(countMatrix=simData,
                                readNumber=c(unlist(genomeReadMatrix[genomeName,])))
  reducedSimData[is.na(reducedSimData)] <- 0
  # make final object collecting all the info from the simulation of this genome
  simulatedGenome=list(simulationData=reducedSimData, DEgenes=DEgenesNames,
                       numSamples=nSamples, numControls=nControls)
  ##TODO: check seeds and include in the info for reproducibility
  return(simulatedGenome)
}








#' Calculate reads for one genome, for all samples
#'
#' \code{simulateMetaTranscriptome} simulates a gene count matrix for an entire
#'    metatranscriptome
#' @param genomeFileDir Character string indicating the location of the fasta files
#'    for all genomes to be included in the metatranscriptome simulation. The basenames
#'    of these fasta files must match the rownames of the genomeReadMatrix composition
#'    matrix. See details
#' @param genomeReadMatrix Microbial composition matrix containing the number of reads
#'    per genome and per sample. Can be obtained using the function
#'    \code{compositionGenomesMetaT}
#' @param modelMatrix A composition matrix of gene expression, in which rows represent
#'    genes and columns represent replicates. User can provide one of their own,
#'    otherwise the matrix from the Pasilla dataset will be used. It's used to fit a
#'    zero-inflated negative binomial and set the parameters to randomly assign gene
#'    expression to the genes from the microbial genome.
#' @param DE Logical, whether or not to simulate differential expression between cases
#'    and controls (defaults to FALSE)
#' @param foldChanges Numeric vector, containing the fold changes to simulate. It should
#'    contain the value 1, for genes which are not differentially expressed. Required if
#'    DE set to TRUE
#' @param foldProbs Numeric vector, containing the probabilities for each of the fold-
#'    changes specified in the parameter \code{foldChanges}. Required if DE is set to
#'    TRUE. See examples
#' @param nSamples An integer, must be specified if DE is set to TRUE. Number of cases
#'    in the simulated experiment. nSamples + nControls must be equal to the number of
#'    columns in the composition matrix \code{genomeReadMatrix}
#' @param nControls An integer, must be specified if DE is set to TRUE. Number of controls
#'    in the simulated experiment. nSamples + nControls must be equal to the number of
#'    columns in the composition matrix \code{genomeReadMatrix}
#' @param seed An integer, sets the random seed for the read distribution.
#' @details This function iterates over all the genomes present in the composition matrix
#'    and simulates their corresponding gene expression matrix, putting them all together
#'    Valid fasta extensions for the fasta files located in  \code{genomeFileDir}:
#'    *.fa, *.fasta, *.fna, *.genes.fa, *.genes.fasta, *.genes.fna
#' @references
#'    - Huber W, Reyes A (2018). pasilla: Data package with per-exon and per-gene read
#'    counts of RNA-seq samples of Pasilla knock-down by Brooks et al. R package version
#'    1.8.0
#'    - Alyssa C. Frazee, Andrew E. Jaffe, Rory Kirchner and Jeffrey T. Leek (2018).
#'    polyester: Simulate RNA-seq reads. R package version 1.16.0.
#' @import stats
#' @import utils
#' @import polyester
#' @import seqinr
#' @export
#' @return A list, containing the following elements:
#'    - simulationData: a data.frame with the read counts for each gene and each sample.
#'      Each row represents a gene and each column a sample. If there is differential
#'      expression, column names indicate whether each sample is a case or a control
#'    - numSamples: if DE is set to TRUE, the number of cases specified, otherwise NULL
#'    - numControls: if DE is set to TRUE, the number of controls, otherwise NULL
#'    - DEgenes: if DE is set to TRUE, a two-column data.frame, the first column
#'      indicating gene names and the second column the fold change applied to each gene
#' @examples
#' # First, define a list of genomes to simulate. The names of these genomes need to match
#' # the names of the fasta files (without the extension). The genomes used are:
#' # - F. prausnitzii
#' # - R. intestinalis
#' # - L. johnsonii
#' # - E. faecalis
#' # - B. obeum
#' genomesToSimulate <- c("fprausnitzii", "rintestinalis", "ljohnsonii", "efaecalis",
#'                        "bobeum")
#'
#' # Then, obtain the empirical composition matrix for this 5 species
#' compMatrix <- compositionGenomesMetaT(composition="empirical", empiricalSeed=1,
#'                                    genomes=genomesToSimulate, nReads=500000,
#'                                    nReplicates=10)
#'
#'
#' # Obtain the gene expression matrix for the full community (metatranscriptome)
#' # In this case, there is no differential expression in any of the bacteria.
#' # No composition matrix is provided, so the one from the pasilla dataset will be used.
#' # For this, first indicate the location of the fasta files
#' genomesFolder = system.file("extdata", package = "metaester", mustWork = TRUE)
#' metatranscriptome <- simulateMetaTranscriptome(genomeFileDir=genomesFolder,
#'                                                genomeReadMatrix=compMatrix)
#'
#' # Obtain the gene expression matrix for the full community (metatranscriptome)
#' # incorporating differential expression: 10% genes (in each bacterium) have a 2-fold
#' # overexpression and 10% have a 0.5-fold depletion.
#' # No composition matrix is provided, so the one from the pasilla dataset will be used.
#' # As there are 10 samples in the count matrix, we assign 5 cases and 5 controls.
#' metatranscriptome <- simulateMetaTranscriptome(genomeFileDir=genomesFolder,
#'                                                genomeReadMatrix=compMatrix, DE=TRUE,
#'                                                foldChanges=c(0.5,1,2),
#'                                                foldProbs=c(10,80,10),
#'                                                nSamples=5, nControls=5)
simulateMetaTranscriptome <- function(genomeFileDir, genomeReadMatrix, modelMatrix=NULL,
                                      DE=F, foldChanges=NULL, foldProbs=NULL,
                                      nSamples=NULL, nControls=NULL, seed=42){
  simulatedDataSet <- list(simulationData=NULL, DEgenes=NULL, nSamples=nSamples,
                           nControls=nControls)
  ## if not all the genomes in genomeReadMatrix are used, remove the genome from matrix
  genomeReadMatrix[,"sum"] <- rowSums(genomeReadMatrix)
  genomeReadMatrix <- genomeReadMatrix[genomeReadMatrix$sum != 0,-ncol(genomeReadMatrix)]
  for(i in 1:nrow(genomeReadMatrix)){
    genomeName <- rownames(genomeReadMatrix)[i]
    # check all possible fasta file extensions
    fastaFile <- c(file.path(genomeFileDir,
                             paste(rownames(genomeReadMatrix)[i],".fa", sep="")),
                  file.path(genomeFileDir,
                          paste(rownames(genomeReadMatrix)[i],".fasta", sep="")),
                  file.path(genomeFileDir,
                          paste(rownames(genomeReadMatrix)[i],".fna", sep="")),
                  file.path(genomeFileDir,
                          paste(rownames(genomeReadMatrix)[i],".genes.fa", sep="")),
                  file.path(genomeFileDir,
                          paste(rownames(genomeReadMatrix)[i],".genes.fasta", sep="")),
                  file.path(genomeFileDir,
                          paste(rownames(genomeReadMatrix)[i],".genes.fna", sep="")))
    fastaFile <- fastaFile[file.exists(fastaFile)] ## select only the ones existing
    singleGenome <- simulateSingleGenome(genomeName, fasta=fastaFile, genomeReadMatrix,
                                         modelMatrix, DE, foldChanges, foldProbs,
                                         nSamples, nControls, seed)
    simulatedDataSet$simulationData <- rbind(simulatedDataSet$simulationData,
                                             singleGenome$simulationData)
    simulatedDataSet$DEgenes <- rbind(simulatedDataSet$DEgenes, singleGenome$DEgenes)
  }
  return(simulatedDataSet)
}






# TODO: set class "simulatedDataSet" to store the info of the simulation and replace the
# list structure // check how to initialize an empty instance of this class
# TODO: allow to generate the simulations with the entire genome and the *gff annotation
# file too (currently only with the genes fasta). Does polyester have this function?
# If so, write a wrapper for it.
##  write wrapper for the simulate_experiment_countmat function in the polyester package

#' Generate fasta files for a metatranscriptome gene expression matrix
#'
#' \code{simulateFastaReads} generates the fasta files for a given metatranscriptome
#'    gene expression matrix, calculated with \code{simulateMetaTranscriptome}. The
#'    fasta files are written to a directory. This is a wrapper to the
#'    \code{simulate_experiment_countmat} function from the polyester package
#' @param genomeFileDir Character string indicating the location of the fasta files
#'    for all genomes to be included in the metatranscriptome simulation. The basenames
#'    of these fasta files must match the rownames of the genomeReadMatrix composition
#'    matrix. See details
#' @param simulatedDataSet Metatranscriptome gene expression matrix containing the
#'    read counts per gene and per sample. Can be obtained using the function
#'    \code{simulateMetaTranscriptome}
#' @param outdir Path to the directory where the fasta files with the simulated reads
#'    should be written. Created if not exists.
#' @param paired Logical, whether to generate paired-end fasta reads (defaults to TRUE)
#' @param seed An integer, sets the random seed for the fasta generation
#' @param distr Parameter to be passed to the \code{simulate_experiment_countmat}
#'    function. Controls the distribution of fragment length for paired-end reads.
#'    Defaults to "empirical". Check documentation in the polyester package.
#' @param error_model Parameter to be passed to the \code{simulate_experiment_countmat}
#'    function. Used to simulate error rates from the sequence. Defaults to "illumina5".
#'    Check documentation in the polyester package.
#' @param bias Parameter to be passed to the \code{simulate_experiment_countmat}
#'    function. Controls the distribution of reads within the gene/transcript, which
#'    may be different depending on the fragmentation step of the library preparation
#'    protocol (leading to increased coverage at the 5' or 3' of the transcript).
#'    Defaults to "rnaf". Check documentation in the polyester package.
#' @param strand_specific Logical. Parameter to be passed to the
#'    \code{simulate_experiment_countmat} function. Should the reads be strand-specific?
#'    Check documentation in the polyester package.
#' @details  This function is a wrapper to the \code{simulate_experiment_countmat}
#'    function from the polyester package. Reads can be simulated from a multifasta
#'    file containing one sequence for each gene. Functionality to work with full genomes
#'    + GTF annotation files is not developed. Meanwhile, there are several (external)
#'    tools that allow to convert a genome fasta + a GTF file into a transcript fasta,
#'    such as \code{gffread} from Cufflinks.
#' @references
#'    - Alyssa C. Frazee, Andrew E. Jaffe, Rory Kirchner and Jeffrey T. Leek (2018).
#'    polyester: Simulate RNA-seq reads. R package version 1.16.0.
#' @import stats
#' @import utils
#' @import polyester
#' @import seqinr
#' @export
#' @return Writes the following files into the specified output directory:
#'    - Fasta files (*.fasta), one for each sample, containing the reads per gene
#'      specified in the metatranscriptomics gene count matrix (simulatedDataSet object)
#'    - Gene count matrix (countMatrix.txt): matrix with the counts per gene and per
#'      sample, useful to compare the results of mapping/differential expression tools to
#'    - Differentially expressed genes list (DEgenes.txt): if DE is set to TRUE, a
#'      two-column matrix, the first column indicating gene names and the second column
#'      the fold change applied to each gene. Otherwise the file will be empty.
#' @examples
#' # First, define a list of genomes to simulate. The names of these genomes need to match
#' # the names of the fasta files (without the extension). The genomes used are:
#' # - F. prausnitzii
#' # - R. intestinalis
#' # - L. johnsonii
#' # - E. faecalis
#' # - B. obeum
#' genomesToSimulate <- c("fprausnitzii", "rintestinalis", "ljohnsonii", "efaecalis",
#'                        "bobeum")
#'
#' # Then, obtain the empirical composition matrix for this 5 species
#' compMatrix <- compositionGenomesMetaT(composition="empirical", empiricalSeed=1,
#'                                    genomes=genomesToSimulate, nReads=500000,
#'                                    nReplicates=10)
#'
#' # Obtain the gene expression matrix for the full community (metatranscriptome)
#' # incorporating differential expression: 10% genes (in each bacterium) have a 2-fold
#' # overexpression and 10% have a 0.5-fold depletion.
#' # No composition matrix is provided, so the one from the pasilla dataset will be used.
#' # As there are 10 samples in the count matrix, we assign 5 cases and 5 controls.
#' genomesFolder = system.file("extdata", package = "metaester", mustWork = TRUE)
#' metatranscriptome <- simulateMetaTranscriptome(genomeFileDir=genomesFolder,
#'                                                genomeReadMatrix=compMatrix, DE=TRUE,
#'                                                foldChanges=c(0.5,1,2),
#'                                                foldProbs=c(10,80,10),
#'                                                nSamples=5, nControls=5)
#'
#' # Finally, generate the fasta files and write them to the output directory
#' \donttest{simulateFastaReads(genomeFileDir=genomesFolder, simulatedDataSet=metatranscriptome,
#'                    outdir=".")}
simulateFastaReads <- function(genomeFileDir, simulatedDataSet, genomeReadMatrix, outdir=".",
                               paired=T, seed=42, distr="empirical", error_model="illumina5",
                               bias="rnaf", strand_specific=T){
  ## if not all the genomes in genomeReadMatrix are used, remove the genomes from matrix
  genomeReadMatrix[, "sum"] <- rowSums(genomeReadMatrix)
  genomeReadMatrix <- genomeReadMatrix[genomeReadMatrix$sum != 0, -ncol(genomeReadMatrix)]
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

  # create output directory if not exists
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  # simulate fasta readsand write output to outdir
  simulate_experiment_countmat(fasta=fastaFullFile, readmat=countMatrix, outdir=outdir,
                               paired=paired, seed=seed, distr=distr,
                               error_model=error_model, bias=bias, strand_specific=T)
  system(paste("rm", fastaFullFile, sep=" "))
  out_countMatrix <- file.path(outdir, "countMatrix.txt")
  out_DEgenes <- file.path(outdir, "DEgenes.txt")
  write.table(simulatedDataSet$simulationData, file=out_countMatrix, col.names=T,
              row.names=T, quote=F, sep="\t")
  write.table(simulatedDataSet$DEgenes, file=out_DEgenes, col.names=F, row.names=F,
              quote=F, sep="\t")
}





# function to calculate reads per genome per sample
# TODO: important - make nReads possibly a vector, so that experiment can have a different sequencing depth

compositionGenomesMetaT <- function(composition, empiricalSeed=NULL, genomes, compositionMatrix=NULL, nReads=1000000, nReplicates=10, seed=42){
  genomeReadMatrix <- data.frame(row.names = genomes)

  # generate or check that the user has provided a composition matrix ----
  if(composition == "custom" & is.null(compositionMatrix)){
    print("Please provide a composition matrix, with rows for each genome and columns for each sample")
  } else if(composition == "empirical"){
    if(is.null(empiricalSeed)){
      empiricalSeed <- rep(1, times=nReplicates)
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

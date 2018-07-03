# creates data object for the species distribution

speciesDistribution=get(load("inst/extdata/speciesDist.RData"))
devtools::use_data(speciesDistribution, compress = "xz", internal = TRUE, overwrite = T)


# creates data objects for each of the 5 example genomes

efaecalis=read.fasta("inst/extdata/efaecalis.genes.fna",as.string = TRUE)
fprausnitzii=read.fasta("inst/extdata/fprausnitzii.genes.fna",as.string = TRUE)
rintestinalis=read.fasta("inst/extdata/rintestinalis.genes.fna",as.string = TRUE)
ljohnsonii=read.fasta("inst/extdata/ljohnsonii.genes.fna",as.string = TRUE)
bobeum=read.fasta("inst/extdata/bobeum.genes.fna",as.string = TRUE)
devtools::use_data(efaecalis, fprausnitzii, rintestinalis, ljohnsonii, bobeum, overwrite = T)

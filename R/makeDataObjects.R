# creates data object for the species distribution

speciesDistribution=get(load("data-raw/speciesDist.RData"))
devtools::use_data(speciesDistribution, compress = "xz", internal = TRUE, overwrite = T)


# creates data objects for each of the 5 example genomes
efaecalis=read.fasta("data-raw/efaecalis.genes.fna",as.string = TRUE)
fprausnitzii=read.fasta("data-raw/fprausnitzii.genes.fna",as.string = TRUE)
rintestinalis=read.fasta("data-raw/rintestinalis.genes.fna",as.string = TRUE)
ljohnsonii=read.fasta("data-raw/ljohnsonii.genes.fna",as.string = TRUE)
bobeum=read.fasta("data-raw/bobeum.genes.fna",as.string = TRUE)
devtools::use_data(efaecalis, fprausnitzii, rintestinalis, ljohnsonii, bobeum, overwrite = T)

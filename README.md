# Introduction

Metaester is an R package to simulate metatranscriptomics datasets. 

Metaester is built from [Polyester](https://github.com/alyssafrazee/polyester), an R package to simulate RNA-seq experiments, hence its name. Many of the functions in Metaester are wrappers to Polyester functions. However, we have not exploited all the possible functionalities of polyester yet. 

Metaester allows to simulate an entire microbial community, by first distributing the RNA-seq reads among the different genomes present in the simulated sample and then distributing the RNA-seq reads among the different genes, within each simulated genome.


# Installation

Metaester is built upon the devel version of Polyester (as of July 2018, not all changes in the devel version were pushed to the bioC version). Therefore, to use Metaester, we strongly recommend to install the devel version of Polyester:
```{r installPolyester, eval=FALSE, echo=TRUE}
devtools::install_github("alyssafrazee/polyester")
```

To install Metaester, run the following code:
```{r installMetaester, eval=FALSE, echo=TRUE}
devtools::install_github("vllorens/metaester")
```

# The metaester workflow

Metaester is built with the aim of simulating metatranscriptomics datasets as simply as possible. Only 4 steps are required to simulate an experiment:

- Step 0: Define the parameters of the experiment: number of reads, replicates, number of cases and controls for differential expression, etc.
- Step 1: Generate genome read distribution: distribute the total reads *per sample* among the genomes to simulate
- Step 2: For each genome, generate gene read distribution: distribute the total reads *per genome* among the genes of that genome
- Step 3: Generate the fasta files with the raw reads of the simulation


# Metaester examples

See all available examples in the Metaester [vignette](https://github.com/vllorens/metaester/blob/master/vignettes/metaester.Rmd)

# Contribute

External contributions are highly appreciated. If you feel like a new functionality could be added, or found a bug in the code, please feel free to open an issue.

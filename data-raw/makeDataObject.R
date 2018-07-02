# creates data object for the species distribution

speciesDistribution=get(load("data-raw/speciesDist.RData"))
devtools::use_data(speciesDistribution, compress = "xz", internal = TRUE, overwrite = T)

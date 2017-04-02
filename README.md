
# Required Libraries


# Install required libraries in this suggested order:

source("http://bioconductor.org/biocLite.R")

biocLite("GenomicRanges")

biocLite("GenomicFeatures")

biocLite("rtracklayer")

biocLite("Gviz")

biocLite("biomaRt")


# Load libraries

library(shiny)

library(GenomicFeatures)

library(rtracklayer)

library(Gviz)

library(biomaRt)

runGitHub("test1","Marisol-AMtz")

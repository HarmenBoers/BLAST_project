# Parse arguments
args = commandArgs(TRUE)
uniprotfile = normalizePath(args[1])
directory = dirname(uniprotfile)
setwd(directory)
# Set seed for random selection for reproducible results
set.seed(1)
#Select a random set of GO terms using the 'sample' function in R with a list of uniprot identifiers
uniprot <- read.table(uniprotfile)
random <- sample(uniprot$V1, 1000)
write(as.vector(random), file="./random_1000_ids.txt")

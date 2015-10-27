source("../ChainRank.R")
outfile <- "ChainRank_results2315.txt"

network <- read.delim("MuscleSpecCase_Network.txt",stringsAsFactors=FALSE)
elements <- read.delim("MuscleSpecCase_Scores.txt",stringsAsFactors=FALSE)
candidates <- read.delim("MuscleSpecCase_Candidates.txt", header=F, dec=",", colClasses = "character",stringsAsFactors=FALSE)[,1]
targets <- read.delim("MuscleSpecCase_Targets.txt", header=F, dec=",", colClasses = "character",stringsAsFactors=FALSE)[,1]
validation <- read.delim("MuscleSpecCase_GoldStandard.txt", header=F, dec=",", colClasses = "character",stringsAsFactors=FALSE)[,1]
maxLength <- 5

results <- RunChainSearch(network,elements,candidates,targets,maxDepth=maxLength,outfile,RetRanks=T)

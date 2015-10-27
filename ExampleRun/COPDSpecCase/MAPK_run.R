source("ChainRank.R")
options(stringsAsFactors = FALSE)

outfile <- "ChainRank"
network <- read.delim("COPDSpecCase_Network.txt",header=T,colClasses = "character")
elements <- read.delim("COPDSpecCase_Scores.txt")
candidates <- read.delim("COPDSpecCase_Candidates.txt", header=F, dec=",", colClasses = "character")[,1]
targets <- read.delim("COPDSpecCase_Targets.txt", header=F, dec=",", colClasses = "character")[,1]
validation <- read.delim("COPDSpecCase_GoldStandard.txt", header=F, dec=",", colClasses = "character")  [,1]
maxLength <- 7

out <- RunChainSearch(network,elements,candidates,targets,maxDepth=maxLength,outfile,RetRanks=TRUE)

mainDir <- ""
source(paste(mainDir,"ChainRank.R",sep=""))
outfile <- paste(mainDir,"Chain_Search",sep="")
network <- read.delim(paste("IgfAkt_HubReducing_Network.txt",sep=""))
elements <- read.delim(paste("IgfAkt_HubReducing_Scores.txt",sep=""))
pathways <- elements[,c(1,9)]
candidates <- read.delim("IgfAkt_Candidates.txt", header=T, dec=",", colClasses = "character")
targets <- read.delim("IgfAkt_Targets.txt", header=T, dec=",", colClasses = "character")
validation <- read.delim("IgfAkt_GoldStandard.txt", header=T, dec=",", colClasses = "character")  

maxdepth <- 6
networkDim <- c(3,5)
pathwayDim <- 2
scoreDim <- 4:6

# preprocess data, missingness
# can't delete element because it's in network, wouldn't find it. Instead = 0
missing <- which(elements[,scoreDim] == "" | is.na(elements[,scoreDim]),arr.ind=TRUE)
if(nrow(toMatrix(missing))){
  for(i in 1:nrow(toMatrix(missing))){
    col <- (scoreDim[1] + missing[i,2] - 1)
    elem[missing[i,1],col] <- 0.0
  }
}

# delete non protein interactions
network <- network[network[,4] == "Protein" & network[,6] == "Protein",]

out <- RunChainSearch(network,elements,pathways,candidates[,2],targets[,2],
                      networkDim,pathwayDim,scoreDim,outfile,maxDepth=maxdepth)

results <- out[[1]]

# Attach candidates and targets to the chain
results[,2] <- paste(results[,1],results[,2],results[,3],sep="|")


# Default validation and analysis results (network ,protein frequency, combined rank of the top 10 chains)
ShortValidationChain(results[,c(2,5:7)],file=mainDir,validation[,1],pathways,results[,c(2,5:7)])

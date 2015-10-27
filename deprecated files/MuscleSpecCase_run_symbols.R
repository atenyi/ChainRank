source("../ChainRank.R")
outfile <- "ChainRank"
# network <- read.delim("IgfAkt_HubEnriching_Network.txt")
# elements <- read.delim("IgfAkt_HubEnriching_Scores.txt")
network <- read.delim("MuscleSpecCase_Network.txt",stringsAsFactors=FALSE)
elements <- read.delim("MuscleSpecCase_Scores.txt",stringsAsFactors=FALSE)
candidates <- read.delim("MuscleSpecCase_Candidates.txt", header=F, dec=",", colClasses = "character",stringsAsFactors=FALSE)[,1]
targets <- read.delim("MuscleSpecCase_Targets.txt", header=F, dec=",", colClasses = "character",stringsAsFactors=FALSE)[,1]
validation <- read.delim("MuscleSpecCase_GoldStandard.txt", header=F, dec=",", colClasses = "character",stringsAsFactors=FALSE)[,1]
maxLength <- 5

set.seed(314)
out <- RunChainSearch(network,elements,candidates,targets,maxDepth=maxLength,outfile,RetRanks=T)

#############################
# Extended validation
#############################
results <- read.delim("ChainRank_results.txt")
results <- out[[1]]

# reverse Connectivity ranking so the smaller ranked chains are preferred
revTopResults <- data.frame(Chain=results[,2],Connectivity.score.pval=((results[,5]*-1)+1),results[,6:9])

# Validation of scores
set.seed(19880301)
ShortValidationChain(Chains=revTopResults[,1],Pvals=revTopResults[,-1],goldStandard=validation,results[,c(2,5:9)])
# ShortValidationChain(revTopResults,file=paste(dir,subdir,sep=""),validation[,1],pw,results[,c(2,5:7)])

############################
# Combined scores

pvalThreshold <- 0.05

# Filtering: Connectity
rankedResult <- FilterByRank(revTopResults,order=c(2),type=1,threshold=pvalThreshold)
ShortValidationChain(Chains=rankedResult[,1],Pvals=rankedResult[,-1],goldStandard=validation,valData=results[,c(2,5:9)])

# Intersection: localisation + connectivity
rankedResult <- FilterByRank(revTopResults,order=c(2,4),type=2,threshold=pvalThreshold)  
ShortValidationChain(Chains=rankedResult[,1],Pvals=rankedResult[,-1],goldStandard=validation,valData=results[,c(2,5:9)])
# Intersection: localisation + connectivity + relevance
rankedResult <- FilterByRank(revTopResults,order=c(2:4),type=2,threshold=pvalThreshold)
ShortValidationChain(Chains=rankedResult[,1],Pvals=rankedResult[,-1],goldStandard=validation,valData=results[,c(2,5:9)])
 

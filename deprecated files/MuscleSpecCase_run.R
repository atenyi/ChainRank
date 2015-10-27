source("../ChainRank.r")
outfile <- "ChainRank"
# network <- read.delim("IgfAkt_HubEnriching_Network.txt")
# elements <- read.delim("IgfAkt_HubEnriching_Scores.txt")
network <- read.delim("MuscleSpecCase_Network.txt")
elements <- read.delim("MuscleSpecCase_scores.txt")
candidates <- read.delim(paste(mainDir,"MuscleSpecCase_Candidates.txt",sep=""), header=T, dec=",", colClasses = "character")
targets <- read.delim(paste(mainDir,"MuscleSpecCase_Targets.txt",sep=""), header=T, dec=",", colClasses = "character")
validation <- read.delim(paste(mainDir,"MuscleSpecCase_GoldStandard.txt",sep=""), header=T, dec=",", colClasses = "character")
maxLength <- 6

out <- RunChainSearch(network,scores,candidates,targets,maxDepth=maxLength,outfile)

#############################
# Extended validation
#############################
# results <- read.delim("ChainRank_results.txt")
results <- out[[1]]

# reverse Connectivity ranking so the smaller ranked chains are preferred
revTopResults <- data.frame(Chain=results[,2],Topology.score.pvals=((results[,5]*-1)+1),results[,6:7])

# Validation of scores
ShortValidationChain(Chains=revTopResults[,1],Pvals=revTopResults[,-1],goldStandard=validation,results[,c(2,5:7)])
# ShortValidationChain(revTopResults,file=paste(dir,subdir,sep=""),validation[,1],pw,results[,c(2,5:7)])

############################
# Combined scores

pvalThreshold <- 0.05

# Filtering: Connectity
rankedResult <- FilterByRank(revTopResults,order=c(2),type=1,threshold=pvalThreshold)
ShortValidationChain(Chains=rankedResult[,1],Pvals=rankedResult[,-1],goldStandard=validation,valData=results[,c(2,5:7)])

# Intersection: localisation + connectivity
rankedResult <- FilterByRank(revTopResults,order=c(2,4),type=2,threshold=pvalThreshold)  
ShortValidationChain(Chains=rankedResult[,1],Pvals=rankedResult[,-1],goldStandard=validation,valData=results[,c(2,5:7)])
# Intersection: localisation + connectivity + relevance
rankedResult <- FilterByRank(revTopResults,order=c(2:4),type=2,threshold=pvalThreshold)
ShortValidationChain(Chains=rankedResult[,1],Pvals=rankedResult[,-1],goldStandard=validation,valData=results[,c(2,5:7)])
 

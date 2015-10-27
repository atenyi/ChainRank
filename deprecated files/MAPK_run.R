mainDir <- "../../"
source(paste(mainDir,"ChainRank.R",sep=""))
options(stringsAsFactors = FALSE)

outfile <- "ChainRank"
network <- read.delim(paste("COPDSpecCase_Network_symbol.txt",sep=""),header=T,colClasses = "character")
elements <- read.delim(paste("COPDSpecCase_Scores_symbol.txt",sep=""))
candidates <- read.delim("COPDSpecCase_Candidates_symbol.txt", header=F, dec=",", colClasses = "character")[,1]
targets <- read.delim("COPDSpecCase_Targets_symbol.txt", header=F, dec=",", colClasses = "character")[,1]
validation <- read.delim("COPDSpecCase_GoldStandard_symbol.txt", header=F, dec=",", colClasses = "character")  [,1]
maxLength <- 7

set.seed(314159)
out <- RunChainSearch(network,elements,candidates,targets,maxDepth=maxLength,outfile,RetRanks=TRUE)
# out <- RunChainSearch(network,scores,pathways,candidates[,1],targets[,1],outfile,maxDepth=maxdepth)

#############################
# Extended validation
#############################
results <- read.delim("ChainRank_results.txt",stringsAsFactors=FALSE)
# results <- out[[1]]

# ShortValidationChain(results[,c(2,5:8)],file="",validation[,1],results[,c(2,5:8)],topFunction=10)

# change the ordering of Topological ranking so the smaller ranked chains are prefered
revTopResults <- data.frame(Chain=results[,2],Topology.score.pvals=((results[,5]*-1)+1),results[,6:9])
revTopResults <- data.frame(revTopResults,(results[,10]-max(results[,10]))*-1+1,results[,11:14])
names(revTopResults) <- names(results)[-c(1,3:4)]

# Reverse topology
#ShortValidationChain(revTopResults[,1],revTopResults[,-1],file="pval",validation,results[,c(2,9:12)])
# ShortValidationChain(revTopResults,file=paste(dir,subdir,sep=""),validation[,1],pw,results[,c(2,5:7)])
# set.seed(3141593)
# ShortValidationChain(Chains=revTopResults[,1],Pvals=revTopResults[,2:6],goldStandard=validation,revTopResults[,1:6)])

# Filtering: topol reversed
# rankedResult <- FilterByRank(revTopResults,c(2),1,threshold=0.02)    # topol increasing
# ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="topolFilt_pval",validation,results[,c(2,5:8)])
# # rankedResult <- FilterByRank(results[,c(2,5:8)],c(2),1,threshold=.05)    # topol increasing
# # Filtering: topol reversed + expr
# rankedResult <- FilterByRank(revTopResults,c(2,4),1,threshol=0.02)    # topol increasing
# ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="allFilt_pval",validation,results[,c(2,5:8)])
# # Filtering: expr
# rankedResult <- FilterByRank(revTopResults,c(4),1,threshold=0.05)    # topol increasing
# ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="exprFilt_pval",validation,results[,c(2,10:14)])


# # Intersection: expr + topol reversed
# # rankedResult <- FilterByRank(revTopResults,c(4,2),2)    # topol increasing
# # ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="topExprIntersect",validation,results[,c(2,5:8)],Ranks=rankedResult[,6:9],threshold=100)
# # Intersection: expr + topol reversed + relev
# rankedResult <- FilterByRank(revTopResults,c(2:4),2,threshold=0.02)    # topol increasing
# ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="allIntersect_pval",validation,results[,c(2,5:8)])
# # Intersection: expr + relev
#rankedResult <- FilterByRank(revTopResults,c(3:4),2,threshold=0.02)    # topol increasing
#ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="exprIntersect_pval",validation,results[,c(2,5:8)])



#######################################
## Rank validation
set.seed(3141593)
ShortValidationChain(revTopResults[,1],revTopResults[,2:5],file="ranked",validation,results[,c(2,9:12)],Ranks=revTopResults[,6:9])

# Filtering: topol reversed
#rankedResult <- FilterByRank(revTopResults,c(6),1,threshold='q')    # topol increasing
#ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="topolFilt_ranked",validation,results[,c(2,9:12)],Ranks=rankedResult[,6:9])
# rankedResult <- FilterByRank(results[,c(2,5:8)],c(2),1,threshold=.05)    # topol increasing
# Filtering: topol reversed + expr
#rankedResult <- FilterByRank(revTopResults,c(6,8),1,threshold='q')    # topol increasing
#ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="allFilt_ranked",validation,results[,c(2,9:12)],Ranks=rankedResult[,6:9])
# Filtering: expr
set.seed(3141593)
rankedResult <- FilterByRank(revTopResults,c(9),1,threshold='q')    # topol increasing
ShortValidationChain(rankedResult[,1],rankedResult[,2:6],file="exprFilt_ranked",validation,revTopResults[,c(1,7:11)],Ranks=rankedResult[,7:11])

# Filtering: expr by pVal then rank Relevance
set.seed(3141593)
rankedResult <- FilterByRank(revTopResults,c(9),1,threshold=0.05)    # topol increasing
ShortValidationChain(rankedResult[,1],rankedResult[,2:6],file="exprFiltpVal_ranked",validation,revTopResults[,c(1,7:11)],Ranks=rankedResult[,7:11])


# Intersection: expr + topol reversed
# rankedResult <- FilterByRank(revTopResults,c(4,2),2)    # topol increasing
# ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="topExprIntersect",validation,results[,c(2,5:8)],Ranks=rankedResult[,6:9],threshold=100)
# Intersection: expr + topol reversed + relev
#rankedResult <- FilterByRank(revTopResults,c(6:8),2,threshold='q')    # topol increasing
#ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="allIntersect_ranked",validation,results[,c(2,9:12)],Ranks=rankedResult[,6:9])
# Intersection: expr + relev
set.seed(3141593)
rankedResult <- FilterByRank(revTopResults,c(8:9),2,threshold='q')    # topol increasing
ShortValidationChain(rankedResult[,1],rankedResult[,2:6],file="exprIntersect_ranked",validation,revTopResults[,c(1,7:11)],Ranks=rankedResult[,7:11])


# ShortValidationChain(rankedResult,file="",validation,results[,c(2,5:8)],topFunction=0.05)
# ShortValidationChain(rankedResult[,1],rankedResult[,2:5],file="ranked",validation,results[,c(2,5:8)],Ranks=rankedResult[,6:9],threshold=100)

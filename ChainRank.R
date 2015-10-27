library(igraph)
library(data.table)
library(plyr)
library(parallel)

options(stringsAsFactors = FALSE)

global.separator <- "|"
mc.cores <- 1

## Returns an weighted undirected igraph graph. The edge weights are the sum of the costs of the two vertices.
## The function transforms the rankings to costs: it inverts them, normalize them to 0-1 interval, then
## multiplies it by 5. Therefore the weight of the edges after summation are between 0.01-10.  edgelist - ...
## rankings - list of proteins and their rankings, the higher the better. Rows must be named by the the protein
## ids, same as used in edgelist.  weights - weight of the ranks. There is a possibility to sum up the
## different rankings, using a weight for each.(not tested). e.x. weight=c(1,0,0)

# Example: convert network (sNw) with BioXM names to entrez id (alias)
#   SNW=data.frame(mapName(sNw[,1],alias),mapName(sNw[,2],alias))
mapName <- function(names,map){
  row.names(map) <- map[,1]
  return(map[names,])
}

# rbind two matrices (m2 to m1), if column number differs it completes the smaller with 0s
SmartRbind <- function(m1, m2) {
  if (dim(m2)[2] == dim(m1)[2]) {
    m1 <- rbind(m1, m2)
  } else {
    
    if (dim(m2)[2] > dim(m1)[2]) {
      add <- matrix(0, dim(m1)[1], dim(m2)[2] - dim(m1)[2])
      m1 <- cbind(m1, add)
      m1 <- rbind(m1, m2)
    } else {
      add <- matrix(0, dim(m2)[1], dim(m1)[2] - dim(m2)[2])
      m2 <- cbind(m2, add)
      m1 <- rbind(m1, m2)
    }
  }
  if (length(m1[rowSums(m1) > 0, ]) == 0) {
    return(matrix(0, 1, 1))
  } else {
    if (is.vector(m1[rowSums(m1) > 0, ])) {
      return(as.matrix(t(m1[rowSums(m1) > 0, ])))
    } else {
      return(as.matrix(m1[rowSums(m1) > 0, ]))
    }
  }
}

# NULL test
isOk <- function(x){
  if(length(x) == 0){
    return(FALSE)
  }
  else return(TRUE)
}

## normalize vector or matrix (globally, not row or column)
normalize <- function(x) {
  if((max(x,na.rm=TRUE) - min(x, na.rm=TRUE))!=0){
    (x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x, na.rm=TRUE))
  } else {
    x - x
  } 
}

RemoveSemicol <- function(string) {
  string = gsub("\\;{2,}", ";", string)  #removes multiple semicolon from the string
  string = gsub("^\\;?", "", string)  #removes semicolon from the beginning of the string
  string = gsub("\\;?$", "", string)  #removes semicolon from the end of the string
  return(string)
}

RemoveGlobalSeparator <- function(string) {
  regExp = paste("\\", global.separator, "{2,}", sep = "")
  string = gsub(regExp, global.separator, string)  #removes multiple semicolon from the string
  
  regExp = paste("^\\", global.separator, "?", sep = "")
  string = gsub(regExp, "", string)  #removes semicolon from the beginning of the string
  
  regExp = paste("\\", global.separator, "?$", sep = "")
  string = gsub(regExp, "", string)  #removes semicolon from the end of the string
  return(string)
}

# Returns a vector containing the chain scores of the attrName attribute of the network
GetChainScores <- function(graph, chains,attrName) {
  chainScores <- c(0)
  scores <- get.vertex.attribute(graph, attrName, index=V(graph))
  for (i in 1:dim(chains)[1]) {
    chain <- chains[i,]
    trueLength <- sum(chain!=0)
    chainScores[i] <- sum(scores[chain]) / trueLength
  }
  chainScores <- data.frame(chainScores)
  names(chainScores) <- attrName
  return(chainScores)
}

# Computes scores and p-values of chains. 
# Args:
#   graph:
#   sores:
#   paths:
#   graph.rand:
#   paths.rand:
#   RetRank: Boolean variable. True means that scores are converted to ranks and returned. Otherwise raw scores are returned.
#
# Returns:
#   Data.frame containing p-values and the score or ranks of the chains.
#
GetChainScoresPval <- function(graph,scores,paths,graph.rand,paths.rand,RetRanks=FALSE){
  ## Add score attributes to igraph object
  SetScores <- function(graph,scores){
    name <- get.vertex.attribute(graph, "name")
    scores.names <- names(scores)[-1]
    scores.ordered <- mapName(name,scores)
    scores.ordered <- scores.ordered[!is.na(scores.ordered[,1]),-1]
    for(i in 1:length(scores.names)){
      graph <- set.vertex.attribute(graph,scores.names[i],,scores.ordered[,i])
    }
    return(graph)
  }
  
  ## Calculate chain scores of network
  CalcChainScores <- function(graph,paths,score.names){
    CScores <- NULL
    for (i in score.names) {
      if(!isOk(CScores)){
        CScores <- data.frame(GetChainScores(graph, paths,i))
      }else{
        CScores <- data.frame(CScores, GetChainScores(graph, paths,i))  
      }
    }
    return(CScores)
  }
  
  ## Ranking scores
  GetChainRanks <- function(Scores){
    Ranks <- data.frame(rank(Scores[,1]))
    for(i in 2:ncol(Scores[,])){
      score <- (Scores[,i] - max(Scores[,i]))*-1  #invert score
      rank <- rank(score)
      Ranks <- data.frame(Ranks,rank)
    }
    return(Ranks)
  }
  
  ## Create graph with  randomized scores, GloSS
  GetRandomGraph <- function(graph,scores){
    scores.rand <- data.frame(scores[,1]) #names
    for(i in 2:ncol(scores))
      scores.rand <- data.frame(scores.rand,sample(scores[,i])) #random scores
    names(scores.rand) <- names(scores)
    graph.rand <- SetScores(graph,scores.rand)
    return(graph.rand)
  }
  
  ## Calculate p-value Scott et al. 2011
  GetPval <- function(CScores.orig, CScores.rand){
    pVal <- NULL
    CScores.pval <- NULL
    for(i in 1:length(names(CScores.orig))){
      ## calculate p-value: % of scores in null model > than jth score
      for(j in 1:nrow(CScores.orig)){
        pVal[j] <- sum(CScores.orig[j,i] <= CScores.rand[,i]) / nrow(CScores.rand)
      }
      
      if(!isOk(CScores.pval)){
        CScores.pval <- data.frame(pVal)
      }else{
        CScores.pval <- data.frame(CScores.pval,pVal)
      }
    }
    names(CScores.pval) <- names(CScores.orig)
    return(CScores.pval)
  }
  
  
  ###
  score.names <- names(scores)[-1]
  
  ## Create graph with the ORIGINAL scores
  graph.orig <- SetScores(graph,scores)
  
  ## Calculate chain scores of original network
  CScores.orig <- CalcChainScores(graph.orig,paths,score.names)
  
  if(RetRanks){
    print("Calculating ranks...")
    CScores.ranks <- GetChainRanks(CScores.orig)
    names(CScores.ranks) <- paste(names(scores[,-1]),".rank",sep="")
  } else {
    CScores.ranks <- CScores.orig
    names(CScores.ranks) <- names(scores[,-1])
  }
    
  ## Calculate chain scores of random network
#   graph.rand <- SetScores(graph.rand,scores)
  graph.rand <- GetRandomGraph(graph.rand,scores)
  CScores.rand <- CalcChainScores(graph.rand,paths.rand,score.names)

  ## Calculate p-value
  CScores.pval <- GetPval(CScores.orig,CScores.rand)

  #naming
  
  names(CScores.pval) <- paste(names(scores[,-1]),".pval",sep="")

  return(data.frame(CScores.pval,CScores.ranks))
}

SynDFS <- function(graph, start, end, maxDepth, includeStart = FALSE) {
  
  GetLength <- function(graph, path) {
    Length <- c(0)
    for (i in 1:dim(path)[1]) {
      l = 0
      for (j in 1:dim(path)[2]) {
        if (path[i, j] == 0) {
          break
        }
        l = l + 1
      }
      Length[i] = l
    }
    return(Length)
  }
    
  # Recursive depth First Search algorithm
  RecDFS <- function(node, parents, graph, root, end, depth, maxDepth) {
      if (length(parents) == 1) {
          end = end[end != node]
      }
      maxLength = 1
      pathway = matrix(0, 1, maxLength)
      child = neighbors(graph, node)
      child = child[child != root & !(child %in% parents) & !(child %in% start)]	# simple path criteria: remove start nodes, already visited nodes 
      if (length(child) > 0) {
          for (i in 1:length(child)) {
              path = "search"
              
              # dead end
              if (depth == maxDepth | length(child[i]) < 1) {
                path = as.matrix(0)
              }
              
              # target node
              if (child[i] %in% end) {
                path = as.matrix(t(c(node, child[i])))
              }
              
              # not target or dead end
              if (path == "search") {
                tail = RecDFS(child[i], append(parents, node), graph, root, end, depth + 1, maxDepth)
                # cutting the rows where the sum is 0, therefore maxDepth or dead end reached tail =
                # tail[(1:dim(tail)[1])[rowSums(tail)!=0],]
                if (rowSums(tail) > 0) {
                  tail = tail[rowSums(tail) > 0, ]
                }
                if (is.vector(tail)) {
                  tail = as.matrix(t(tail))
                }
                nodes = matrix(node, dim(tail)[1], 1)
                path = cbind(nodes, tail)
              }
              
              pathway <- SmartRbind(pathway, path)
              maxLength <- dim(pathway)[2]
              
              # leave just those rows which contains end vertex exclusion of paths going through more than 1 end vertex is
              # done here
              containsEnd = NULL
              for (i in 1:dim(pathway)[1]) {
                row = pathway[i, ]
                row = row[which(row != 0)]
                containsEnd[i] = row[length(row)] %in% end
              }
              if (length(pathway[containsEnd, ]) == 0) {
                pathway = matrix(0, 1, maxLength)
              } else {
                if (is.vector(pathway[containsEnd, ])) {
                  pathway = as.matrix(t(pathway[containsEnd, ]))
                } else {
                  pathway = pathway[containsEnd, ]
                }
              }
          }
      }
      colnames(pathway) <- paste0("s",1:ncol(pathway))
      results <- as.matrix(pathway)
      # results <-list(as.matrix(pathway),scores)
      return(results)
  }
    
  maxDepth = maxDepth - 2  # depth correction. Without it start and end points are not part of depth.
  Paths <- matrix(0, 1, 1)
  Scores <- numeric(0)
  Lengths <- numeric(0)
  
  # run search algorithm for every start element
  # parallelized
  Paths.list <- mclapply(start, function(start, graph, end, maxDepth){
    RecDFS(start, 0, graph, start, end, 0, maxDepth)
    },graph=graph, end=end, maxDepth=maxDepth,mc.cores=mc.cores)
  
  Paths <- rbind.fill.matrix(Paths.list)
  Paths[is.na(Paths)] <- 0
  Paths <- unique(Paths)
  Lengths <- GetLength(graph, Paths)
  
  return(list(Paths, Lengths))
}

GenerateOutput <- function(chains, length, scores, name) {
    name <- as.vector(name)
    # create the proper chain representation
    start <- name[chains[, 1]]
    end <- c(0)
    chainElements <- c(0)
    for (i in 1:dim(chains)[1]) {
        e <- name[chains[i, 2]]
        ce <- ""
        for (j in 3:dim(chains)[2]) {
            if (chains[i, j] == 0) {
                break
            }
            ce <- paste(ce, e, sep = global.separator)
            e <- name[chains[i, j]]
        }
        chainElements[i] = RemoveGlobalSeparator(ce)
        end[i] <- e
    }
    
    # end <- name[end] merge all into a data frame
    data <- data.frame(start, chainElements, end, length, scores)
    # name columns
    names(data) <- c("Start protein", "Chain", "End protein", "Length", names(scores))
    
    return(data)
}

## Run Chain Search 
# Args:
#   Network: input network, edge list of IDs, 2 x n matrix
#   NetworkScore: score(s) of the nodes, n x m matrix where 1st col is ID
#   Candidates: ID(s) of candidates
#   Targets: ID(s) of candidates
#   maxDepth: maximal depth of the depth limited search
#   file: output file name
#   RetRanks - Boolean, if true results are returned as ranks instead of scores.
# 
# Returns:
#   A data frame with the Chains and their chain scores
RunChainSearch <- function(Network, NetworkScores, Candidates, Targets, maxDepth = 5,
    file=NULL, RetRanks = FALSE) {
    result = tryCatch({
        
        startList <- Candidates
        endList <- Targets
        sNw <- as.matrix(Network)
        sNw <- sNw[!(sNw[, 1] == "" | sNw[, 2] == ""), ]  # delete empty rows, deprecated, proper file should be provided
        
        ## Create graph and 
        myGraph <- graph.edgelist(sNw, directed = FALSE)
        name <- get.vertex.attribute(myGraph, "name")
        startPos <- which(name %in% startList)
        endPos <- which(name %in% endList)
        
        ## check score list - network size mismatch
        if (nrow(NetworkScores) < length(name)) {
          print(paste("Network - Protein list dimension errror. There are ", nrow(NetworkScores) - length(name), " more proteins in the network: "))
          print(name[!(name %in% NetworkScores[,1])])
          stop("Dimension mismatch")
        }
        
        print(date())
        print("Path search is starting...")
        system.time(Paths <- SynDFS(myGraph, startPos, endPos, maxDepth))
        if (Paths[[1]] == 0) {
            stop(paste("No paths were found with length less than ", maxDepth))
        } else {
            print("Paths were found")
        }
        
        
        ## Retrieve p-values
        ## Scott et al. 2005 type p-val calc.
        ## shuffle edges of the graph (keeping degree of nodes intact), shuffle weights
        print(date())
        print("Random network generation")
        myGraph.rand <- degree.sequence.game(degree(myGraph),method="vl")
        myGraph.rand <- set.vertex.attribute(myGraph.rand,"name",value=get.vertex.attribute(myGraph,"name"))
        print("Random network path search is starting...")
#         system.time(Paths.rand <- SynDFS(myGraph.rand, startPos, endPos, maxDepth))
        Paths.rand <- Paths
        print(date())
        print("Calculating p-values...")
        system.time(Scores <- GetChainScoresPval(myGraph,NetworkScores,Paths[[1]],myGraph.rand,Paths.rand[[1]],RetRanks=RetRanks))
        
        ## Output
        print(date())
        print("Generating output...")

        output <- GenerateOutput(Paths[[1]], Paths[[2]], Scores, name)
        
        if(isOk(file))
          write.table(output, file = file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
        
        print("Chain search algorithm finished!")
#         return(list(output, Paths[[1]], Paths[[2]], name))
        return(output)
        
    }, error = function(err) {
        print(paste("ERROR in execution :", err))
        return(NULL)
    })
    return(result)
}
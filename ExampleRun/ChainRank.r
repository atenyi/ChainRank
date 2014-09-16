library(igraph)
library(data.table)

global.separator = "|"
global.pathwSep = "||"

## Returns an weighted undirected igraph graph. The edge weights are the sum of the costs of the two vertices.
## The function transforms the rankings to costs: it inverts them, normalize them to 0-1 interval, then
## multiplies it by 5. Therefore the weight of the edges after summation are between 0.01-10.  edgelist - ...
## rankings - list of proteins and their rankings, the higher the better. Rows must be named by the the protein
## ids, same as used in edgelist.  weights - weight of the ranks. There is a possibility to sum up the
## different rankings, using a weight for each.(not tested). e.x. weight=c(1,0,0)
WeigthMx <- function(edgelist, rankings, weights) {
    # normalize function (should be removed, already in hasznosScriptek.r)
    normalize <- function(x) {
        if ((max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) != 0) {
            (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        } else {
            x - x
        }
    }
    
    rankNorm <- rankings %*% weights
    
    Cost = rowSums(rankNorm)
    
    # create graph from edges
    myGraph <- graph.edgelist(edgelist, directed = FALSE)
    
    # adjecency matrix
    adjMx <- get.adjacency(myGraph, "both")
    for (i in 1:dim(edgelist)[1]) {
        edgeCost = (Cost[edgelist[i, 1]] + Cost[edgelist[i, 2]])/2
        if (edgeCost == 0) {
            edgeCost = 0.01
        }
        tryCatch({
            adjMx[edgelist[i, 1], edgelist[i, 2]] = edgeCost
            adjMx[edgelist[i, 2], edgelist[i, 1]] = edgeCost
        }, error = function(e) {
            browser()
        })
    }
    # add weights to the graph
    myGraph <- graph.adjacency(adjMx, "undirected", weighted = TRUE)
    return(myGraph)
}

#
# Example: convert network (sNw) with BioXM names to entrez id (alias)
#   SNW=data.frame(mapName(sNw[,1],alias),mapName(sNw[,2],alias))
mapName <- function(names,map){
  row.names(map) <- map[,1]
  return(map[names,2])
}

# NULL test
isOk <- function(x){
  if(length(x) == 0){
    return(FALSE)
  }
  else return(TRUE)
}

# select rows in a matrix that contains element
rowContainsElement <- function(element,matrix){
  rows = matrix[unique(which(matrix == element,arr.ind=TRUE)[,1]),]
  return(rows)
}


# select columns in a matrix that contains element
columnContainsElement <- function(element,matrix){
  columns = matrix[,unique(which(matrix == element,arr.ind=TRUE)[,2])]
  return(columns)
}

toMatrix <- function(v,mode=NULL){
  if(!isOk(mode)){
    if(is.vector(v)){
      v=as.matrix(t(v))
    }
  }else{
    if(is.vector(v)){
      v=as.matrix(v)
    }
  }
  return(v)
}

getEnd <- function(vector){
  end = 0
  len = length(vector)
  while(end == 0){
    end = vector[len]
    len = len -1 
  }
  return(end)
}

getMxDiff <- function(m1,m2){
  
  l1=mclapply(seq_len(nrow(m1)), function(i) m1[i,])
  l2=mclapply(seq_len(nrow(m2)), function(i) m2[i,])
  
  getNotContainsEnd <- function(row,end){
    chain=row[which(row!=0)] 
    if(sum(end %in% chain[1:(length(chain)-1)])==0)
    {     
      return(row)
    }
  }
  
  lc=mclapply(l1,getNotContainsEnd,end=endPos)
  
  lcc=lc[!sapply(lc, is.null)]
  
  wow=do.call(rbind,lapply(lcc,matrix,ncol=ncol(m1),byrow=TRUE))
  
  W <- matrix(unlist(wow), ncol=ncol(wow), 
              dimnames=list(NULL, colnames(wow)))
  
  diff=na.omit(m2[W,which=TRUE])
}


#normalize vector or matrix function
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

GetScores <- function(graph, chains) {
  GetEdgeWeight <- function(graph, v1, v2) {
    if (length(v1) == 0 | length(v2) == 0) {
      return(0)
    } else {
      if (v1 == 0 | v2 == 0) {
        return(0)
      } else {
        ei <- get.edge.ids(graph, c(v1, v2))
        return(E(graph)$weight[ei])
      }
    }
  }
  
  scores <- c(0)
  for (i in 1:dim(chains)[1]) {
    score = 0
    g = 0
    for (j in 1:dim(chains)[2]) {
      if (chains[i, j] == 0) {
        break
      }
      score = score + GetEdgeWeight(graph, chains[i, j], chains[i, j - 1])
      g = g + 1
    }
    scores[i] = score/(g - 1)
  }
  # scores = scores
  return(scores)
}

## pathwList - matrix of protein names as line names and connected pathways (n x 1 matrix) chain - matrix of
## protein vectors, lines represents chains between start and end protein, proteins are represented by their
## number returns a vector with the protein chains and the connected pathways
GetPathways <- function(pathwList, chain) {
  # removes unnecesarry semicolons.
  pathways = c(0)
  for (i in 1:dim(chain)[1]) {
    pathways[i] = paste(pathwList[chain[i, ], 1], collapse = global.pathwSep)
    pathways[i] = RemoveSemicol(pathways[i])
  }
  return(pathways)
}

SynDFS <- function(graph, start, end, maxDepth, includeStart = FALSE) {
  
  GetLength <- function(graph, pathway) {
    Length <- c(0)
    for (i in 1:dim(pathway)[1]) {
      l = 0
      for (j in 1:dim(pathway)[2]) {
        if (pathway[i, j] == 0) {
          break
        }
        l = l + 1
      }
      Length[i] = l
    }
    return(Length)
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
    
    # Recursive depth First Search algorithm
    RecDFS <- function(node, parents, graph, root, start, end, depth, maxDepth) {
        if (length(parents) == 1) {
            end = end[end != node]
        }
        maxLength = 1
        pathway = matrix(0, 1, maxLength)
        child = neighbors(graph, node)
        child = child[child != root & !(child %in% parents) & !(child %in% start)]
        # cat('- Node:', node, ' | Child:',child,'\n')
        
        if (length(child) > 0) {
            for (i in 1:length(child)) {
                # cat('-- Node:', node, ' | FOR Child:',child[i],'\n')
                path = "search"
                
                # dead end
                if (depth == maxDepth | length(child[i]) < 1) {
                  path = as.matrix(0)
                }
                
                # target node
                if (child[i] %in% end) {
                  path = as.matrix(t(c(node, end[end == child[i]])))
                }
                
                # not target or dead end
                if (path == "search") {
                  tail = RecDFS(child[i], append(parents, node), graph, root, start, end, depth + 1, maxDepth)
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
        results <- as.matrix(pathway)
        # results <-list(as.matrix(pathway),scores)
        return(results)
    }
    
    maxDepth = maxDepth - 2  # depth correction. Without it start and end points are not part of depth.
    Paths <- matrix(0, 1, 1)
    Scores <- numeric(0)
    Lengths <- numeric(0)
    if (!includeStart) {
        inclStart <- start
    } else {
        inclStart <- c()
    }
    for (i in 1:length(start)) {
        resPaths <- RecDFS(start[i], 0, graph, start[i], inclStart, end, 0, maxDepth)
        
        Paths <- SmartRbind(Paths, resPaths)
        
        cat("Paths: ", i/length(start) * 100, "%,\n")
    }
    Paths <- unique(Paths)
    Lengths <- GetLength(graph, Paths)
    
    return(list(Paths, Lengths))
}

GenerateOutput <- function(chains, length, scores, pathways, name) {
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
    data <- data.frame(start, chainElements, end, length, scores, pathways)
    # name columns
    names(data) <- c("Start protein", "Chain", "End protein", "Length", names(scores), "Pathways")
    
    return(data)
}

## RunChainSearch @validation - vector, containing element name used for validation @return - output, Paths by
## numbering, Paths length
RunChainSearch <- function(synergyNetwork, synergyNetworkProt, synergyNetworkPath, Candidates, Targets, nwDim, 
    pathDim, scoreDim, file, extendedStat = FALSE, alias = 0, maxDepth = 5, validation = NULL) {
    result = tryCatch({
        
        startList <- unlist(strsplit(Candidates, split = ";"))  #in case more proteins are expressed by one gene
        endList <- unlist(strsplit(Targets, split = ";"))
        sNw <- as.matrix(synergyNetwork[, nwDim])
        sNw <- sNw[!(sNw[, 1] == "" | sNw[, 2] == ""), ]  # delete empty rows
        myGraph <- graph.edgelist(sNw, directed = FALSE)
        name <- get.vertex.attribute(myGraph, "name")
        startPos <- which(name %in% startList)
        endPos <- which(name %in% endList)
        
        # Pathways
        
        sNwPath <- as.matrix(synergyNetworkPath[, pathDim])
        attributes(sNwPath)$dimnames[[1]] = synergyNetworkPath[, 1]
        # sNwPath <- as.matrix(sNwPath[name %in% row.names(sNwPath),])
        sNwPath <- as.matrix(sNwPath[row.names(sNwPath) %in% name, ])
        if (nrow(sNwPath) < length(name)) {
            print(paste("Network - Protein list dimension errror. There are ", nrow(sNwPath) - length(name), " more proteins in the network: "))
            print(name[!(name %in% row.names(sNwPath))])
            stop("Dimension mismath")
        }
        
        print("Path search is starting...")
        Paths <- SynDFS(myGraph, startPos, endPos, maxDepth)
        if (Paths[[1]] == 0) {
            stop(paste("No paths were found with length less than ", maxDepth))
        } else {
            print("Paths were found")
        }
        
        
        # Pathways
        sNwPath <- as.matrix(synergyNetworkPath[, pathDim])
        attributes(sNwPath)$dimnames[[1]] = synergyNetworkPath[, 1]
        # sNwPath <- as.matrix(sNwPath[name %in% row.names(sNwPath),])
        sNwPath <- as.matrix(sNwPath[row.names(sNwPath) %in% name, ])
        Pathways <- GetPathways(sNwPath, Paths[[1]])
        
        # Scores
        print("Calculating ranks...")
        Ranking <- as.matrix(synergyNetworkProt[, scoreDim])
        Ranking <- apply(Ranking, 2, as.numeric)
        row.names(Ranking) = synergyNetworkProt[, 1]
        # Ranking <- as.matrix(Ranking[name,])
        Ranking <- as.matrix(Ranking[row.names(Ranking) %in% name, ])
        weights <- vector("integer", 3)
        weights[1] = 1
        myGraph <- WeigthMx(sNw, Ranking, weights)  #create graph from edges
        
        Scores <- GetScores(myGraph, Paths[[1]])
        # Scores <- GetScores(sNw,Ranking,weights,Paths[[1]])
        if (dim(Ranking)[2] > 1) {
            for (i in 2:dim(Ranking)[2]) {
                weights <- vector("integer", dim(Ranking)[2])
                weights[i] = 1
                myGraph <- WeigthMx(sNw, Ranking, weights)  #create graph from edges
                Scores <- cbind(Scores, GetScores(myGraph, Paths[[1]]))
            }
        }
        print("Generating output...")
        Scores <- as.data.frame(Scores)
        names(Scores) <- names(synergyNetworkProt)[scoreDim]
        if (length(alias) > 1) {
            name <- alias[name, 2]
        }
        
        output <- GenerateOutput(Paths[[1]], Paths[[2]], Scores, Pathways, name)
        
        write.table(output, file = paste(file, "results.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
        
        if (extendedStat) {
            GetSummaryStatExtended(output, Paths[[1]], file, name, validation)
        }
        print("Chain search algorithm finished!")
        return(list(output, Paths[[1]], Paths[[2]], name))
        
    }, error = function(err) {
        print(paste("ERROR in execution :", err))
        return(NULL)
    })
    return(result)
}

ShortValidationChain <- function(rankAndName, file, validation, pathways = NULL, valData, weights = NULL) {
    
    # RankPathwayChain <- function(var) {
        # # Pathways stat
        # pw <- unlist(strsplit(as.character(var), split = ";"))
        # pw.table <- table(pw)
        # result <- as.data.frame(pw.table)
        # return(result[order(result[, 2], decreasing = TRUE), ])
    # }
    
    ExhValidation <- function(sudt, rankI, file, validation, pathway, random) {
        avgImprovement <- 0  # average Improvement
        N <- length(c(8:12))
        startLength <- 8
        endLength <- 12
        
        for (i in c(startLength:endLength)) {
            topGrp10 <- unique(sudt[, 1])[1:i]
            topRank10 <- min(sudt[sudt[, 1] %in% topGrp10, rankI])
            top10 <- sudt[sudt[, rankI] >= topRank10, 2]
            t10 <- as.data.frame(table(top10))
            
            # validation
            if (!is.null(validation)) {
                top10elem <- unique(top10)
                top10v <- table(c(as.character(top10elem), as.character(validation)))
                top10v <- as.data.frame(top10v - 1)
                top10v <- top10v[top10v[, 2] > 0, ]
                t10v <- t10[which(t10[, 1] %in% top10v[, 1]), ]
                
                # smta <- merge(smta,top10v,by.x = 'Name',by.y='Var1',all=TRUE) smta[is.na(smta)]=0
                # names(smta)[length(names(smta))] <- (paste(rank,' validation'))
                
                TP <- sum(t10v[, 2])
                FP <- sum(t10[, 2]) - TP
                Pr <- TP/(TP + FP)
                
                nTP <- nrow(top10v)
                nFP <- length(top10elem) - nTP
                nPr <- nTP/(nTP + nFP)
                
                # random set validation
                vTP <- random$vTP[i - (startLength - 1)]
                vFP <- random$vFP[i - (startLength - 1)]
                vnTP <- random$vnTP[i - (startLength - 1)]
                vnFP <- random$vnFP[i - (startLength - 1)]
                
                vPr <- vTP/(vTP + vFP)
                vnPr <- vnTP/(vnTP + vnFP)
                
                # calculating improvement
                imprGross <- Pr/vPr
                imprNet <- nPr/vnPr
                avgImprovement <- avgImprovement + imprNet
                
                # printing out results
                out = data.frame(Name = c("Gross", "Gross random", "Net", "Net random"), TruePositive = c(TP, 
                  vTP, nTP, vnTP), FalsePositive = c(FP, vFP, nFP, vnFP), Precision = c(Pr, vPr, nPr, vnPr), Improvement = c("", 
                  imprGross, "", imprNet))
                cat(paste("\nTest ", i, "\n"), file = paste(file, "validation.txt"), append = TRUE)
                write.table(out, file = paste(file, "validation.txt"), append = TRUE, row.names = FALSE, quote = FALSE, 
                  sep = "\t", col.names = TRUE)
                
                cat(paste("\n", "Test ", i, " ", names(sudt)[rankI], "\t AvgFreq = ", mean(t10v[, 2]), "\n"), 
                  file = paste(file, "topProteins.txt"), append = TRUE)
                write.table(data.frame(Proteins = t10v), file = paste(file, "topProteins.txt"), append = TRUE, 
                  row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
                # cat(paste('\n RANDOM\t AvgFreq = ',mean(v10v[,2]),'\n'),file=paste(file,'topProteins.txt'),append=TRUE)
                # write.table(data.frame(Proteins=v10v),file=paste(file,'topProteins.txt'),append=TRUE,row.names=FALSE,
                # quote=FALSE, sep='\t',col.names=TRUE)
                
                if (i == 10) {
                  # write.table(v10, file=paste(file,'randomProteins.txt'),append=FALSE,row.names=FALSE, quote=FALSE,
                  # sep='\t',col.names=TRUE)
                }
                
                print("Validation done")
            }
            
            # pathways on net top10
            # if (i == 10 & !is.null(pathways)) {
                # top <- as.data.frame(unique(top10))
                # top <- merge(top, pathways, by.x = names(top)[1], by.y = names(pathways)[1])
                # write.table(RankPathwayChain(top), file = paste(file, "pathway", names(sudt)[rankI], ".txt"), 
                  # row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
            # }
        }
        cat(paste("\nAvg net improvement \t", avgImprovement/N, "\n"), file = paste(file, "validation.txt"), append = TRUE)
    }
    
    ## Output table containing the frequency of proteins that are in the top quantile chains (default is the last rank)
    ## topFunction 	- 'q' -> top quantile is selected
	##				- number -> top number is selected
    RankFreqChain <- function(rankAndName, file, rankI = 4, topFunction) {
        # rankI <- 4 top <- quantile()
        sudt = rankAndName[order(rankAndName[, rankI], decreasing = TRUE), ]
        if (isOk(topFunction) & topFunction == "q") {
            topRank = quantile(sudt[, rankI])[4]
        }
        if (isOk(topFunction) & is.numeric(topFunction)) {
            sudt.ordered <- sudt[order(sudt[, rankI], decreasing = TRUE), ]
            topRank <- sudt.ordered[topFunction, rankI]
        }
        
        topquant <- sudt[sudt[, rankI] >= topRank, ]
        # separate csv format to list
        apl = lapply(topquant, function(o) {
            strsplit(as.character(o), paste("\\", global.separator, sep = ""))
        })
        # 
        apf = data.frame(index = rep(seq_along(apl[[1]]), lapply(apl[[1]], length)), Chain = unlist(apl[[1]]), 
            Topology = rep(as.numeric(unlist(apl[[2]])), lapply(apl[[1]], length)), Relevance = rep(as.numeric(unlist(apl[[3]])), 
                lapply(apl[[1]], length)), Expression = rep(as.numeric(unlist(apl[[4]])), lapply(apl[[1]], length)))
        # network representation of chains
        chain <- apl[[1]]
        listToNetwork <- function(list) {
            network <- data.frame(v1 = list[1:(length(list) - 1)], v2 = list[2:length(list)])
        }
        chain.network <- lapply(chain, listToNetwork)  # network of top chains in list of chain networks format
        chain.network <- do.call(rbind.data.frame, chain.network)  # network of top chains in bulk data frame format
        
        tq = as.data.frame(table(apf[, 2]))
        write.table(tq, file = paste(file, "top", names(sudt[rankI]), "Proteins.txt"), append = FALSE, row.names = FALSE, 
            quote = FALSE, sep = "\t", col.names = TRUE)
        write.table(chain.network, file = paste(file, "top", names(sudt[rankI]), "Network.txt"), append = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
        write.table(rankAndName, file = paste(file, "top", names(sudt[rankI]), "Results.txt"), append = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
    }
    
    ## Returns the chain and the row wise sum of the ranks
    MakeCombinedRank <- function(rankAndName, weights) {
        ## Add Combined rank combines the three ranking into one
        ranks.norm <- apply(rankAndName[, -1], 2, normalize)
        if (!isOk(weights)) {
            weights <- c(1, 1, 1)
        }
        ranks.rowsum <- ranks.norm %*% weights
        # ranks.rowsum <- rowSums(ranks.norm)
        return(data.frame(Chain = rankAndName[, 1], Combined = ranks.rowsum))
    }
    
    RandomRanking <- function(valData, validation) {
        n <- 100
        vTP <- NULL
        vFP <- NULL
        vnTP <- NULL
        vnFP <- NULL
        vTP_c <- NULL
        vFP_c <- NULL
        vnTP_c <- NULL
        vnFP_c <- NULL
        for (i in 8:12) {
            for (j in 1:n) {
                index <- unique(valData[, 1])
                valGrp10 <- sample(index, i)
                val10 <- valData[valData[, 1] %in% valGrp10, 2]
                
                val10elem <- unique(val10)
                val10v <- table(c(as.character(val10elem), as.character(validation)))
                val10v <- as.data.frame(val10v - 1)
                val10v <- val10v[val10v[, 2] > 0, ]
                v10 <- as.data.frame(table(val10))
                v10v <- v10[which(v10[, 1] %in% val10v[, 1]), ]
                
                vTP <- vTP + sum(v10v[, 2])
                vFP <- vFP + sum(v10[, 2]) - sum(v10v[, 2])
                vTN <- length(validation) - sum(v10v[, 2])
                
                vnTP <- vnTP + nrow(v10v)
                vnFP <- vnFP + length(val10elem) - nrow(v10v)
                
                vTP_c <- c(vTP_c, vTP/n)
                vFP_c <- c(vFP_c, vFP/n)
                vnTP_c <- c(vnTP_c, vnTP/n)
                vnFP_c <- c(vnFP_c, vnFP/n)
            }
        }
        
        return(data.frame(vTP = vTP_c, vFP = vFP_c, vnTP = vnTP_c, vnFP = vnFP_c))
    }
    
    # Preprocess results
    combinedRank <- MakeCombinedRank(rankAndName, weights)  # create combined rank
    rankAndName.combined <- merge(rankAndName, combinedRank, by = "Chain")
    apl = lapply(rankAndName.combined, function(o) {
        strsplit(as.character(o), paste("\\", global.separator, sep = ""))
    })
    
    # Map chain id and ranks to proteins
    apf = data.frame(index = rep(seq_along(apl[[1]]), lapply(apl[[1]], length)), Chain = unlist(apl[[1]]), Topology = rep(as.numeric(unlist(apl[[2]])), 
        lapply(apl[[1]], length)), Relevance = rep(as.numeric(unlist(apl[[3]])), lapply(apl[[1]], length)), Expression = rep(as.numeric(unlist(apl[[4]])), 
        lapply(apl[[1]], length)), Combined = rep(as.numeric(unlist(apl[[5]])), lapply(apl[[1]], length)))
    
    aplVal = lapply(valData, function(o) {
        strsplit(as.character(o), paste("\\", global.separator, sep = ""))
    })
    apfVal = data.frame(index = rep(seq_along(aplVal[[1]]), lapply(aplVal[[1]], length)), Chain = unlist(aplVal[[1]]), 
        Topology = rep(as.numeric(unlist(aplVal[[2]])), lapply(aplVal[[1]], length)), Relevance = rep(as.numeric(unlist(aplVal[[3]])), 
            lapply(aplVal[[1]], length)), Expression = rep(as.numeric(unlist(aplVal[[4]])), lapply(aplVal[[1]], 
            length)))
    
    names(apf)[-1] <- names(rankAndName.combined)
    
    # Validation
    append = FALSE
    RandomRanking <- RandomRanking(apfVal, validation)
    
    for (i in 3:length(apf)) {
        if (i != 3) {
            append = TRUE
        }
        sudt = apf[order(apf[, i], decreasing = TRUE), ]
        cat(paste("\n\n", names(apf)[i], " \n"), file = paste(file, "topProteins.txt"), append = append)
        cat(paste("\n\n", names(apf)[i], " \n"), file = paste(file, "validation.txt"), append = append)
        ExhValidation(sudt, i, file, validation, pathways, RandomRanking)
        RankFreqChain(rankAndName.combined, file, i - 1, 10)  # output frequency table of top proteins
    }
}

# $data - data frame containing chain and connected ranks, consequently $order - using this variable the order
# of the ranking can be change from the default order $type - can take value {1,2}.  1 filtering ranking 2
# intersection ranking
FilterByRank <- function(data, order = NULL, type = NULL) {
    if (!isOk(order)) {
        order <- 2:ncol(data)
    }
    nRank <- length(order) + 1
    
    ranks <- data.frame(index = 1:nrow(data), data[, order])
    intersection <- ranks[, 1]
    for (i in 2:nRank) {
        quant <- quantile(ranks[, i])[4]
        ranks.quant <- ranks[ranks[, i] >= quant, ]
        
        ## Filtering ranking filter the top quantile using the 1st ranking criteria in the order, then using the
        ## results it rank again using the 2nd rank, etc.
        if (!isOk(type) | type == 1) {
            ranks <- ranks.quant
        }
        
        ## Intersection ranking filter using the intersection of the top quantile chains of all of the rankings
        if (isOk(type) & type == 2) {
            intersection <- intersect(intersection, ranks.quant[, 1])
            if (i == nRank) {
                ranks <- ranks[ranks[, 1] %in% intersection, ]
            }
        }
    }
    return(data[ranks[, 1], ])
} 

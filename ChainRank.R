library(igraph)
library(data.table)
library(plyr)
library(parallel)

global.separator <- "|"
global.pathwSep <- "||"
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
    chainScores[i] <- sum(scores[chains[i,]]) / trueLength
  }
  chainScores <- data.frame(chainScores)
  names(chainScores) <- attrName
  return(chainScores)
}

## Returns the p-values of the chains scores
GetChainScoresPval_Ideker <- function(graph,scores,paths,graph.rand,paths.rand){
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
  
  print("Calculating ranks...")
  CScores.ranks <- GetChainRanks(CScores.orig)
  
  ## Calculate chain scores of random network
#   graph.rand <- SetScores(graph.rand,scores)
  graph.rand <- GetRandomGraph(graph.rand,scores)
  CScores.rand <- CalcChainScores(graph.rand,paths.rand,score.names)

  ## Calculate p-value
  CScores.pval <- GetPval(CScores.orig,CScores.rand)

  #naming
  setnames(CScores.ranks,names(CScores.pval))
  setnames(CScores.pval,paste(names(CScores.pval),".pvals",sep=""))

  return(data.frame(CScores.pval,CScores.ranks))
}

## pathwList - matrix of protein names as line names and connected pathways (n x 1 matrix) chain - matrix of
## protein vectors, lines represents chains between start and end protein, proteins are represented by their
## number returns a vector with the protein chains and the connected pathways
GetPathways <- function(pathwList, chain) {
#   removes unnecesarry semicolons.
  pathways = c(0)
  for (i in 1:dim(chain)[1]) {
    pathways[i] = paste(pathwList[chain[i, ], 1], collapse = global.pathwSep)
    pathways[i] = RemoveSemicol(pathways[i])
  }
  return(pathways)
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

GenerateOutput <- function(chains, length, scores, name, pathways=NULL) {
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
    
    if(isOk(pathways)){          
    # end <- name[end] merge all into a data frame
    data <- data.frame(start, chainElements, end, length, scores, pathways)
    # name columns
    names(data) <- c("Start protein", "Chain", "End protein", "Length", names(scores), "Pathways")
    }else{
    # end <- name[end] merge all into a data frame
    data <- data.frame(start, chainElements, end, length, scores)
    # name columns
    names(data) <- c("Start protein", "Chain", "End protein", "Length", names(scores))
    }
    
    return(data)
}

## RunChainSearch 
## @return - output, Paths by numbering, Paths length
## @pval - 1: GloSS like p-value calculation
##        2: Scott et al. 2011 like p-value calculation
## @NetworkPath - possible future 
RunChainSearch <- function(Network, NetworkScores, Candidates, Targets, maxDepth = 5,
    file, alias=0, NetworkPath=NULL) {
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
        if (nrow(NetworkScores) != length(name)) {
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
        system.time(Paths.rand <- SynDFS(myGraph.rand, startPos, endPos, maxDepth))
        
        print(date())
        print("Calculating p-values...")
        system.time(Scores <- GetChainScoresPval_Ideker(myGraph,NetworkScores,Paths[[1]],myGraph.rand,Paths.rand[[1]]))
        
        ## Pathways
        if(isOk(NetworkPath)){
          sNwPath <- as.matrix(NetworkPath[, -1])
          attributes(sNwPath)$dimnames[[1]] = NetworkPath[, 1]
          sNwPath <- as.matrix(sNwPath[row.names(sNwPath) %in% name, ])
          Pathways <- GetPathways(sNwPath, Paths[[1]])
        }else{Pathways<-NULL}
        
        ## Output
        print(date())
        print("Generating output...")
        if (length(alias) > 1) {
            name <- alias[name, 2]
        }
        output <- GenerateOutput(Paths[[1]], Paths[[2]], Scores, name, Pathways)
        
        write.table(output, file = paste(file, "_results.txt",sep=""), row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
        
        print("Chain search algorithm finished!")
        return(list(output, Paths[[1]], Paths[[2]], name))
        
    }, error = function(err) {
        print(paste("ERROR in execution :", err))
        return(NULL)
    })
    return(result)
}

## weights - weight vector, if 0 is the input value then 
ShortValidationChain <- function(Chains,Pvals, goldStandard, valData, weights = NULL,Ranks=NULL,threshold=NULL,file="") {
    ## Splits the 1st column of data into separate rows by the separator argument.
    ## Corresponding indexes of data's rows are returned in the 1st column, can be used as groups
    SplitToGroups <- function(data,separator){
      apl = lapply(data, function(o) {
        strsplit(as.character(o), paste("\\", separator, sep = ""))
      })
      
      apf = data.frame(index = rep(seq_along(apl[[1]]), lapply(apl[[1]], length)), Chain = unlist(apl[[1]]))
      
      for(i in 2:length(apl)){
        apf <- data.frame(apf,rep(as.numeric(unlist(apl[[i]])),lapply(apl[[1]], length)))
      }
      
      setnames(apf,c("index",names(data)))
      return(apf)
    }
    
    ## topFunction   - 'q' -> top quantile is selected
    ##				      - integer -> top number is selected
    ##              - 0<floating<1 -> p-value cut off
    ExhValidation <- function(sudt, rankI, file, validation, valData) {
        if (isOk(topFunction) & topFunction == "q") {
          chainRepresenter<-sudt[!duplicated(sudt$index),]
          topRank <- quantile(chainRepresenter[,rankI])[4]
          topquant <- sudt[sudt[, rankI] >= topRank, ]
        }
        if (isOk(topFunction) & is.numeric(topFunction)) {
          if(0<=topFunction & topFunction<1){
            topRank <- topFunction
            topquant <- sudt[sudt[, rankI] <= topRank & !duplicated(sudt$index), ]
          }
          else{
            sudt.ordered <- sudt[order(sudt[, rankI], decreasing = TRUE), ]
            sudt.ordered.unique <- sudt.ordered[!duplicated(sudt.ordered$index),]
            topRank <- sudt.ordered.unique[topFunction,rankI] 
            topquant <- unique(sudt[sudt[, rankI] >= topRank & !duplicated(sudt$index), ])
          }
        }
        
        N <- nrow(topquant)
        startLength <- N
        endLength <- N
        avgImprovement <- 0  # average Improvement
        
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
                random <- RandomRanking(valData,validation,i)
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
    
    
    # Calculates cummulative precision and recall rates for each chain
    # $data - data frame containing information on chains. A list of protein with chain identifier and p-value.
    #         format: ChainIndex, Protein name, p-value
    CumStats <- function(data,valData,goldStandard,file,extRank=NULL){
      ##TODO investigate why so slow
      getCumStats <- function(pVal,data,goldStandard){
        names(data)[1:2] <- c("index","ID")
        
        top.elem <- data[data[,3] <= pVal,]
        top.chain.n <- sum(!duplicated(top.elem$index)) # Number of top chains
        top.elem.freq <- as.data.frame(table(top.elem$ID))  #top.elem$IDs might be factors -> top.gs can contain all the factors (with 0 frequency those that doesn't occure)
        top.gs <- top.elem.freq[which(top.elem.freq[, 1] %in% goldStandard), ]
        top.gs <- top.gs[top.gs[,2]>0,]
        
        # Precision from frequency values
        TP <- sum(top.gs[, 2])
        FP <- sum(top.elem.freq[, 2]) - TP
        Precision <- TP/(TP + FP)
        
        # Precision, Recall from unique occurances
        nTP <- nrow(top.gs)
        nFP <- sum(top.elem.freq[,2]>0) - nTP
        nFN <- length(goldStandard) - nTP
        nPrecision <- nTP/(nTP + nFP)
        nRecall <- nTP/length(goldStandard)
        
        if(ncol(data)>3){
          pVal <- data.frame(pVal,data[which(data[,3]==pVal)[1],4])
          # Distribution of chains
          Distrib <- 0
        }else{
          # Distribution of chains
          Distrib <- top.chain.n/length(unique(data$index))
        }
        
        out <- data.frame(pVal,top.chain.n,nTP,nFP,nFN,nPrecision,nRecall,TP, FP, Precision,Distrib)
        return(out)
      }
      
      getRandomRank <- function(data,valData){
        random.chains <- sample(valData[,1],nrow(data))
        random.data <- data.frame(random.chains,data[,2])
        protList <- SplitToGroups(random.data,global.separator)
      }

      if(isOk(extRank)){
        data <- data.frame(data,extRank)
        names(data)[3] <- "p-value"
      }
      
      protList <- SplitToGroups(data,global.separator)
      
      pVals <- sort(unique(protList[,3]))
      
      goldStandard.contained <- goldStandard[goldStandard %in% unique(protList[,2])]
      Out <- sapply(pVals,getCumStats,data=protList,goldStandard=goldStandard.contained)
      Out <- as.data.frame(t(Out))
      
      # Simulated results
      n <- 10
      Random.n <- 0
      for(i in 1:n){
        protList.random <- getRandomRank(data,valData)
        Random.mx <- sapply(pVals,getCumStats,data=protList.random,goldStandard=goldStandard.contained)
        Random.mx <- as.data.frame(t(Random.mx))[,c(-2,-11)]  # remove N chains, distribution columns
        Random.n <- Random.n + data.matrix(Random.mx)
      }
      Random <- as.data.frame(Random.n/n)

      names(Out) <-c(names(data)[-1],"N.chains","net.TP","net.FP","net.FN","net.Precision","net.Recall","gross.TP", "gross.FP", "gross.Precision","Distribution")
      names(Random) <-c("R.p-values","R.net.TP","R.net.FP","R.net.FN","R.net.Precision","R.net Recall","R.gross.TP", "R.gross.FP", "R.gross.Precision")
      output <- data.frame(lapply(data.frame(Out,Random[,-1]), as.numeric), stringsAsFactors=FALSE)
      
      # Calculate improvement
      grossImprovement <- output$gross.Precision/output$R.gross.Precision
      netImprovement <- output$net.Precision/output$R.net.Precision

      # Calculate F-measure
      Fmeasure <- 2*output[,5]*output[,6]/(output[,5]+output[,6])
      
      output <- data.frame(output,gross.Improvement=grossImprovement,net.Improvement=netImprovement,Fmeasure=Fmeasure)
      
      write.table(output, file = paste("CummulativeStats_",file,".txt",sep=""), append = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
      
      print("Cummulative statistics calculated")
    }

    ## Returns the chain and the row wise weighted sum of the ranks
    MakeCombinedRank <- function(rankAndName, weights) {
        ## Add Combined rank combines the three ranking into one
        ranks.norm <- apply(rankAndName[,2:4], 2, normalize)
        if (!isOk(weights)) {
            weights <- c(1, 1, 1)
        }
        ranks.rowsum <- (ranks.norm %*% weights)/length(weights)
        # ranks.rowsum <- rowSums(ranks.norm)
        return(data.frame(Chain = rankAndName[, 1], Combined = ranks.rowsum))
    }
    
    RandomRanking <- function(valData, validation, topN) {
          n <- 100
          vTP <- 0
          vFP <- 0
          vnTP <- 0
          vnFP <- 0
          vTP_c <- 0
          vFP_c <- 0
          vnTP_c <- 0
          vnFP_c <- 0
          for (j in 1:n) {
              index <- unique(valData[, 1])
              valGrp10 <- sample(index, topN)
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
          }
#           vTP_c <-  c(vTP_c,vTP/n)
#           vFP_c <- c(vFP_c, vFP/n)
#           vnTP_c <- c(vnTP_c, vnTP/n)
#           vnFP_c <- c(vnFP_c, vnFP/n)
        
        return(data.frame(vTP = vTP/n, vFP = vFP/n, vnTP = vnTP/n, vnFP = vnFP/n))
    }
    
    ####
    if(nrow(Pvals)==0){
      print(paste("No chains were recovered for ",file,sep=""))
      return(0)
    }
    if(!isOk(Ranks)){
      rankAndName <- data.frame(Chains,Pvals)
      extRank <- NULL
    }else{
      rankAndName <- data.frame(Chains,Ranks)
      extRank <- Pvals
    }
  
    # Combined rank 
#     combinedRank <- MakeCombinedRank(rankAndName, weights)  # create combined rank
#     rankAndName.combined <- merge(rankAndName, combinedRank, by=1)
    rankAndName.combined <- rankAndName    

    apf <- SplitToGroups(rankAndName.combined,global.separator)
    apfVal <- SplitToGroups(valData,global.separator)
    
    # Reduce gold standards to the ones that are contained in the network
    

    # Validation
    append = FALSE
    for (i in 2:ncol(rankAndName.combined)) {
#       for (i in 5:5) {
#         if (i != 3) {
#             append = TRUE
#         }
#         sudt = apf[order(apf[, i], decreasing = TRUE), ]
#         cat(paste("\n\n", names(apf)[i], " \n"), file = paste(file, "topProteins.txt"), append = append)  # print file headers
#         cat(paste("\n\n", names(apf)[i], " \n"), file = paste(file, "validation.txt"), append = append) # print file headers
#         ExhValidation(sudt, i, file, validation, apfVal,topFunction)
        CumStats(rankAndName.combined[,c(1,i)], valData[,c(1,i)],goldStandard,paste(names(rankAndName.combined)[i],file,sep="_"),as.data.frame(extRank[,i-1]))
    }
}

ExtendedResults <- function(rankAndName, file, weights = NULL,topFunction){
  RankFreqChain <- function(rankAndName, file, rankI = 4, topFunction) {
    # rankI <- 4 top <- quantile()
    sudt = rankAndName[order(rankAndName[, rankI], decreasing = TRUE), ]
    if (isOk(topFunction) & topFunction == "q") {
      topRank = quantile(sudt[, rankI])[4]
      topquant <- sudt[sudt[, rankI] >= topRank, ]
    }
    if (isOk(topFunction) & is.numeric(topFunction)) {
      if(0<topFunction & topFunction<1){
        topRank <- topFunction
        topquant <- sudt[sudt[, rankI] <= topRank, ]
      }
      else{
        sudt.ordered <- sudt[order(sudt[, rankI], decreasing = TRUE), ]
        topRank <- sudt.ordered[topFunction, rankI]  
        topquant <- sudt[sudt[, rankI] >= topRank, ]
      }
    }
    
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
  
  ## Returns the chain and the row wise weighted sum of the ranks
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
  
  # Preprocess results
#   combinedRank <- MakeCombinedRank(rankAndName, weights)  # create combined rank
#   rankAndName.combined <- merge(rankAndName, combinedRank, by = "Chain")
rankAndName.combined <- rankAndName
  for (i in 2:ncol(rankAndName.combined)) {
    ## Output table containing the frequency of proteins that are in the top quantile chains (default is the last rank)
    ## topFunction   - 'q' -> top quantile is selected
    ##				      - integer -> top number is selected
    ##              - 0<floating<1 -> p-value cut off
    RankFreqChain(rankAndName.combined, file, i, topFunction)  # output frequency table of top proteins
  }
}

# $data - data frame containing chain and connected ranks, consequently 
# $order - using this variable the order of the ranking can be change from the default order 
# $type - can take value {1,2}.  1 filtering ranking, 2 intersection ranking
FilterByRank <- function(data, order = NULL, type = NULL,threshold = 0.05) {
    if (!isOk(order)) {
        order <- 2:ncol(data)
    }
    nRank <- length(order) + 1
    
    ranks <- data.frame(index = 1:nrow(data), data[, order])
    intersection <- ranks[, 1]
    for (i in 2:nRank) {
      if(threshold=='q'){  
        quant <- quantile(ranks[, i])[2]
        ranks.quant <- ranks[ranks[, i] <= quant, ]
      }else{
        ranks.quant <- ranks[ranks[, i] <= threshold, ]
      }
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


ChainRank
=========

ChainRank, a chain search based method for prioritization and contextualization of biological subnetworks

Motivation: Advances in high throughput technologies and growth of biomedical knowledge has contributed to an expo-nential increase in associative data. These can be represented in the form of complex networks of biological associations which are suitable for systems analyses. However, these net-works usually lack both, context specificity in time and space as well as the distinctive borders usually assigned in the clas-sical pathway view of molecular events (e.g. signal transduc-tion). This complexity and high interconnectedness call for automated techniques that can identify smaller targeted sub-networks specific to a given research context (e.g. a disease scenario).
Results: Our method, named ChainRank, finds relevant sub-networks by identifying and scoring chains of interactions that link specific network components. Scores can be generated from integrating multiple general and context specific measures. The performance of the novel ChainRank method was evaluated on recreating specific signalling pathways from a human protein interaction network. The analysis showed that ChainRank can identify main mediators of context specific molecular signalling.

=========

The ChainRank.R file contains the ChainRank algorithm. An example run is also available to test the algorithm. The ExampleRun folder contains the necessary data files and the MuscleSpecCase_run.R file which contains the script for the example. The scenario aims to recreate the IGF-Akt pathway explained in (Schiaffino and Mammucari, 2011) and reported in (Tenyi et al., 2015).

  Usage:

   RunChainSearch(Network, NetworkScores, Candidates, Targets, maxDepth = 5,file=NULL, RetRanks = FALSE)
  
  Args:
  
    Network: input network, edge list of IDs, 2 x n matrix
    NetworkScore: score(s) of the nodes, n x m matrix where 1st col is ID
    Candidates: ID(s) of candidates
    Targets: ID(s) of candidates
    maxDepth: maximal depth of the depth limited search
    file: output file name
    RetRanks - Boolean, if true results are returned as ranks instead of scores.
  
  Returns:
  
    A data frame with the Chains and their chain scores
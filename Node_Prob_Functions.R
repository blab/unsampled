# Packages
require(ape)
require(hash)
require(phytools)
require(ctv)
require(phylobase)
require(adephylo)
require(dplyr)
require(stringr)
require(ggplot2)

# Given n, the number of sampled tips, and popsize, the population size, this coalescent simulator steps back in time randomly generating coalescent events and returning a list of lengths, descendants, number of lineages from each node, a hash table mapping parent nodes to descendants, total generations (total tree depth), n, and the intervals between each branch.

cosim <-function(n, popsize)
  {
  #initialization
  #lengths of each branch
  blength <- rep(0,n)
  #number of descendants each node has
  descendants <- rep(1, n)
  #the node number of the parent, 0 if at the top of the tree
  parent <- rep(0, n)
  #set of current lineages that can be chosen to coalesce
  currentLins <- 1:n
  nextLin <- n + 1
  #total generations (counter)
  iter <- 0
  h <- hash()
  
  last.event <- 0
  interval.list <- rep(0,n-1)
  i=1
  
  while(length(currentLins)>2)
    {
    #condition of no coalescence
    if(runif(1) > choose(length(currentLins),2)/popsize) #currentLins*(currentLins-1)/(2*popsize)
      {
      #increases each branch length by 1 and increases the generation
      blength[currentLins] <- blength[currentLins] + 1
      iter <- iter + 1
      }
    #coalescence
    else
      {
      blength[currentLins] <- blength[currentLins] + 1
      iter = iter + 1
      interval.list[i] <- iter - last.event
      last.event <- iter
      i <- i + 1
      #randomly picking the nodes to coalesce
      rs <- sample(currentLins, 2, replace=FALSE)
      #debugging
      # cat(rs[1], " and ", rs[2], "coalesce after ", iter, " iterations\n")
      
      #set the new node to the next lineage number
      parent[rs] <- nextLin
      .set(h, keys=make.keys(nextLin), values=rs)
    
      #remove the two lineages that coalesced from the current list
      currentLins <-  setdiff(currentLins, rs)
      
      #record the new (parent) lineage back in at the end
      parent[nextLin] <- 0
      blength[nextLin] <- 0
      
      #update available lineages and update the next lineage number
      descendants[nextLin] <- descendants[rs[1]] + descendants[rs[2]]
      currentLins <- c(currentLins, nextLin)
      nextLin <- nextLin + 1
      
      #debugging
      #cat("branch lengths: ", blength, "\n")
      #cat("descendants: ", descendants, "\n")
      #cat("parents: ", parent, "\n")
      #cat("current lineages: ", currentLins, "\n")
      #cat("next lineage:", nextLin, "\n")
      #cat("generation: ", iter, "\n")
      }
  }
  
  if(length(currentLins)==2){
    
    while(runif(1) > choose(length(currentLins),2)/(2*popsize))
    {
      blength[currentLins] <- blength[currentLins] + 1
      iter <- iter + 1
    }
    blength[currentLins] <- blength[currentLins] + 1
    iter = iter + 1
    interval.list[i] <- iter - last.event
    
    #debugging
    # cat(currentLins[1], " and ", currentLins[2], "coalesce after ", iter, " iterations\n")
    #descendants[nextLin] <- descendants[currentLins[1]] + descendants[currentLins[2]]
    
     .set(h, keys=make.keys(0), values=currentLins[1:2])
    currentLins <- setdiff(currentLins, rs)
    
    #debugging
    # cat("branch lengths: ", blength, "\n")
    # cat("descendants: ", descendants, "\n")
    # cat("parents: ", parent, "\n")
    # cat("generation: ", iter, "\n")  
  }
  result <- list(blength=blength, descendants=descendants, parent=parent, iter=iter, h=h, n=n, interval.list = rev(interval.list))
  return(result)
}

# Traverse takes in a hash table and an optional starting root point (defaulted at the root of the tree) and does a recursive in-order tree traversal

traverse <- function(xhash, node=0){
  root <- node
  #cat("Root:", root, "\n")
  visited <<- visited[-match(root, visited)]
  leftChild <- min(values(xhash, keys=as.character(root)))
  #cat("left child:", leftChild, "\n")
  
  #to add tip labels for left children
  if(leftChild <=n){
      tip.label[cTip] <<- leftChild
      edge[cRow,1] <<- cN
      edge[cRow,2] <<- cTip
      
      #cN, nN stay the same
      #update new nodes
      newNodes$Freq[which(newNodes==cN)] <<- 1
      paste(newNodes, "\n")
      nN <<- min(which(newNodes==2))
      
      cRow <<- cRow + 1
      cTip <<- cTip+1
      #cat("ADD LEFT TIP \n")
  }
  else{
    edge[cRow,1] <<- cN
    edge[cRow,2] <<- nN
    newNodes$Freq[which(newNodes==cN)] <<- 1
    paste(newNodes, "\n")
    cN <<- nN
    nN <<- nN + 1
    #cat("current node: ", cN, " next node: ", nN, "\n")
    #cat("Root", root, " to left child", leftChild," \n")
    cRow <<- cRow + 1
  }
  
  if(length(which(visited==leftChild))>1){
    visited <<- visited[-match(leftChild, visited)]
    #cat("visited: pt 2 (LC)", visited, "\n")
    traverse(xhash, leftChild)
  }
  rightChild <- max(values(xhash, keys=as.character(root)))
   #cat("right child:", rightChild, "\n")
   
    #add tip labels for right children
   if(rightChild <=n){
    tip.label[cTip]<<-rightChild
    edge[cRow,1] <<- cN
    edge[cRow,2] <<- cTip
    newNodes$Freq[which(newNodes==cN)] <<- 0
    paste(newNodes, "\n")
    cN <<- max(which(newNodes==1))
    nN <<- min(which(newNodes==2))

    #cat("current node: ", cN, " current Tip: ", cTip, "\n")
    cRow <<- cRow +1
    cTip <<- cTip+1
  }
  else{
    edge[cRow,1] <<- cN
    edge[cRow,2] <<- nN
    newNodes$Freq[which(newNodes==cN)] <<- 0
    paste(newNodes, "\n")
    cN <<- nN
    nN <<- cN + 1
    
    #cat("current node: ", cN, " next node: ", nN, "\n")
    cRow <<- cRow +1
  }

  if(length(which(visited==rightChild))>1){
    visited <<- visited[-match(rightChild, visited)]
    #cat("visited: pt 3: (RC)", visited, "\n")
    traverse(xhash, rightChild)
  }
}

# Given the number of sampled tips and a population size, create.phylo() uses cosim() and traverse() to return a corresponding phylo object. Want n << N.

create.phylo <- function(n, popsize){
  options(warn=-1)
  s <- cosim(n, popsize)

  cTip=1
  cedge <- 1
  newNodes <- as.data.frame(table(rep((n+1):(n*2), each=2)))
  visited <- c(rep(0,1), rep((n+1):(2*n-2), each=2))
  cN <- min(which(newNodes==2))
  nN <- cN + 1
  cRow <- 1
  
  tip.label <- rep(0,n)
  edge.length <- rep(0,2*n-2)
  edge <- matrix(nrow=2*n-2, ncol=2)
  nNode <- n-1
  
  traverse <- function(xhash, node=0){
    root <- node
    #cat("Root:", root, "\n")
    visited <<- visited[-match(root, visited)]
    leftChild <- min(values(xhash, keys=as.character(root)))
    #cat("left child:", leftChild, "\n")
    
    #add current left child edge
    edge.length[cedge] <<- s$blength[leftChild]
    cedge <<- cedge + 1
    #to add tip labels for left children
    if(leftChild <=n){
        tip.label[cTip] <<- leftChild
        edge[cRow,1] <<- cN
        edge[cRow,2] <<- cTip
        
        #cN, nN stay the same
        #update new nodes
        newNodes$Freq[which(newNodes==cN)] <<- 1
        paste(newNodes, "\n")
        nN <<- min(which(newNodes==2))
        
        cRow <<- cRow + 1
        cTip <<- cTip+1
        #cat("ADD LEFT TIP \n")
    }
    else{
      edge[cRow,1] <<- cN
      edge[cRow,2] <<- nN
      newNodes$Freq[which(newNodes==cN)] <<- 1
      paste(newNodes, "\n")
      cN <<- nN
      nN <<- nN + 1
      #cat("current node: ", cN, " next node: ", nN, "\n")
      #cat("Root", root, " to left child", leftChild," \n")
      cRow <<- cRow + 1
    }
    
    if(length(which(visited==leftChild))>1){
      visited <<- visited[-match(leftChild, visited)]
      #cat("visited: pt 2 (LC)", visited, "\n")
      traverse(xhash, leftChild)
    }
    rightChild <- max(values(xhash, keys=as.character(root)))
     #cat("right child:", rightChild, "\n")
     
    #add current right child edge
    edge.length[cedge] <<- s$blength[rightChild]
    cedge <<- cedge + 1
     
      #add tip labels for right children
     if(rightChild <=n){
      tip.label[cTip]<<-rightChild
      edge[cRow,1] <<- cN
      edge[cRow,2] <<- cTip
      newNodes$Freq[which(newNodes==cN)] <<- 0
      paste(newNodes, "\n")
      cN <<- max(which(newNodes==1))
      nN <<- min(which(newNodes==2))
  
      #cat("current node: ", cN, " current Tip: ", cTip, "\n")
      cRow <<- cRow +1
      cTip <<- cTip+1
    }
    else{
      edge[cRow,1] <<- cN
      edge[cRow,2] <<- nN
      newNodes$Freq[which(newNodes==cN)] <<- 0
      paste(newNodes, "\n")
      cN <<- nN
      nN <<- cN + 1
      
      #cat("current node: ", cN, " next node: ", nN, "\n")
      cRow <<- cRow +1
    }
  
    if(length(which(visited==rightChild))>1){
      visited <<- visited[-match(rightChild, visited)]
      #cat("visited: pt 3: (RC)", visited, "\n")
      traverse(xhash, rightChild)
    }
  }
  traverse(s$h)
  
  result <- list(edge=edge, edge.length=edge.length, tip.label = tip.label, Nnode = nNode)
  class(result) <- "phylo"
  return(result)
}

# Takes in a tree and finds the length of each interval. Returns a matrix with the first column the length of the current interval, column 2 as the origin internal node of that interval and column 3 the number of active lineages for that interval.

interval.sort <- function(tr){
  #node to tip distance
  br <- branching.times(tr)
  actual.tbr <- matrix(NA,nrow = length(br), ncol=3)
  
  
  for(k in 1:(length(br))){
    oldest.event <- as.numeric(max(br))
    actual.tbr[k,2] <- as.numeric(names(br)[which(br==oldest.event)])
    actual.tbr[k,3] <- k+1
    br <- br[-which(br==max(br))]
    if(length(br)>=1){
      actual.tbr[k,1] <- (oldest.event - max(br))
    }
    else{
      actual.tbr[k,1] <- oldest.event
    }
  }
  return(actual.tbr)
}

# Returns the probability of coalescence within an interval for each internal node

interval.pcoal.cont <- function(tr, popsize){
  ints <- interval.sort(tr)
  pcoal <- rep(0, length(ints[,1]))
  
  for(i in 1:length(ints[,1])) { 
    pcoal[i] <- 1-exp(-(ints[i,3]/popsize)*ints[i,1])
  }
  return(pcoal)
}

one.lin.p <- function(t, popsize){

  full <- interval.sort(t)
  one.lin.prob <- interval.pcoal.cont(t, popsize)
  nocoal <- 1-one.lin.prob
  one.lin.prob[length(one.lin.prob)] <- one.lin.prob[length(one.lin.prob)]/full[length(full[,3]), 3]
  
  for(i in 1: length(one.lin.prob)-1){
    for(j in (i+1):length(one.lin.prob)){
      one.lin.prob[i] <- one.lin.prob[i]*nocoal[j]
    }
    one.lin.prob[i] <- one.lin.prob[i]/full[i,3]
  }
  return(one.lin.prob)
}

# For an internal node, calculates the probability of coalescing under this subtree rooted at the internal node

node.prob <- function(t, subtree, popsize){
  full<-interval.sort(t)
  sub.inter <- interval.sort(subtree)
  
  #probability for a single lineage segment in a given interval
  #probability of coalescence for each interval over the number of active lineages in that interval multiplied by the probability of not coalescing in any of the intervals before the current one
  lin.p <- one.lin.p(t,popsize)
  
  #initialize the node probability
  #start at the row in the full tree corresponding to the root of the subtree
  prob <- 0
  current <- which(full[,2]==sub.inter[1,2])
  k=1

  while(k <= subtree$Nnode){
    full.tree.sum <- full[current,1]
    c.s.inter <- sub.inter[k,1]

    while(c.s.inter >= full.tree.sum && current <= length(full[,1])){
      prob <- prob + sub.inter[k,3]*lin.p[current]
      current <- current + 1
      if(current <= length(full[,1])){
        full.tree.sum <- full.tree.sum + full[current,1]
      }
    }
    k <- k+1
  }
  return(prob)
}


# For Non-Ultrametric Trees:
# These functions travel tip to root. Intervals with nodes labels less than the number of tips are intervals that are started by tips on the right side of the  interval.
# You can't travel from internal node to internal node, you must travel by intervals (toward the root)

# Takes in a tree (allows for non-ultrametric) and returns a matrix of distance from each node or tip to the root, ordered by decreasing length

dist.to.root <- function(t){
  m <- matrix(NA, nrow=length(t$tip.label)+t$Nnode, ncol=2)
  colnames(m) <- c("dist","node")
  
  #Fill in the distance matrix for the tips using distRoot()
  a <- distRoot(t)
  for(i in 1:length(a)){
    m[i,2] <- i
    m[i,1] <- a[[i]]
  }
  
  #Fill in the root length as 0
  root <- length(t$tip.label) + 1
  m[root,2] <- root
  m[root,1] <- 0
  
  #Fill in the internal nodes in the matrix
  #Check first to see if you have a tree with internal nodes
  if(length(m[,1])>3){
    #starting from the internal node after the root and going to the last internal node
    for(node in (length(t$tip.label)+2):(length(t$tip.label)+t$Nnode)){
      #find the shortest path between nodes to get from the root to the node of interest. If there are internal nodes between the node of interest, add their lengths
      #otherwise just record the length of the node of interest
      p <- shortestPath(t, root, node)
      if(length(p)>0){
        length <- t$edge.length[which(t$edge[,2]==node)]
        for(j in 1:length(p)){
          length <- length + t$edge.length[which(t$edge[,2]==p[[j]])]
        }
      }
      else{
        length <- t$edge.length[which(t$edge[,2]==node)]
      }
      
      m[node, 2] <- node
      m[node, 1] <-length
    }
  }
  #return an ordered distance from root, node paired matrix
  return(m[order(m[, "dist"], decreasing=TRUE),])
}

# Like the ultrametric interval sort, takes in a tree and returns a matrix with the interval length in the first column, the node associated with that length in the second, and the number of active lineages in the third

nu.interval.sort <- function(t){
  x <- dist.to.root(t)
  intervals <- matrix(NA, nrow = length(x[,1])-1, ncol=3)
  colnames(intervals) <- c("len","node", "lineages")
  lins <- 0
  
  #For each row in the distance to the root, takes the difference between consecutive rows
  for(i in 1:(length(x[,1])-1)){
    intervals[i,1] <- x[i,1] - x[(i+1),1]
    intervals[i,2] <- x[i,2]
    
    #if its a tip, add a new lineage to the current number of lineages
    #otherwise, it's a node and you should subtract a lineage from the list
    if(x[i,2] <= length(t$tip.label)){
      lins <- lins + 1
    }
    else{
      lins <- lins - 1
    }
    intervals[i,3] <- lins
  }
  return(intervals)
}

# Returns the probability of a coalescent event in a specific interval using the continuous time (exponentially distributed) probability of the coalescent math

nu.pcoal.cont <- function(t, popsize){
  ints <- nu.interval.sort(t)
  pcoal <- rep(0, length(ints[,1]))
  #for each row in the sorted interval, find p(Coal) = 1-e^(-k/N*t)
  for(i in 1:length(ints[,1])) { 
    pcoal[i] <- 1-exp(-(ints[i,3]/popsize)*ints[i,1])
  }
  return(pcoal)
}

# Returns the probability of coalescing onto one specific lineage in an interval. Takes into account all of the previous intervals seen, i.e. the probability of coalescing in interval i given that the event did not occur in the intervals i-1 to 1.

nu.one.lin.p <- function(t, popsize){
  full <- nu.interval.sort(t)
  one.lin.prob <- nu.pcoal.cont(t, popsize)
  nocoal <- 1-one.lin.prob
  one.lin.prob[1] <- one.lin.prob[1]/full[1,3]  
  
  for(i in 2: length(one.lin.prob)){
    for(j in (i-1):1){
      one.lin.prob[i] <- one.lin.prob[i]*nocoal[j]
    }
    one.lin.prob[i] <- one.lin.prob[i]/full[i,3]
  }
  return(one.lin.prob)
}

# For an internal node, calculates the probability of coalescing under this subtree rooted at the internal node. Takes in a non-ultrametric tree. 
# Note: every previous function for non-ultrametric trees works from tip to root, but this function takes that information and works backwards (essentially mimicing the root to tip traveling in the ultrametric tree interval traversal)

nu.node.prob <- function(t, subtree, popsize){
  f.sort <- nu.interval.sort(t)
  s.sort <- nu.interval.sort(subtree)
  lin.p <- nu.one.lin.p(t,popsize)
  prob <- 0
  #end.interval <- which(distRoot(x)==max(distRoot(x)))
  #initialize start points in both trees
  if(as.numeric(subtree$name)==(length(t$tip.label)+1)){
    #deals with starting on the full tree (the interval sort does not have a row for the root because the distance to the root is 0)
    f.index <- length(f.sort[,2])
  }
  else{
    f.index <- which(f.sort[,2]==as.numeric(subtree$name))-1
  }
  s.index <- length(s.sort[,2])
  # cat("Initialization- Full tree index: ", f.index, "\n")
  # cat("Initialization- Subtree index: ", s.index, "\n")

  #Must add the probability of coalescing in each interval defined by the full tree. However, the number of lineages in the full tree is not necessarily the same as in the subtree. Thus, we sum branches of the full tree until we reach the length of an interval in the subtree. Then we update the number of current lineages (by the number in the subtree)
  while(s.index > 0){
    f.len <- f.sort[f.index, 1]
    s.inter <- s.sort[s.index, 1]
    # cat("Outer while S.INDEX: ", s.index, "\n")
    # cat("Full Length: ", f.len, "\n")
    # cat("subtree interval length: ", s.inter, "\n")
    
    while((f.index>0)&&(f.len<=s.inter || abs(f.len-s.inter) <= 1e-15)){
      # cat("Inner while F.INDEX: ", f.index, "\n")
      # cat("Inner while F.LEN: ", f.len, "\n")
      # cat("Inner while S.INTER: ", s.inter, "\n")
      #update probability to add the number of lineages in the current interval times the probability of coalescence in that interval for one lineage
      prob <- prob + s.sort[s.index,3]*lin.p[f.index]
      # debugging
      # cat("Probability: ", prob, "\n")
      # cat("Added prob: ", s.sort[s.index,3]*lin.p[f.index], "\n")
      

      f.index <- f.index - 1
      f.len <- f.len + f.sort[f.index, 1]
      # cat("Full Tree Length:", f.len, "\n")
      # cat("Full index:", f.index, "\n")
      
    }
    s.index <- s.index - 1
    # cat("Subtree index:", s.index, "\n")
  }
  return(prob)
}

# Data distribution across nodes simulation

# Creates a hash table to find the lineage names in each interval. The interval is marked by the internal node at the left of the interval (the key), the lineages in it are denoted by the first node or tip to the right of the interval start (keys)

lins.in.interval <- function(t, ints=interval.sort(t)){
  h <- hash()
  .set(h, keys=ints[1,2], values=t$edge[which(t$edge[,1]==ints[1,2]), 2])
  
  for(inter in 2:length(ints[,1])){
    vals <- values(h, keys=ints[inter-1,2])
    vals <- vals[-which(vals==ints[inter,2])]
    vals <- c(vals, t$edge[which(t$edge[,1]==ints[inter,2]), 2])
    .set(h, keys=ints[inter, 2], values=vals)
  }
  return(h)
}

# Creates a matrix where each column is an internal node and each row is a simulated run. Randomly adds a data point to the tree based on the coalescent and records where the nodes is falls under at each run. Returns an average of where the added data point falls.  This corresponds to the algorithmic model of where each node falls in node.prob()

interval.sim <- function(t, popsize, numsim=1000){
  sims <- matrix(0, nrow=numsim, ncol=t$Nnode)
  colnames(sims) <- c((length(t$tip.label)+1):(2*length(t$tip.label)-1))
  ints <- interval.sort(t)
  p.coal <- interval.pcoal.cont(t, popsize)
  active.lineages <- lins.in.interval(t, ints)
  
  
  for(i.sim in 1:numsim){
    interval <- length(p.coal)
    coal <- FALSE
    while(!coal){
      if(interval >= 1){
        if(runif(1) < p.coal[interval]){
          #coalescent event in the interval DOES occur
          lin <- sample(values(active.lineages, keys=as.character(ints[interval, 2])), 1)
          currentNode <- t$edge[which(t$edge[,2]==lin), 1]
          #sims[i.sim, t$Nnode+1] <- currentNode
          root <- t$edge[1,1]
          
          while(currentNode != root){
            sims[i.sim,currentNode - length(t$tip.label)] <- 1
            currentNode <- t$edge[which(t$edge[,2]==currentNode),1]
          }
          sims[i.sim,root - length(t$tip.label)] <- 1
          coal <-TRUE
        }
      }
      else{
        coal <- TRUE
      }
    interval <- interval - 1    
    }
  }  
  return(colMeans(sims))
}

# drop.t takes in a tree and the number of tips you want to randomly drop and plots the corresponding dropped tips as well as returning the resulting tree and the dropped tips in a list. drop.t also takes in an optional parameter p with the default of true which detotes whether to plot the changes in the tree.

drop.t <- function(tr, ndtips, p=TRUE){
  if(ndtips < length(tr$tip.label)){
    dtips <- sample(tr$tip.label, ndtips)
  }
  if(p){
    lty<-rep(1,nrow(tr$edge))
    clr <- rep("black", nrow(tr$edge))
    tip.clr <- rep("black", length(tr$tip.label))
    tip.clr[dtips] <- "seagreen"
    for(i in 1:length(dtips)){
      lty[which(tr$edge[,2]==dtips[i])]<-2
      clr[which(tr$edge[,2]==dtips[i])]<- "springgreen3"
    }
    plot(tr, edge.lty=lty, edge.color=clr, tip.color=tip.clr, edge.width=1.4)
  }
  dt <- drop.tip(tr, dtips)
  #par(new=TRUE)
  #plot(dt, edge.color="cornflowerblue", tip.color="dodgerblue4")
  result <- list(tr=dt, dropped=tr$tip.label[dtips])
  return(result)
}

# tree.overlay takes in a tree and the number of tips to drop and plots the new phylogenetic tree over the original one in different colors.

tree.overlay <- function(tr, ndtips){
  if(ndtips < length(tr$tip.label)){
    dtips <- sample(tr$tip.label, ndtips)
  }
  lty<-rep(1,nrow(tr$edge))
  clr <- rep("black", nrow(tr$edge))
  tip.clr <- rep("black", length(tr$tip.label))
  tip.clr[dtips] <- "seagreen"
  for(i in 1:length(dtips)){
    lty[which(tr$edge[,2]==dtips[i])]<-2
    clr[which(tr$edge[,2]==dtips[i])]<- "springgreen3"
  }
  plot(tr, edge.lty=lty, edge.color=clr, tip.color=tip.clr, edge.width=1.4)
  dt <- drop.tip(tr, dtips)
  par(new=TRUE)
  plot(dt, edge.color="cornflowerblue", tip.color="dodgerblue4")
  # result <- list(tr=dt, dropped=tr$tip.label[dtips])
  # return(result)
}

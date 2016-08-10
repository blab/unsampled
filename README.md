# Unsampled

## Introduction

In viral outbreaks, some infections are likely to be missed. Finding as many cases as possible is an important foundation for developing better treatment and containing the spread of pathogens. This project is ultimately designed to highlight clades containing possible undersampled pockets. 

To achieve this, we first calculate the proportion of data expected 
under each internal node in a given tree. We extend this method to target nodes and branches within the tree expected to be over or under sampled based upon z-scores of these nodal probabilities.  For a given tree structure, conditional probabilities describing the chance of seeing a new, unknown sample coalescing under a given node in the tree are calculated based on effective population size, interval lengths (captured by sampling and coalescent events), and number of lineages within intervals. 

## Files:

### Node_Prob_Functions.Rmd

This file is a compilation of methods which include a coalescence simulator for trees sampled at the same time, conversion of the simulator into a phylo object accepted by apein R, simulation of where an additional unknown sample would fall for an ultrametric tree,and the proportion of samples expected under any internal node in a given phylogenetic tree with a known effective population size.  

Methods Included:
#### For Ultrametric Trees 
cosim(n, popsize)  
  *Description: Uses the coalescent framework to simulate an ultrametric tree structure. Returns an object with a list of branch lengths, descendants, number of lineages from each node, hash table mapping parent nodes to descendants, total generations, and intervals between branches.  
  *Arguments: n-number of samples, popsize-effective population size
traverse(xhash, node=0)  
  *Description: Recursively performs in-order tree traversal  
  *Arguments: xhash-the hash table with key value pairs of nodes and their descendants, node-node to start traversal at, default at 0    
create.phylo(n, popsize)  
  *Description: Creates a phylo object accepted by ape using cosim() and traverse() functions for ultrametric trees.  
  *Arguments: n-number of samples, popsize-effective population size  
interval.sort(tr)  
  *Description: returns a matrix with three columns corresponding to length of current interval, origin node for that interval, and number of lineages in that interval respectively.  
  *Arguments: tr-a phylo object tree  
interval.pcoal.sort(tr, popsize)  
  *Description: Returns the probability of coalescence within an interval for each internal node using an exponential distribution to model probabilities of coalescence  
  *Arguments: tr-a phylo object tree, popsize-effective population size    
one.lin.p(t, popsize)  
  *Description: Returns a list of conditional probabilities for each interval of the probability of landing on a specific lineage in an interval given that it didn't coalesce to the tree in a previous interval  
  *Arguments: t-a phylo object tree, popsize- effective population size   
node.prob(t, subtree, popsize)  
  *Description: for an internal node, returns the proportion of samples expected to coalesce under this subtree rooted at the internal node  
  *Arguments: t-a phylo object tree, subtree- a phylo object subtree rooted at the node of interest, popsize- effective population size 
drop.t(tr, ndtips, p=TRUE)
  *Description: Returns a phylo object tree with ndtips randomly dropped tips. If p=TRUE, plots the difference between trees.
  *Arguments: t-a phylo object tree, ndtips- number of tips to drop, p-whether to plot the change between trees  
tree.overlay(tr, ndtips)
  *Description: Plots the overerlay between a tree and the resulting tree with ndtips fewer tips
  *Arguments: tr- a phylo object tree, ndtips- number of tips to drop   
##### For Simulation (of Ultrametric Trees)  
lins.in.interval(t, ints=interval.sort(t))  
  *Description: Returns a hash table with the lineage names in each interval. The interval serves as the key and is denoted by the node at the left boundary of the interval. The lineages in it are values, denoted to the first top or node to the right of the interval boundary.  
  *Arguments: t-a phylo object tree, ints-a matrix of sorted intervals returned by the interval.sort(t) method    
interval.sim(t, popsize, numsim=1000)  
  *Description: Returns a list of average proportion of times a new, unknown sample would land under a particular node. Does this through storing a matrix of numsim rows and number of nodes columns and recording where a randomly added sample falls.
  *Arguments: t-a phylo object tree, popsize- effective population size, numsim- number of simulations desired  

#### For Time Sampled Trees  
dist.to.root(t)
  *Description: Returns a matrix of node or tip label and the corresponding distance to the root of the tree, ordered by decreasing distance.
  *Arguments: t-a phylo object tree  
nu.interval.sort(t)
  *Description: returns a matrix with three columns corresponding to length of current interval, origin node for that interval, and number of lineages in that interval respectively.  
  *Arguments: tr-a phylo object tree
nu.pcoal.cont(t, popsize) 
  *Description: Returns the probability of coalescence within an interval for each internal node using an exponential distribution to model probabilities of coalescence  
  *Arguments: tr-a phylo object tree, popsize-effective population size 
nu.one.lin.p(t, popsize)
  *Description: Returns a list of conditional probabilities for each interval of the probability of landing on a specific lineage in an interval given that it didn't coalesce to the tree in a previous interval  
  *Arguments: t-a phylo object tree, popsize- effective population size   
nu.node.prob(t, subtree, popsize)  
  *Description: for an internal node, returns the proportion of samples expected to coalesce under this subtree rooted at the internal node  
  *Arguments: t-a phylo object tree, subtree- a phylo object subtree rooted at the node of interest, popsize- effective population size


### Ebola_Example_Visualizations.Rmd

Uses data from the Ebola Sequence Consortium to demonstrate the methodology introduced
in Node_Prob_Functions, quantifying proportion of samples expected under internal nodes.
Includes visualizations of the labeled Ebola tree as well as calculation and visualization of 
Z-scores measuring the difference between proportion of tips seen versus expected under a node.
Z-scores serve as a metric for whether a node is predicted to be undersampled.  
Requires methods from Node_Prob_Functions and the Conakry Ebola subtree Nexus file to run.

#### Visualizations:  
Ebola tree with labeled node proportions   
  *Colored by proportions    
Histogram of node Z-scores showing outliers 
  *Corresponding histogram without tip 38   
Corresponding Ebola tree with colored outliers  
  *Corresponding tree without tip 38
Scatterplot comparing predicted to actual proportion of tips under each node, colored by clade, sized by abs(z-score)
  *Corresponding tree colored by clade
Toy tree plot and node calculations
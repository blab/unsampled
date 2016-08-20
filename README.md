# Using the coalescent to detect missed infections in phylogenetic trees

#### Stephanie Stacy<sup>1,2</sup>, Allison Black<sup>1,3</sup>, Gytis Dudas<sup>1</sup>, Trevor Bedford<sup>1</sup>

<sup>1</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA, <sup>2</sup>Williams College, Williamstown, MA, USA, <sup>3</sup>Department of Epidemiology, University of Washington, Seattle, WA, USA

## Introduction

In viral outbreaks, some infections are likely to be missed. Finding as many cases as possible is an important foundation for developing better treatment and containing the spread of pathogens. This project is ultimately designed to highlight clades containing possible undersampled pockets.

To achieve this, we first calculate the proportion of data expected 
under each internal node in a given tree. We extend this method to target nodes and branches within the tree expected to be over or under sampled based upon z-scores of these nodal probabilities.  For a given tree structure, conditional probabilities describing the chance of seeing a new, unknown sample coalescing under a given node in the tree are calculated based on effective population size, interval lengths (captured by sampling and coalescent events), and number of lineages within intervals.

## Build

Install necessary R packages with:

    Rscript install.R

`.Rmd` files can be knit in `.md` files with:

    Rscript -e "library(knitr); knit('Toy_Tree_Example.Rmd')"
    Rscript -e "library(knitr); knit('Ebola_Example_Visualizations.Rmd')"

## Coalescent algorithms

Core functions are in the file [`Node_Prob_Functions.R`](Node_Prob_Functions.R). The file [`Node_Prob_Functions.md`](Node_Prob_Functions.md) describes these functions.

More mathematical detail can be found in the [SURP poster](poster/surp_poster.png)

## Toy example

Code is in the file [`Toy_Tree_Example.Rmd`](Toy_Tree_Example.Rmd) and knitted results are in the file [`Toy_Tree_Example.Rmd`](Toy_Tree_Example.Rmd).

## Ebola analysis

Code is in the file [`Ebola_Example_Visualizations.Rmd`](Ebola_Example_Visualizations.Rmd) and knitted results are in the file [`Ebola_Example_Visualizations.Rmd`](Ebola_Example_Visualizations.Rmd).

This uses data from the Ebola Sequence Consortium to demonstrate the methodology introduced in `Node_Prob_Functions`, quantifying proportion of samples expected under internal nodes. Includes visualizations of the labeled Ebola tree as well as calculation and visualization of  Z-scores measuring the difference between proportion of tips seen versus expected under a node. Z-scores serve as a metric for whether a node is predicted to be undersampled. Requires methods from `Node_Prob_Functions` and the Conakry Ebola subtree file [`EBOV_Conakry_subtree.tree`](EBOV_Conakry_subtree.tree) to run.

The Ebola subtree of the Conakry lineage was extracted from the BEAST analysis conducted by [Dudas et al. 2016](https://github.com/ebov/space-time).

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

    Rscript -e "library(knitr); setwd('toy-tree'); knit('Toy_Tree.Rmd', 'README.md')"
    Rscript -e "library(knitr); setwd('ebola-tree'); knit('Ebola_Tree.Rmd', 'README.md')"

## [Coalescent functions](functions/)

These are the core functions of the method. More mathematical detail can be found in the [SURP poster](poster/surp_poster.png).

## [Toy example](toy-tree/)

This is a simple tree example to demonstrate the algorithm.

## [Ebola analysis](ebola-tree/)

This is an example using a lineage of Ebola circulating in Conakry during the West African outbreak.

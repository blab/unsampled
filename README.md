# Using the coalescent to detect missed infections in phylogenetic trees

#### Stephanie Stacy<sup>1,2</sup>, Allison Black<sup>1,3</sup>, Gytis Dudas<sup>1</sup>, Trevor Bedford<sup>1</sup>

<sup>1</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA, <sup>2</sup>Williams College, Williamstown, MA, USA, <sup>3</sup>Department of Epidemiology, University of Washington, Seattle, WA, USA

## Introduction

In viral outbreaks, some infections are likely to be missed. Finding as many cases as possible is an important foundation for developing better treatment and containing the spread of pathogens. This project is ultimately designed to highlight clades containing possible undersampled pockets.

To achieve this, we first calculate the proportion of data expected
under each internal node in a given tree. We extend this method to target nodes and branches within the tree expected to be over or under sampled based upon z-scores of these nodal probabilities.  For a given tree structure, conditional probabilities describing the chance of seeing a new, unknown sample coalescing under a given node in the tree are calculated based on effective population size, interval lengths (captured by sampling and coalescent events), and number of lineages within intervals.

## Organization

* Analysis code in [analysis](analysis/)
* Manuscript in [manuscript](manuscript/)

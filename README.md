
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package graphppm

<!-- badges: start -->
<!-- badges: end -->

Graph product partition models (GraphPPM) is a probabilistic clustering
method for graph structured data with node-attributes. Useful functions
include counting and generating uniform spanning trees using Wilson’s
algorithm.

## Installation

You can install the development version of graphppm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("changwoo-lee/graphppm")
```

## Example 1: counting and sampling uniform spanning trees

``` r
library(graphppm)
library(igraph)
# base graph
graph0 = igraph::make_lattice(c(5,5))
mylayout = as.matrix(expand.grid(1:5,1:5))
plot(graph0, layout = mylayout)

# count the total number of spanning trees using matrix-tree theorem
nsptrees_igraph(graph0, log = F)
# cross-check with, e.g. https://oeis.org/A007341

# sample from the uniform spanning tree 
sptree = runifsptree(graph0, return.igraph = T) 
# Our code is much faster than igraph::sample_spanning_tree(); see below benchmark

plot(graph0, layout = mylayout)
plot(sptree, layout = mylayout, vertex.color = NA, vertex.label = NA, edge.color = "red", edge.width = 2,  add = T)

# sampling from the uniform spanning tree, given the partition 
z = c(rep(1,15),rep(2,10))
sptree_given_z = runifsptree_given_z(graph0, z, return.igraph = T) # output: (n-1) x 2 matrix of edge list

plot(graph0, vertex.color = z, layout = mylayout)
plot(sptree_given_z, layout = mylayout, vertex.color = NA, vertex.label = NA, edge.color = "red", edge.width = 2, add = T)
```

## Example 2:

This is a basic example with Zachary’s karate club dataset.

``` r

# prior analysis with spanning tree based cohesion function 
library(graphppm)
graph0 = make_graph("Zachary")
# cohesion function examples
source("https://raw.github.com/changwoo-lee/graphppm/main/R/logcohesions.R")

cohesion_param = list(alpha = 1, psi = 1)
# using importance sampler
fit = graphppm_prior_impsamp(graph0, 
                             logcohesion = logcohesion_vetau, 
                             cohesion_param = cohesion_param, 
                             nsave = 10000)
# using split-merge MCMC
fit2 = graphppm_prior(graph0, 
                             logcohesion = logcohesion_vetau, 
                             cohesion_param = cohesion_param, 
                             nsave = 10000, nburn = 1000, nthin = 1)
par(mfrow = c(1,2))
plot(fit$kpmf[1:5])
plot(fit2$kpmf[1:5])

library(fields)
image.plot(fit$simmatrix)
image.plot(fit2$simmatrix)
```

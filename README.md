
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package graphppm

<!-- badges: start -->
<!-- badges: end -->

Graph product partition models (GraphPPM) is a probabilistic clustering
method for graph structured data with node-attributes.

## Installation

You can install the development version of graphppm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("changwoo-lee/graphppm")
```

## Example

This is a basic example with Zachary’s karate club dataset.

``` r
library(graphppm)
library(igraph)
graph0 = make_graph("Zachary")

# prior analysis with spanning tree based cohesion function 
cohesion_param = list(alpha = 1, psi = 1)

fit = graphppm_prior_impsamp(graph0, 
                             logcohesion = logcohesion_vetau, 
                             cohesion_param = cohesion_param, 
                             nsave = 5000)
```

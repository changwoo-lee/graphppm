
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
# see R/logcohesions.R for examples
logcohesion_vetau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  
  nj = ncol(A)
  if(minus_logtau && (psi == 1)){  # no need to calculate tau(G_j)
    out = log(alpha) -0.5*log(nj) + log(sum(A)/2)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = -0.5*log(nj) + log(sum(A)/2) + logtau
    if(minus_logtau){
      out = log(alpha) + psi*logc0 - logtau
    }else{
      out = log(alpha) + psi*logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}


cohesion_param = list(alpha = 1, psi = 1)

fit = graphppm_prior_impsamp(graph0, 
                             logcohesion = logcohesion_vetau, 
                             cohesion_param = cohesion_param, 
                             nsave = 5000)
```


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
graph0 = make_graph("Zachary")

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

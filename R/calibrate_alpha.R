#' Calibration of alpha based on k, the number of clusters
#'
#' @param fit fitted object using importance sampler. Must have "save_logweight", "save_k", and "cohesion_param$alpha"
#' @param alpha vector of alpha values. 
#'
#' @return 
#' @export
#'
#' @examples
calibrate_alpha_k <- function(fit, alpha){
  if(any(alpha <=0)) stop("alpha must be positive")
  logweight = fit$save_logweight
  save_k = fit$save_k
  nsave = length(fit$save_k)
  nalpha = length(alpha)
  
  alpha_base = fit$cohesion_param$alpha
  if(alpha_base != 1){
    print("alpha used in the importance sampler is not 1, adjusted accordingly")
    logweight = logweight - save_k*log(alpha_base)
  }
  mat1 = matrix(logweight, nrow = nsave, ncol = nalpha) + # logweight
    t(log(alpha)*t(matrix(save_k, nrow = nsave, ncol = nalpha))) # save_k*log(alpha)
  mat1_normalized = t(t(mat1) - matrixStats::colLogSumExps(mat1))
  mat2 = matrix(log(save_k), nrow = nsave, ncol = nalpha) # log(save_k) (functional)
  k_alpha = exp(matrixStats::colLogSumExps(mat1 + mat2) - matrixStats::colLogSumExps(mat1))
  
  # (9.9) of owen importance sampling paper
  k_alpha_var = exp(matrixStats::colLogSumExps( 2*mat1_normalized + 2*log(abs(exp(mat2) - matrix(k_alpha, nrow = nsave, ncol = nalpha, byrow = T))))) 
  out = list()
  out$ncluster_mean = k_alpha
  out$ncluster_mean_sd = sqrt(k_alpha_var)
  out$nsave = nsave
  out
}
# 
# library(graphppm)
# library(igraphdata)
# data("karate")
# plot(karate)
# load("data/HVR6_connected.rda")
# g = HVR6_connected
# # cohesion function examples
# source("https://raw.github.com/changwoo-lee/graphppm/main/R/logcohesions.R")
# par(mfrow = c(1,2))
# cohesion_param = list(alpha = 1, psi = 1)
# fit_tau = graphppm_prior(g,
#                                    logcohesion = logcohesion_tau,
#                                    cohesion_param = cohesion_param,
#                                    nsave = 10000)
# 
# out = calibrate_alpha_ncluster(fit_tau, alpha = seq(0.01, 1, length = 20))
# 
# 
# 
# 
# fit_vetau = graphppm_prior_impsamp(g,
#                                  logcohesion = logcohesion_vetau,
#                                  cohesion_param = cohesion_param,
#                                  nsave = 10000)
# 
# fit_vetau$save_logweight
# out2 = calibrate_alpha_ncluster(fit_vetau, alpha = seq(0.01, 1, length = 20))
# 
# par(mfrow = c(1,1))
# plotrix::plotCI(seq(0.01, 1, length = 20),
#                 out$ncluster_mean,
#                 ui = out$ncluster_mean + 2*out$ncluster_mean_sd,
#                 li = pmax(1, out$ncluster_mean - 2*out$ncluster_mean_sd), 
#                 xlab = "alpha", ylab = "ncluster", main = "tau(G_j)")
# 
# plotrix::plotCI(seq(0.01, 1, length = 20),
#                 out2$ncluster_mean,
#                 ui = out2$ncluster_mean + 2*out2$ncluster_mean_sd,
#                 li = pmax(1, out2$ncluster_mean - 2*out2$ncluster_mean_sd), 
#                 xlab = "alpha", ylab = "ncluster", main = "|e_j| tau(G_j)", add = T, col = 2)
# 
# 

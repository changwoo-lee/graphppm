#' GraphPPM prior analysis with importance sampler
#'
#' @param graph0 igraph object, base graph
#' @param logcohesion logcohesion function
#' @param cohesion_param cohesion function parameters
#' @param nsave number of samples
#' @param kvals range of k
#' @param kprobs unnormalized probability of q(k)
#' @param nCores default 1, set 2 or more lead to parallel run
#' @param ... 
#'
#' @return
#' @import igraph matrixStats
#' @export
#'
#' @examples
graphppm_prior_impsamp <- function(graph0, logcohesion, cohesion_param = NULL, nsave = 2000,
                           kvals = NULL, kprobs = NULL, nCores = 1, ...){
  
  if(is.null(cohesion_param)){
    cohesion_param = list(alpha = 1, psi = 1)
  }
  
  
  if(nCores >= 2) print(paste("Run parallel, nCores set as ",nCores))
  if(nsave %% nCores !=0){
    warning("nsave is not a multiple of nCores. reduce nsave a bit..")
    nsave = nparticle - (nparticle %% nCores)
    print(paste("set nsave = ",nsave))
  }
  if(nCores-1 > parallel::detectCores()){
    stop("too large nCores, should be less than or equal to detectCores() - 1")
  }
  
  
  n = vcount(graph0)
  m = ecount(graph0)
  
  if(is.null(kvals)) kvals = 1:(n-1)
  if(is.null(kprobs)) kprobs = 0.7^kvals
  
  E0 = as_edgelist(graph0, names =F)
  A0 = as_adj(graph0, sparse = F, names = F) # sparse = F, lookup speed is faster in non-sparse form
  g0 = graph_from_edgelist(E0, directed = F)
  deg0 = degree(g0) # external 
  
  logncluster = -Inf
  lognclustersq = -Inf
  
  save_k = integer(nsave)
  save_logweight = numeric(nsave)
  
  logsimmatrix = matrix(-Inf, n, n)
  logsimmatrixsq = matrix(-Inf, n, n)
  
  
  # select good parameters for wilson algorithm..
  Wilson_root = which.max(deg0)
  Wilson_order = as.integer(bfs(g0,  Wilson_root)$order)
  Wilson_maxiter = as.integer(n^3)
  #browser()
  
  ########### non-parallel setting ###############
  if(nCores == 1){
    for(i in 1:nsave){
      
      # 1. Draw UST cut partition 
      # 1-1. unif spanning tree
      edge_list = graphppm:::runcppWilson(n, E0, Wilson_root, Wilson_order, Wilson_maxiter)
      # 1-2. k
      kprobs = kprobs/sum(kprobs)
      kidx = sample.int(length(kvals), size = 1, prob = kprobs)
      qk = kprobs[kidx]
      k = kvals[kidx]
      # 1-3. cut unif spanning tree
      if(k > 1){
        sptree_cutted = igraph::make_undirected_graph(t(edge_list[-sample.int(n-1, size = k-1),,drop = F]), n = n)
        z = igraph::components(sptree_cutted)$membership
      }else{ # k=1
        z = rep(1, n)
      }
      Gquotient = make_quotient_graph(g0, z, return.igraph = T)
      # q(z), except product of tau(G_j) terms and tau(G) term 
      logq = log(qk)- lchoose(n-1, k-1) + graphppm::nsptrees_igraph(Gquotient, log = T)
      
      # p(z), where cohesions divided by tau(G_j) term
      logp = 0
      for(j in 1:k){ # TODO: avoid for loop?
        idx = which(z==j)
        #browser()
        logp = logp + logcohesion(A = A0[idx, idx, drop = F], 
                                 deg = deg0[idx], 
                                 cohesion_param = cohesion_param,
                                 minus_logtau = T) # log cohesionft / tau(G_j)
      }
      # p(z)/q(z)
      logweight = logp - logq 
      save_logweight[i] = logweight
      save_k[i] = k
      
      # h(rho) = 1(z_i = z_j), i.e. pairwise similiarity matrix 
      graphppm:::update_uppertri_logpsm(logsimmatrix, z, n, logweight)
      graphppm:::update_uppertri_logpsm(logsimmatrixsq, z, n, 2*logweight)
      
      if(i %%1000 == 0) cat(paste(i,"th iteration done\n"))
    }
  }else{# parallel settings
  }
  
  
  #logsumpoverq = matrixStats::logSumExp(savelogp - savelogq)
  logweightsum = logsimmatrix[1,1]
  logweightsqsum = logsimmatrixsq[1,1]
  
  # distribution of k 
  logdenom = matrixStats::logSumExp(save_logweight)
  lognumer = save_logweight
  ncluster = as.factor(save_k)
  kpmf = by(lognumer, ncluster, function(x) exp(matrixStats::logSumExp(x - logdenom)))
  attr(kpmf, "call") = NULL
  kpmf = as.table(kpmf)
  
  # pairwise
  simmatrix = exp(logsimmatrix - logweightsum)
  simmatrix = (simmatrix + t(simmatrix))
  diag(simmatrix) <- diag(simmatrix)/2
  
  # effective sample size
  effsize = exp(2*logsimmatrix - logsimmatrixsq)
  effsize[lower.tri(effsize)] = t(effsize)[lower.tri(effsize)]
  
  out = list()
  #out$ncluster = ncluster
  #out$ncluster_sd = ncluster_sd
  #out$ncluster_effsize = ncluster_effsize
  out$kpmf = kpmf
  out$simmatrix = simmatrix
  out$simmatrix_effsize = effsize
  out$save_logweight = save_logweight
  out$save_k = save_k
  return(out)
}

# example
# 


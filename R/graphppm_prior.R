#' GraphPPM prior analysis using MCMC
#' 
#' @param graph0 *igraph* igraph object, base graph
#' @param logcohesion *function* logcohesion function 
#' @param cohesion_param *list* cohesion function parameters
#' @param nsave *int* number of samples
#' @param nburn *int* number of burn-in iterations
#' @param nthin *int* number of thin-in rate. Total iterations: nburn + nsave*nthin
#' @param z_init *vector&lt;int&gt; (n)*, initial partition
#'
#' @return
#' @import igraph salso 
#' @export
#'
#' @examples
graphppm_prior <- function(graph0, logcohesion, cohesion_param = NULL, 
                           nsave = 1000, nburn = 1000, nthin = 1, 
                           z_init = NULL){
  niter = nburn + nthin*nsave 
  
  # input check for graph0
  if(!igraph::is.igraph(graph0)) stop("graph must be an igraph object")
  E0 = igraph::as_edgelist(graph0, names =F)
  n = igraph::vcount(graph0)
  m = igraph::ecount(graph0)
  rownames(E0) = 1:m
  A0 = igraph::as_adj(graph0, sparse = F, names = F) # sparse = F, lookup speed is faster in non-sparse form
  g0 = igraph::graph_from_edgelist(E0, directed = F)
  if(igraph::components(g0)$no > 1) stop("base graph must be connected")
  #adjlist0 = lapply(igraph::as_adj_list(g0), as.integer)
  deg0 = igraph::degree(g0) 
  V(g0)$vid = 1:n
  E(g0)$eid = 1:m
  
  # initialize patition, default one-block cluster
  if(is.null(z_init)) z_init = rep(1, n)
  z = partition_sort(z_init)
  k = length(unique(z))
  # check contiguity and calculate initial cohesion 
  for(j in 1:k){
    Gj = igraph::induced_subgraph(g0, vids = which(z==j))
    if(igraph::components(Gj)$no > 1) stop("z_init provided is invaild. Provide vaild contiguous clustering")
  }
  
  # initialize spanning tree
  # sptree is  (n-1) x 3 matrix with last column corresponds to is.crossing status (1 if crossing edge, 0 if not)
  sptree = runifsptree_given_z(g0, z, return.igraph = F)
  
  # save object
  counter = 0
  move_acc = numeric(4)
  move_cnt = numeric(4)
  save_z = matrix(0L, nrow = nsave, ncol = n)
  pd = 0.1

  
  # for loop
  for(iter in 1:niter){
    
    ## Step 1 ----------------
    if(k == 1) {pa = 0.9; pb = 0; pc = 0 # only split step allowed
    } else if(k == n) {pa = 0; pb = 0.9; pc = 0 # death/change/hyper allowed
    } else {pa = 0.3; pb = 0.3; pc = 0.3}
    move = sample(4, 1, prob = c(pa, pb, pc, pd))
    move_cnt[move] = move_cnt[move] + 1
    
    #cat(iter, move, log_like,'\n', sep = ' ')
    if(move == 1) { 
      
      ### Step 1(a), split -------
      
      # split cluster by choosing a within-cluster edge
      if(sum(sptree[,"iscrossing"]==0) > 1){
        cutidx = sample(which(sptree[,"iscrossing"]==0), size = 1)
      }else{
        cutidx = which(sptree[,"iscrossing"]==0)
      }
      
      # choose one cut_edge uniformly at random among non-crossing edges
      cutidx = sample(which(sptree[,"iscrossing"]==0), size = 1)
      # identify splitting cluster
      splitting_j = z[sptree[cutidx,1:2]]
      if(splitting_j[1] != splitting_j[2]) stop("wrong splitting edge")
      splitting_j = splitting_j[1]
      idx_j = which(z==splitting_j) # splitting cluster's vid
      
      treecutidx = c(cutidx, which(sptree[,"iscrossing"]==1)) # including already cutted edges
      # proposed split
      z_star = igraph::components(igraph::make_undirected_graph(t(sptree[-treecutidx,1:2]), n = n))$membership
      # identify splitted cluster
      splitted_j1j2 = z_star[sptree[cutidx,1:2]]
      
      idxstar1 = which(z_star==splitted_j1j2[1])
      idxstar2 = which(z_star==splitted_j1j2[2])
      
      loglik_ratio = logcohesion(A0[idxstar1,idxstar1, drop = F], deg0[idxstar1], cohesion_param) +
                     logcohesion(A0[idxstar2,idxstar2, drop = F], deg0[idxstar2], cohesion_param) - 
                     logcohesion(A0[idx_j, idx_j, drop = F], deg0[idx_j], cohesion_param) + # can be precalculated
                     nsptrees_igraph(make_quotient_graph(g0, z, return.igraph = T), log = T) -  # can be precalculated
                     nsptrees_igraph(make_quotient_graph(g0, z_star, return.igraph = T), log = T)
      
      if(k == n-1) {
        pb_new = 0.9
      } else {pb_new = 0.3}
      
      logprop_ratio = log(pb_new) - log(pa) - log(k) + log(n-k)
      
      #acceptance probability
      logacc_prob = min(0, loglik_ratio + logprop_ratio)
      
      if(runif(1) < exp(logacc_prob)){
        # accept
        move_acc[1] = move_acc[1] + 1
        z = partition_sort(z_star)
        k = k + 1
        sptree[cutidx, "iscrossing"] = 1
      }
    }
    
    if(move == 2) { 
      ### Step 1(b), merge-------------------
      # merge two existing clusters
      # Given spanning tree, merge cluster by choosing a cross-cluster edge
      if(sum(sptree[,"iscrossing"]==1) > 1){
        pasteidx = sample(which(sptree[,"iscrossing"]==1), size = 1)
      }else{
        pasteidx = which(sptree[,"iscrossing"]==1)
      }
      # identify merging cluster
      merging_j1j2 = z[sptree[pasteidx,1:2]]
      idx_j1 = which(z==merging_j1j2[1])
      idx_j2 = which(z==merging_j1j2[2])
      
      # simply combine idx_j1 and idx_j2 to merge
      z_star = z
      z_star[idx_j1] = min(merging_j1j2)
      z_star[idx_j2] = min(merging_j1j2)
      z_star = partition_sort(z_star)
      idx_jstar = union(idx_j1, idx_j2)
      
      loglik_ratio = logcohesion(A0[idx_jstar, idx_jstar, drop = F], deg0[idx_jstar], cohesion_param) - # can be precalculated
        logcohesion(A0[idx_j1, idx_j1, drop = F], deg0[idx_j1], cohesion_param) - 
        logcohesion(A0[idx_j2, idx_j2, drop = F], deg0[idx_j2], cohesion_param) + 
        nsptrees_igraph(make_quotient_graph(g0, z, return.igraph = T), log = T) -  
        nsptrees_igraph(make_quotient_graph(g0, z_star, return.igraph = T), log = T)
      
      # # compute log-proposal ratio
      if(k == 2) {pa_new = 0.9
      }else {pa_new = 0.3}
      
      logprop_ratio = log(pa_new) - log(pb) - log(n-k+1) + log(k-1)
      
      #acceptance probability
      logacc_prob = min(0, loglik_ratio + logprop_ratio)
      
      if(runif(1) < exp(logacc_prob)){
        # accept
        move_acc[2] = move_acc[2] + 1
        z = partition_sort(z_star)
        k = k - 1
        sptree[pasteidx, "iscrossing"] = 0
      }
    }
    if(move == 3){ 
      ### step 1(c), change ------
      
      # merge cluster by choosing a cross-cluster edge
      if(sum(sptree[,"iscrossing"]==1) > 1){
        pasteidx = sample(which(sptree[,"iscrossing"]==1), size = 1)
      }else{
        pasteidx = which(sptree[,"iscrossing"]==1)
      }
      
      # identify merging cluster
      merging_j1j2 = z[sptree[pasteidx,1:2]]
      idx_j1 = which(z==merging_j1j2[1])
      idx_j2 = which(z==merging_j1j2[2])
      
      # simply combine idx_j1 and idx_j2 to merge
      z_star = z
      z_star[idx_j1] = min(merging_j1j2)
      z_star[idx_j2] = min(merging_j1j2)
      z_star = partition_sort(z_star)
      idx_jstar = union(idx_j1, idx_j2)
      
      sptree_star = sptree
      sptree_star[pasteidx, "iscrossing"] = 0
      
      loglik_ratio1 = logcohesion(A0[idx_jstar, idx_jstar, drop = F], deg0[idx_jstar], cohesion_param) - # can be precalculated
        logcohesion(A0[idx_j1, idx_j1, drop = F], deg0[idx_j1], cohesion_param) - 
        logcohesion(A0[idx_j2, idx_j2, drop = F], deg0[idx_j2], cohesion_param) + 
        nsptrees_igraph(make_quotient_graph(g0, z, return.igraph = T), log = T) -  
        nsptrees_igraph(make_quotient_graph(g0, z_star, return.igraph = T), log = T)
      
      # split cluster by choosing a within-cluster edge
      if(sum(sptree_star[,"iscrossing"]==0) > 1){
        cutidx = sample(which(sptree_star[,"iscrossing"]==0), size = 1)
      }else{
        cutidx = which(sptree_star[,"iscrossing"]==0)
      }
      
      # identify splitting cluster
      splitting_j = z_star[sptree_star[cutidx,1:2]]
      if(splitting_j[1] != splitting_j[2]) stop("wrong splitting edge")
      splitting_j = splitting_j[1]
      idx_j = which(z_star==splitting_j) # splitting cluster's vid
      
      treecutidx = c(cutidx, which(sptree_star[,"iscrossing"]==1)) # including already cutted edges
      # proposed split
      z_starstar = igraph::components(igraph::make_undirected_graph(t(sptree_star[-treecutidx,1:2]), n = n))$membership
      # identify splitted cluster
      splitted_j1j2 = z_starstar[sptree_star[cutidx,1:2]]
      
      idxstar1 = which(z_starstar==splitted_j1j2[1])
      idxstar2 = which(z_starstar==splitted_j1j2[2])
      
      loglik_ratio2 = logcohesion(A0[idxstar1,idxstar1, drop = F], deg0[idxstar1], cohesion_param) +
        logcohesion(A0[idxstar2,idxstar2, drop = F], deg0[idxstar2], cohesion_param) - 
        logcohesion(A0[idx_j, idx_j, drop = F], deg0[idx_j], cohesion_param) + # can be precalculated
        nsptrees_igraph(make_quotient_graph(g0, z_star, return.igraph = T), log = T) -  # can be precalculated
        nsptrees_igraph(make_quotient_graph(g0, z_starstar, return.igraph = T), log = T)
      
      logprop_ratio = 0
      #acceptance probability
      logacc_prob = min(0, loglik_ratio1 + loglik_ratio2 + logprop_ratio)
      
      if(runif(1) < exp(logacc_prob)){
        # accept
        move_acc[3] = move_acc[3] + 1
        z = partition_sort(z_starstar)
        sptree_star[cutidx, "iscrossing"] = 1
        sptree = sptree_star
      }
    }
    if(move == 4) { 
      ### Step 1(d), hyper ####
      # update UST
      move_acc[4] = move_acc[4] + 1
      #browser()
      sptree = runifsptree_given_z(g0, z, return.igraph = F)
    }
    
    # save loop 
    
    if(iter > nburn & (iter - nburn) %% nthin == 0) {
      counter = counter + 1
      save_z[counter,] <- z
    }
    # if(verbose){
    #   if(imcmc %% 1000 == 0){
    #     cat("iter",imcmc,",",paste(z)," k:",k,"\n")
    #     V(graph0)$color <- z
    #     plot(graph0, layout = coords, main = paste("iteration",iter))
    #     readline(prompt="Press [enter] to continue")
    #   }
    # }
    if(iter %% 1000 == 0) cat(paste("iteration",iter,"done\n"))
  }
  
  out = list()
  out$save_z = save_z
  out$kpmf = table(apply(save_z, 1, function(x) (length(unique(x)))))/nsave
  out$simmatrix = salso::psm(save_z)
  out$move_cnt = move_cnt
  out$move_acc = move_acc
  return(out)
}

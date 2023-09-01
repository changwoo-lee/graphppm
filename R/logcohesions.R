
# modularity 
# log(c) where c = alpha * c0^psi, and
# c0(G_j, \partial G_j) = exp( resolution_parameter * (2|E_j| + partial G_j)^2/(2m)^2 - (2*|E_j|)/(2m)
# default resolution_parameter = 1
# REQUIRE argument m, the number of edges in the original graph
#
logcohesion_modularity <- function(A, deg, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  resolution_parameter = cohesion_param$resolution_parameter
  m = cohesion_param$m
  
  logc0 = sum(A)/(2*m) - resolution_parameter * (sum(deg)/(2*m))^2
  if(minus_logtau){
    logtau = graphppm::nsptrees(A, log = T) 
    out = log(alpha) + psi*logc0 - logtau
  }else{
    out = log(alpha) + psi*logc0
  }
  attr(out,"logc0") = logc0
  return(out)
}

# constant Potts model 
# log(c) where c = alpha * c0^psi, and
# c0(G_j, \partial G_j) = resolution * |V_j|^2 - (2*|E_j|)
# default resolution = quantile(degree(g))[2] / (degree(g) - 1) ? (See ?igraph::cluster_leiden)
# require argument m, the number of edges in the original graph

logcohesion_CPM <- function(A, deg, cohesion_param, minus_logtau = T){
  
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  resolution_parameter = cohesion_param$resolution_parameter
  
  logc0 = sum(A) - resolution_parameter * ncol(A)^2
  if(minus_logtau){
    logtau = graphppm::nsptrees(A, log = T) 
    out = log(alpha) + psi*logc0 - logtau
  }else{
    out = log(alpha) + psi*logc0
  }
  attr(out,"logc0") = logc0
  return(out)
}



# 2. spanning tree based cohesions ##
#
# log(c) where c(G_j) = alpha * (tau(G_j))^psi
# 
# if div_by_nsptree = T, return log(c) - log(tau)

logcohesion_tau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  
  if(minus_logtau && (psi == 1)){  # no need to calculate tau(G_j)
    out = log(alpha)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = logtau
    if(minus_logtau){
      out = log(alpha) + psi*logc0 - logtau
    }else{
      out = log(alpha) + psi*logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}



# log(c) where c(G_j) = alpha * (|E_j|*tau(G_j))^psi
# 
# if div_by_nsptree = T, return log(c) - log(tau)

logcohesion_etau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  
  nj = ncol(A)
  if(nj==1){ ej = 1 }else{  ej = sum(A)/2 }
  
  if(minus_logtau && (psi == 1)){  # no need to calculate tau(G_j)
    out = log(alpha) + log(ej)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = log(sum(A)/2) + logtau
    if(minus_logtau){
      out = log(alpha) + psi*logc0 - logtau
    }else{
      out = log(alpha) + psi*logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}


# log(c) where c(G_j) = alpha * (|V_j|^(-1/2)|E_j|*tau(G_j))^psi
# 
# if div_by_nsptree = T, return log(c) - log(tau)

logcohesion_vetau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  
  nj = ncol(A)
  if(nj==1){ ej = 1 }else{  ej = sum(A)/2 }
  
  if(minus_logtau && (psi == 1)){  # no need to calculate tau(G_j)
    out = log(alpha) -0.5*log(nj) + log(ej)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = -0.5*log(nj) + log(ej) + logtau
    if(minus_logtau){
      out = log(alpha) + psi*logc0 - logtau
    }else{
      out = log(alpha) + psi*logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}



### misc ###############
# log(c) where c(G_j) = alpha * Gamma(n_j)^psi

logcohesion_gamma <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  psi = cohesion_param$psi
  
  nj = ncol(A)
  logc0 = lgamma(nj)
  if(minus_logtau){
    logtau = graphppm::nsptrees(A, log = T) 
    out = log(alpha) + psi*logc0 - logtau
  }else{
    out = log(alpha) + psi*logc0
  }
  attr(out,"logc0") = logc0
  return(out)
}




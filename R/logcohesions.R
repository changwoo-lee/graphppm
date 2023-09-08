#' Modularity cohesion function (loss-based)
#'
#' Return c(G_j, \partial G_j) = exp(-psi*L)
#' where L = resolution * (2|E_j| + partial G_j)^2/(2m)^2 - (2*|E_j|)/(2m)
#' 
#'
#' @param A  *matrix&lt;int&gt; (nj,nj)*, adjacency matrix, can be a sparse format.
#' @param deg *vec&lt;int&gt*, vector of node degrees with length nj
#' @param cohesion_param *list*, must include "psi", "resolution", and "m" parameters. 
#' @param minus_logtau *logical* if true, return log c(Gj) - tau(Gj) 
#'
#' @return
#' @export
#'
#' @examples
logcohesion_modularity <- function(A, deg, cohesion_param, minus_logtau = T){
  
  psi = cohesion_param$psi
  resolution = cohesion_param$resolution
  m = cohesion_param$m
  
  logc0 = sum(A)/(2*m) - resolution * (sum(deg)/(2*m))^2
  if(minus_logtau){
    logtau = graphppm::nsptrees(A, log = T) 
    out =  psi*logc0 - logtau
  }else{
    out =  psi*logc0
  }
  attr(out,"logc0") = logc0
  return(out)
}

 
#' Cohesion based on constant Potts model (loss-based)
#'
#' log(c) where c = alpha * c0^psi, and
#' c0(G_j, \partial G_j) = resolution * |V_j|^2 - (2*|E_j|)
#' default resolution = quantile(degree(g))[2] / (degree(g) - 1) ? (See ?igraph::cluster_leiden)
#' require argument m, the number of edges in the original graph
#'
#' @param A 
#' @param deg 
#' @param cohesion_param 
#' @param minus_logtau 
#'
#' @return
#' @export
#'
#' @examples
logcohesion_CPM <- function(A, deg, cohesion_param, minus_logtau = T){
  
  psi = cohesion_param$psi
  resolution = cohesion_param$resolution
  
  logc0 = sum(A) - resolution * ncol(A)^2
  if(minus_logtau){
    logtau = graphppm::nsptrees(A, log = T) 
    out = psi*logc0 - logtau
  }else{
    out = psi*logc0
  }
  attr(out,"logc0") = logc0
  return(out)
}


#' spanning tree cohesion (tree-based)
#'
#' @param A 
#' @param deg 
#' @param cohesion_param 
#' @param minus_logtau 
#'
#' @return
#' @export
#'
#' @examples
logcohesion_tau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  
  if(minus_logtau){  # no need to calculate tau(G_j)
    out = log(alpha)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = logtau
    if(minus_logtau){
      out = log(alpha) + logc0 - logtau
    }else{
      out = log(alpha) + logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}



#' spanning tree cohesion (tree-based)
#'
#' @param A 
#' @param deg 
#' @param cohesion_param 
#' @param minus_logtau 
#'
#' @return
#' @export
#'
#' @examples
logcohesion_etau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  
  nj = ncol(A)
  if(nj==1){ ej = 1 }else{  ej = sum(A)/2 }
  
  if(minus_logtau){  # no need to calculate tau(G_j)
    out = log(alpha) + log(ej)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = log(sum(A)/2) + logtau
    if(minus_logtau){
      out = log(alpha) + logc0 - logtau
    }else{
      out = log(alpha) + logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}

#' spanning tree cohesion (tree-based)
#'
#' @param A 
#' @param deg 
#' @param cohesion_param 
#' @param minus_logtau 
#'
#' @return
#' @export
#'
#' @examples
logcohesion_vetau <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  
  nj = ncol(A)
  if(nj==1){ ej = 1 }else{  ej = sum(A)/2 }
  
  if(minus_logtau){  # no need to calculate tau(G_j)
    out = log(alpha) -0.5*log(nj) + log(ej)
    attr(out,"logc0") = NA
  }else{ # need to calculate tau(G_j)
    logtau = graphppm::nsptrees(A, log = T) 
    logc0 = -0.5*log(nj) + log(ej) + logtau
    if(minus_logtau){
      out = log(alpha) + logc0 - logtau
    }else{
      out = log(alpha) + logc0
    }
    attr(out,"logc0") = logc0
  }
  return(out)
}


#' misc cohesion 
#'
#' @param A 
#' @param deg 
#' @param cohesion_param 
#' @param minus_logtau 
#'
#' @return
#' @export
#'
#' @examples
logcohesion_gamma <- function(A, deg = NULL, cohesion_param, minus_logtau = T){
  alpha = cohesion_param$alpha
  
  nj = ncol(A)
  logc0 = lgamma(nj)
  if(minus_logtau){
    logtau = graphppm::nsptrees(A, log = T) 
    out = log(alpha) + logc0 - logtau
  }else{
    out = log(alpha) + logc0
  }
  attr(out,"logc0") = logc0
  return(out)
}




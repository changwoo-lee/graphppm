#' Draw a uniform spanning tree using Wilson's algorithm
#'
#' Equivalent algorithm available in igraph package igraph::sample_spanning_tree(g) but this one should be much faster. 
#'
#' @param g0 *igraph object*, Graph.
#' @param root *int*, (optional) Root node. Recommended as a node with max degree
#' @param ordering *vector&lt;int&gt; (n)*, (optional) Ordering of nodes. Recommended as a breadth-first search ordering at the root. 
#' @param return.igraph *logical*, return igraph, otherwise edgelist
#'
#' @return *matrix&lt;int&gt; (n-1,2)*, Edge list of the sampled tree.
#' @import Rcpp igraph
#' @export
#' @useDynLib graphppm
#' 
#' @examples
#' g = igraph::make_graph("Zachary")
#' root = which.max(igraph::degree(g))
#' ordering = as.integer(igraph::bfs(g, root)$order)
#' ust = runifsptree(g, root, ordering)
#' 
#' 
runifsptree <- function(g0, root = NULL, ordering = NULL, maxiter = NULL, return.igraph = F){
  E0 = igraph::as_edgelist(g0, names = F)
  n = igraph::vcount(g0)
  if(is.null(root)) root = 1L
  mode(root) = "integer"
  if(is.null(ordering)) ordering = 1:n
  mode(ordering) = "integer"
  if(is.null(maxiter)) maxiter = max(100, 10*n^3)
  
  edge_list = graphppm:::runcppWilson(n, E0, root, ordering, maxiter)
  if(return.igraph) return(igraph::graph_from_edgelist(edge_list, directed = F))
  return(edge_list)
}



#' Draw a uniform spanning tree conditional on a partition using Wilson's algorithm
#'
#' @param g0 *igraph object*, Graph.
#' @param z *vector&lt;int&gt; (n)*, partition
#' @param return.igraph *logical*, return igraph, otherwise edgelist with THREE columns where last column corresponds to crossing status
#'
#' @return *igraph or matrix&lt;int&gt; (n-1,3)*, 
#' @import Rcpp igraph 
#' @export
#' @useDynLib graphppm
#' 
#' @examples
#' g0 = igraph::make_graph("Zachary")
#' z = c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
#' runifsptree_given_z(g0, z)
#' 
runifsptree_given_z <- function(g0, z, return.igraph = F){
  E0 = igraph::as_edgelist(g0, names = F)
  n = igraph::vcount(g0)
  V(g0)$vid = 1:n
  E(g0)$eid = 1:ecount(g0)
  
  # input check for z
  if(!partition_vaild(z)){
    warning("z is not a vaild partition, converted using partition_sort")
    z = partition_sort(z)
  } 
  k = max(z)
  
  # sample ust for each subgraphs
  edgelist = list()
  for(j in 1:k){
    gj = igraph::induced_subgraph(g0, vids = which(z==j))
    if(!is.connected(gj)) stop("induced subgraph is not connected, invaild z")
    nj = vcount(gj)
    if(nj == 1){
      edgelist[[j]] = matrix(0, 0, ncol = 3)# third column: is.crossing status 
      next;
    }
    vid = V(gj)$vid
    edges = graphppm:::runcppWilson(nj, as_edgelist(gj, names = F), 1, 1:nj, max(100,nj^3))
    edges[,1] = vid[edges[,1]]
    edges[,2] = vid[edges[,2]]
    edgelist[[j]] = cbind(edges, 0) # third column: is.crossing status 
  }
  # sample ust from quotient graph, tracking eids
  qgraph = make_quotient_graph(g0, z, return.igraph = T)
  if(vcount(qgraph) > 1){
    qtree = igraph::sample_spanning_tree(qgraph) # using igraph::sample_spanning_tree to track eid
    edgelist[[k+1]] = cbind(E0[qtree$eid,, drop = F], 1) # third column: is.crossing status
  }
  # list to edgelist
  sptree = do.call(rbind, edgelist)
  # third column
  colnames(sptree) = c("v1","v2","iscrossing")
  # final check
  if( nrow(sptree)!= n-1) stop("error, not a vaild spanning tree")
  
  if(return.igraph){
    sptree.igraph = graph_from_edgelist(sptree[,1:2], directed = F)
    E(sptree.igraph)$iscrossing = sptree[,3]
    return(sptree.igraph)
  }else{
    return(sptree)
  }
}




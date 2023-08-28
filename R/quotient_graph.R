#' Get edge crossing status given a partition
#' 
#' Get logical vector that represents whether an edge is within a cluster or crossing two clusters.
#' 
#' 
#' @param z *vector&lt;int&gt; (n)*, partition, cluster membership vector.
#' @param E0 *matrix&lt;int&gt; (m,2)*, edge list of the graph.
#'
#' @return *vector&lt;logical&gt; (m)*, TRUE if cross-cluster edge, FALSE if between-cluster edge.
#' @export
#'
#' @examples
#' g0 = igraph::make_graph("Zachary")
#' z = c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
#' E0 = igraph::as_edgelist(g0)
#' rownames(E0) = 1:igraph::ecount(g0)
#' is.crossing(E0, z)
is.crossing <- function(E0, z) {
  membership_head = z[E0[, 1]]
  membership_tail = z[E0[, 2]]
  iscrossing = rep(FALSE, nrow(E0))
  iscrossing[membership_head != membership_tail] = TRUE # crossing
  names(iscrossing) = rownames(E0)
  return(iscrossing)
}

#' Make the quotient graph modulo a partition
#' 
#' Make G/z, a quotient graph of G modulo a partition z. 
#' The nodes of G/z corresponds to clusters (blocks) of a partition, and edges corresponds to edges between clusters.   
#' Resulting quotient graph is often a multigraph, i.e. can have more than one edges between two nodes. 
#' 
#' @param z *vector&lt;int&gt; (n)*, cluster membership
#' @param E0 *matrix&lt;int&gt; (m,2)*, edge list of the graph
#' @param return.igraph *logical*, return igraph, otherwise edgelist. Note: igraph output has eid from the original graph 
#'
#' @return *edge list (with 2 columns) or igraph object*, quotient graph modulo a partition. Edge list output contains rownames as an edge id.
#' @import igraph
#' @export
#'
#' @examples
#' g0 = igraph::make_graph("Zachary")
#' z = c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
#' make_quotient_graph(g0, z)
make_quotient_graph <- function(g0, z, return.igraph = F){
  E0 = as_edgelist(g0)
  #rownames(E0) = 1:ecount(g0) # eid
  iscrossing = is.crossing(E0, z)
  # special case
  if(sum(iscrossing)==0){
    if(!return.igraph){
      return(igraph::as_edgelist(igraph::make_empty_graph(1)))
    }else{
      return(igraph::make_empty_graph(1)) # one-block
    }
  }
  edgelist_quotient = cbind(z[E0[iscrossing,1]],z[E0[iscrossing,2]])
  rownames(edgelist_quotient) = which(iscrossing) # eid of crossing edges (WARN: this is saved as char, not num)
  if(!return.igraph){
    return(edgelist_quotient)
  }else{
    qgraph = igraph::graph_from_edgelist(edgelist_quotient, directed = F)
    E(qgraph)$eid = which(iscrossing)
    return(qgraph)
  }
}

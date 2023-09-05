#' The number of spanning trees of a graph
#'
#' Calculate the number of spanning trees of a graph, using adjacency matrix or igraph object
#'
#' @param A *matrix&lt;int&gt; (n,n)*, Adjacency matrix, can be a sparse format.
#' @param log *logical*, If TRUE (default), returns log value. 
#'
#' @return *double*, The (log of) number of spanning trees 
#' @import igraph Matrix
#' @export
#'
#' @examples
#' # Example 1
#' g = igraph::make_graph("Diamond")
#' A = igraph::as_adj(g, sparse = F)
#' nsptrees(A, log = F)
#' nsptrees_igraph(g, log = F)
#' 
#' # Example 2
#' g = igraph::make_graph("Zachary")
#' A = igraph::as_adj(g, sparse = T)
#' nsptrees(A, log = T)
#' nsptrees_igraph(g, log = T, sparse = T)
nsptrees <- function(A, log = T){
  if(dim(A)[1]==1){# singleton
    out = 0
  }else{
    if(is(A,"sparseMatrix")){
      L = Matrix::Diagonal(x = Matrix::colSums(A)) - A
      L = as(L[-1,-1, drop = F], "symmetricMatrix") # make symmetric improves speed a lot
      #R = chol(L)
      #out = 2*sum(log(diag(R)))
      out = as.numeric(Matrix::determinant(L, logarithm = T)$modulus) # Matrix::determinant  
    }else{
      deg = colSums(A)
      L = -A 
      diag(L) <- deg
      out = as.numeric(determinant(L[-1,-1, drop = F], logarithm = T)$modulus)
    }
  } 
  if(!log) out = round(exp(out),0)
  out
}

#' The number of spanning trees of a graph
#'
#' @param g *igraph object*, Graph
#' @param log *logical*, If TRUE (default T), returns log value. 
#' @param sparse *logical*, If TRUE (default F), uses sparse matrix calculation. 
#'
#' @import igraph Matrix
#' @export
#' @describeIn nsptrees Input is replaced with an igraph object g.
#' 
nsptrees_igraph <- function(g, log = T, sparse = F){
  if(igraph::vcount(g)==1){# singleton
    out = 0
  }else{
    if(sparse){
      L = igraph::graph.laplacian(g, weights = NA, sparse = T)
      L = as(L[-1,-1, drop = F], "symmetricMatrix") # make symmetric improves speed a lot
      out = as.numeric(Matrix::determinant(L, logarithm = T)$modulus) # Matrix::determinant  
    }else{
      L = igraph::graph.laplacian(g, weights = NA, sparse = F)
      out = as.numeric(determinant(L[-1,-1, drop = F], logarithm = T)$modulus)
    }
  } 
  if(!log) out =round(exp(out),0)
  out
}





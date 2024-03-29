% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quotient_graph.R
\name{is.crossing}
\alias{is.crossing}
\title{Get edge crossing status given a partition}
\usage{
is.crossing(E0, z)
}
\arguments{
\item{E0}{\emph{matrix<int> (m,2)}, edge list of the graph.}

\item{z}{\emph{vector<int> (n)}, partition, cluster membership vector.}
}
\value{
\emph{vector<logical> (m)}, TRUE if cross-cluster edge, FALSE if between-cluster edge.
}
\description{
Get logical vector that represents whether an edge is within a cluster or crossing two clusters.
}
\examples{
g0 = igraph::make_graph("Zachary")
z = c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)
E0 = igraph::as_edgelist(g0)
rownames(E0) = 1:igraph::ecount(g0)
is.crossing(E0, z)
}

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;

// #include <boost/heap/fibonacci_heap.hpp>
#include "vertex.h"
#include "edge.h"
#include "graph_data_type.h"
#include "graph.h"


// [[Rcpp::export]]
void testVertexEdge() {
  v_ptr v1 = std::make_shared<Vertex>(0);
  v_ptr v2 = std::make_shared<Vertex>(1);
  Rcout << "Vertex 1: " << v1->getVid() << "\n";
  Rcout << "Vertex 2: " << v2->getVid() << "\n";
  
  Edge e1(0, v1, v2);
  Rcout << e1.print() << "\n";
  Rcout << "eid " << e1.getEid() << "\n";
  Rcout << "Endpoint 1: " << e1.getEndpoints().first->getVid() << "\n";
  Rcout << "Endpoint 2: " << e1.getEndpoints().second->getVid() << "\n";
  
  std::pair<const v_ptr, const v_ptr> endpoints(v1, v2);
  Edge e2(1, endpoints);
  Rcout << "eid " << e2.getEid() << "\n";
  Rcout << "Endpoint 1: " << e2.getEndpoints().first->getVid() << "\n";
  Rcout << "Endpoint 2: " << e2.getEndpoints().second->getVid() << "\n";
  Rcout << "Neighbor of endpoint 2: " << e2.getNeighbor(v2)->getVid() << "\n";
  
  // test hash functions
  Rcout << "Hash(v1) = " << VertexHashFunc()(v1) << "\n";
  Rcout << "Hash(v2) = " << VertexHashFunc()(v2) << "\n";
  Rcout << "Hash(e1) = " << EdgeHashFunc()(std::make_shared<Edge>(e1)) << "\n";
  Rcout << "Hash(e2) = " << EdgeHashFunc()(std::make_shared<Edge>(e2)) << "\n";
}

// [[Rcpp::export]]
void testGraph() {
  v_vec vertices;
  for(int i = 0; i < 5; ++i) {
    vertices.push_back(std::make_shared<Vertex>(4 - i));
    Rcout << vertices[i]->getVid() << "\n";
  }
  
  e_vec edges;
  edges.push_back(std::make_shared<Edge>(0, vertices[0], vertices[1]));
  edges.push_back(std::make_shared<Edge>(1, vertices[0], vertices[2]));
  edges.push_back(std::make_shared<Edge>(2, vertices[1], vertices[3]));
  edges.push_back(std::make_shared<Edge>(3, vertices[1], vertices[4]));
  
  // make graph
  Graph g(vertices, edges);
  Rcout << "vcount = " << g.vcount() << "\n";
  Rcout << "ecount = " << g.ecount() << "\n";
  
  Rcout << "Vertex ids: \n";
  for (const auto v : g.getVertices())
    Rcout << v->getVid() << "\n";
  
  const std::vector<v_vec> v_adj_list = g.getVAdjList();
  Rcout << "v_adj_list size = " << v_adj_list.size() << "\n";
  
  const std::vector<e_vec> e_adj_list = g.getEAdjList();
  Rcout << "e_adj_list size = " << e_adj_list.size() << "\n";
  
  // print adjacency list
  Rcout << "Vertex adjacency list: \n";
  for(unsigned int i = 0; i < v_adj_list.size(); ++i) {
    Rcout << g.getVertices()[i]->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(auto v : v_adj_list[i])
      Rcout << v->getVid() << ", ";
    Rcout << '\n';
  }
  
  // print edge adjacency list
  Rcout << "Edge adjacency list: \n";
  for(unsigned int i = 0; i < e_adj_list.size(); ++i) {
    v_ptr v = g.getVertices()[i];
    Rcout << v->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(const auto edge : e_adj_list[v->getVid()])
      Rcout << edge->getEid() << ", ";
    Rcout << "\n";
  }
}

// [[Rcpp::export]]
void testSubgraph() {
  v_vec vertices;
  for(int i = 0; i < 6; ++i) {
    vertices.push_back(std::make_shared<Vertex>(i));
    // Rcout << vertices[i]->getVid() << "\n";
  }
  
  e_vec edges;
  edges.push_back(std::make_shared<Edge>(0, vertices[0], vertices[1]));
  edges.push_back(std::make_shared<Edge>(1, vertices[0], vertices[2]));
  edges.push_back(std::make_shared<Edge>(2, vertices[1], vertices[3]));
  edges.push_back(std::make_shared<Edge>(3, vertices[1], vertices[4]));
  edges.push_back(std::make_shared<Edge>(4, vertices[2], vertices[5]));
  
  // make spanning tree
  Graph st(vertices, edges);
  // make an induced spanning forest subgraph
  v_unordered_set v_sf {vertices[1], vertices[2], vertices[3], vertices[4], vertices[5]};
  Subgraph sf(st, v_sf);
  
  Rcout << "Spanning tree vcount = " << st.vcount() << "\n";
  Rcout << "Spanning tree ecount = " << st.ecount() << "\n";
  Rcout << "Spanning forest vcount = " << sf.vcount() << "\n";
  Rcout << "Spanning forest ecount = " << sf.ecount() << "\n";
  
  // print spanning tree adjacency list
  Rcout << "Spanning tree vertex adjacency list: \n";
  const std::vector<v_vec> v_adj_list_st = st.getVAdjList();
  for(unsigned int i = 0; i < v_adj_list_st.size(); ++i) {
    Rcout << st.getVertices()[i]->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(auto v : v_adj_list_st[i])
      Rcout << v->getVid() << ", ";
    Rcout << '\n';
  }
  
  // print spanning forest adjacency list
  Rcout << "Spanning forest vertex adjacency list: \n";
  const v_unordered_map<v_set> v_adj_map_sf = sf.getVAdjMap();
  for(const auto &adj_pair : v_adj_map_sf) {
    Rcout << adj_pair.first->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(const auto v : adj_pair.second)
      Rcout << v->getVid() << ", ";
    Rcout << '\n';
  }
  
  // find connected components in spanning tree
  Rcout << "Spanning tree connected components after removing edge (0, 1): \n";
  v_unordered_set v_st {vertices[0], vertices[1], vertices[2], vertices[3], vertices[4], vertices[5]};
  Subgraph st_2(st, v_st);
  v_unordered_map<unsigned int> components_st = st_2.getComponents(edges[0]);
  v_unordered_set component0_st;
  st_2.getTreeComponents(edges[0], component0_st);
  for(const auto v_component_pair : components_st)
    Rcout << v_component_pair.first->getVid() << " -> " << v_component_pair.second << "\n";
  Rcout << "Vertices in the component containing vertex " << edges[0]->getEndpoint(0)->getVid() << " :\n";
  for (const auto &v : component0_st)
    Rcout << v->getVid() << "  ";
  Rcout << "\n";
  
  // find connected components in spanning forest
  Rcout << "Spanning forest connected components after removing edge (1, 3): \n";
  v_unordered_map<unsigned int> components_sf = sf.getComponents(edges[2]);
  for(const auto v_component_pair : components_sf)
    Rcout << v_component_pair.first->getVid() << " -> " << v_component_pair.second << "\n";
  
}

// [[Rcpp::export]]
void testMST() {
  // make graph
  v_vec vertices;
  for(int i = 0; i < 6; ++i) {
    vertices.push_back(std::make_shared<Vertex>(i));
  }
  
  e_vec edges;
  edges.push_back(std::make_shared<Edge>(0, vertices[0], vertices[1], 0.1));
  edges.push_back(std::make_shared<Edge>(1, vertices[0], vertices[2], 0.2));
  edges.push_back(std::make_shared<Edge>(2, vertices[1], vertices[3], 0.3));
  edges.push_back(std::make_shared<Edge>(3, vertices[1], vertices[4], 0.5));
  edges.push_back(std::make_shared<Edge>(4, vertices[2], vertices[5], 0.1));
  edges.push_back(std::make_shared<Edge>(5, vertices[2], vertices[4], 0.4));
  edges.push_back(std::make_shared<Edge>(6, vertices[3], vertices[4], 0.2));
  edges.push_back(std::make_shared<Edge>(7, vertices[4], vertices[5], 0.1));
  
  Graph g(vertices, edges);
  
  // get unrooted MST
  Graph mst;
  g.getMST(mst);
  
  Rcout << "Spanning tree vcount = " << mst.vcount() << "\n";
  Rcout << "Spanning tree ecount = " << mst.ecount() << "\n";
  
  Rcout << "Spanning tree edges: \n";
  for(auto e : mst.getEdges())
    Rcout << e->print() << "\n";
  
  // print unrooted spanning tree adjacency list
  Rcout << "Spanning tree vertex adjacency list: \n";
  const std::vector<v_vec> v_adj_list_mst = mst.getVAdjList();
  for(unsigned int i = 0; i < v_adj_list_mst.size(); ++i) {
    Rcout << mst.getVertices()[i]->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(auto v : v_adj_list_mst[i])
      Rcout << v->getVid() << ", ";
    Rcout << '\n';
  }
  
  // get connected components with edges[4] removed
  Rcout << "Clusters after edge (2 ,5) removed: \n";
  v_unordered_set v_mst {vertices[0], vertices[1], vertices[2], vertices[3],
                         vertices[4], vertices[5]};
  Subgraph mst_2(mst, v_mst);
  v_unordered_map<unsigned int> components_mst = mst_2.getComponents(edges[4]);
  for(const auto v_component_pair : components_mst)
    Rcout << v_component_pair.first->getVid() << " -> " << v_component_pair.second << "\n";
  
  v_unordered_set component0_mst;
  mst_2.getTreeComponents(edges[4], component0_mst);
  Rcout << "Vertices in the component containing vertex " << edges[4]->getEndpoint(0)->getVid() << " :\n";
  for (const auto &v : component0_mst)
    Rcout << v->getVid() << "  ";
  Rcout << "\n";
  
  // get rooted MST
  Graph mst_rooted;
  g.getMST(mst_rooted, 1000.0, true);
  
  Rcout << "Rooted spanning tree vcount = " << mst.vcount() << "\n";
  Rcout << "Rooted spanning tree ecount = " << mst.ecount() << "\n";
  
  Rcout << "Rooted spanning tree edges: \n";
  for(auto e : mst_rooted.getEdges()) {
    v_ptr v_head = e->getEndpoints().first;
    v_ptr v_tail = e->getEndpoints().second;
    Rcout << e->getEid() << ": " << v_head->getVid() << "->"
          << v_tail->getVid() << "\n";
  }
  
}

// [[Rcpp::export]]
void testMSF() {
  // make graph
  v_vec vertices;
  for(int i = 0; i < 7; ++i) {
    vertices.push_back(std::make_shared<Vertex>(i));
  }
  
  e_vec edges;
  edges.push_back(std::make_shared<Edge>(0, vertices[0], vertices[1], 0.1));
  edges.push_back(std::make_shared<Edge>(1, vertices[0], vertices[2], 0.2));
  edges.push_back(std::make_shared<Edge>(2, vertices[1], vertices[3], 0.3));
  edges.push_back(std::make_shared<Edge>(3, vertices[1], vertices[2], 0.5));
  edges.push_back(std::make_shared<Edge>(4, vertices[2], vertices[3], 0.1));
  edges.push_back(std::make_shared<Edge>(5, vertices[0], vertices[3], 0.4));
  edges.push_back(std::make_shared<Edge>(6, vertices[4], vertices[6], 0.2));
  edges.push_back(std::make_shared<Edge>(7, vertices[4], vertices[5], 0.1));
  edges.push_back(std::make_shared<Edge>(8, vertices[6], vertices[5], 0.1));
  
  Graph g(vertices, edges);
  
  // get minimal spanning forest
  Graph msf;
  g.getMST(msf);
  
  Rcout << "Spanning forest vcount = " << msf.vcount() << "\n";
  Rcout << "Spanning forest ecount = " << msf.ecount() << "\n";
  
  Rcout << "Spanning forest edges: \n";
  for(const auto e : msf.getEdges())
    Rcout << e->print() << "\n";
  
  // print spanning forest adjacency list
  Rcout << "Spanning forest vertex adjacency list: \n";
  const std::vector<v_vec> v_adj_list_msf = msf.getVAdjList();
  for(unsigned int i = 0; i < v_adj_list_msf.size(); ++i) {
    Rcout << msf.getVertices()[i]->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(auto v : v_adj_list_msf[i])
      Rcout << v->getVid() << ", ";
    Rcout << '\n';
  }
}

// [[Rcpp::export]]
void testRST() {
  /* test random spanning tree */
  
  // make graph
  v_vec vertices;
  for(int i = 0; i < 9; ++i) {
    vertices.push_back(std::make_shared<Vertex>(i));
  }
  
  e_vec edges;
  edges.push_back(std::make_shared<Edge>(0, vertices[0], vertices[1]));
  edges.push_back(std::make_shared<Edge>(1, vertices[0], vertices[2]));
  edges.push_back(std::make_shared<Edge>(2, vertices[1], vertices[3]));
  edges.push_back(std::make_shared<Edge>(3, vertices[1], vertices[4]));
  edges.push_back(std::make_shared<Edge>(4, vertices[2], vertices[5]));
  edges.push_back(std::make_shared<Edge>(5, vertices[2], vertices[4]));
  edges.push_back(std::make_shared<Edge>(6, vertices[3], vertices[4]));
  edges.push_back(std::make_shared<Edge>(7, vertices[4], vertices[5]));
  edges.push_back(std::make_shared<Edge>(8, vertices[0], vertices[6]));
  edges.push_back(std::make_shared<Edge>(9, vertices[2], vertices[6]));
  edges.push_back(std::make_shared<Edge>(10, vertices[5], vertices[6]));
  edges.push_back(std::make_shared<Edge>(11, vertices[1], vertices[2]));
  edges.push_back(std::make_shared<Edge>(12, vertices[4], vertices[7]));
  edges.push_back(std::make_shared<Edge>(13, vertices[5], vertices[7]));
  edges.push_back(std::make_shared<Edge>(14, vertices[1], vertices[8]));
  edges.push_back(std::make_shared<Edge>(15, vertices[3], vertices[8]));
  
  Graph g(vertices, edges);
  
  Rcout << "Graph vcount = " << g.vcount() << "\n";
  Rcout << "Graph ecount = " << g.ecount() << "\n";
  
  /* Cluster 1: vertices[0], vertices[2], vertices[6]
   * Cluster 2: vertices[1], vertices[3], vertices[8]
   * Cluster 3: vertices[4], vertices[5], vertices[7]
   */
  id_vec cluster {0, 1, 0, 1, 2, 2, 0, 2, 1};
  
  // sample compatible spanning tree
  RNG rn;
  Graph st;
  g.sampleSpanningTree(st, cluster, rn, false);
  
  Rcout << "Tree vcount = " << st.vcount() << "\n";
  Rcout << "Tree ecount = " << st.ecount() << "\n";
  
  // print unrooted spanning tree adjacency list
  Rcout << "Spanning tree vertex adjacency list: \n";
  const std::vector<v_vec> v_adj_list_st = st.getVAdjList();
  for(unsigned int i = 0; i < v_adj_list_st.size(); ++i) {
    Rcout << st.getVertices()[i]->getVid() << " ->" << "\n";
    Rcout << "   ";
    for(auto v : v_adj_list_st[i])
      Rcout << v->getVid() << ", ";
    Rcout << '\n';
  }
  
  // print between edges in spanning tree
  Rcout << "Spanning tree between edges: \n";
  e_set edges_btw = st.getBtwEdges(cluster);
  for(const auto e : edges_btw) {
    Rcout << e->print() << "\n";
  }
  
  // sample a wihtin edge
  Rcout << "Random sampling spanning tree within edges: \n";
  for(int i = 0; i < 10; i++)
    Rcout << (st.sampleWithinEdge(edges_btw, rn))->print() << "\n";
  
  // get a subgraph containning Cluster 1 & 2
  v_unordered_set v_subgraph {vertices[0], vertices[2], vertices[6], vertices[1],
                              vertices[3], vertices[8]};
  Subgraph st_subgraph(st, v_subgraph);
  
  // get connected components with between edge (0, 1) removed
  Rcout << "Clusters after between edge (0 -- 1) is removed in subgraph: \n";
  v_unordered_map<unsigned int> components = st_subgraph.getComponents(edges[0]);
  for(const auto v_component_pair : components)
    Rcout << v_component_pair.first->getVid() << " -> " << v_component_pair.second << "\n";
  
  v_unordered_set component0;
  st_subgraph.getTreeComponents(edges[0], component0);
  Rcout << "Vertices in the component containing vertex " << edges[0]->getEndpoint(0)->getVid() << " :\n";
  for (const auto &v : component0)
    Rcout << v->getVid() << "  ";
  Rcout << "\n";
}

// [[Rcpp::export]]
void benchmarkComponent(IntegerMatrix edge_list) {
  int ne = edge_list.nrow();
  int nv = 0;
  for(int i = 0; i < ne; ++i) {
    nv = std::max(edge_list(i, 0), nv);
    nv = std::max(edge_list(i, 1), nv);
  }
  
  // make graph
  v_vec vertices;
  for(int i = 0; i < nv; ++i) {
    vertices.push_back(std::make_shared<Vertex>(i));
  }
  
  e_vec edges;
  for(int i = 0; i < ne; ++i) {
    unsigned int vid_1 = edge_list(i, 0) - 1;
    unsigned int vid_2 = edge_list(i, 1) - 1;
    edges.push_back(std::make_shared<Edge>(0, vertices[vid_1], vertices[vid_2]));
  }
  
  Graph graph(vertices, edges);
  
  // make subgraph
  v_unordered_set v_subgraph;
  // v_subgraph.reserve(nv);
  for(int i = 0; i < nv; ++i)
    v_subgraph.insert(vertices[i]);
  Subgraph subgraph(graph, v_subgraph);
  
  // get connected component
  subgraph.getComponents(edges[3]);
}

// // [[Rcpp::export]]
// void testHeap() {
//   // pair of (edge_weight, e_ptr)
//   typedef std::pair<double, e_ptr> e_weight_pair;
//
//   // make graph
//   v_vec vertices;
//   for(int i = 0; i < 6; ++i) {
//     vertices.push_back(std::make_shared<Vertex>(i));
//   }
//
//   e_vec edges;
//   edges.push_back(std::make_shared<Edge>(0, vertices[0], vertices[1], 0.1));
//   edges.push_back(std::make_shared<Edge>(1, vertices[0], vertices[2], 0.2));
//   edges.push_back(std::make_shared<Edge>(2, vertices[1], vertices[3], 0.3));
//   edges.push_back(std::make_shared<Edge>(3, vertices[1], vertices[2], 0.5));
//   edges.push_back(std::make_shared<Edge>(4, vertices[2], vertices[3], 0.1));
//
//   boost::heap::fibonacci_heap<e_weight_pair> pq;
//   for(auto e : edges)
//     pq.push(e_weight_pair(e->getWeight(), e));
//
//   Rcout << pq.top().first << "\n";
//   Rcout << pq.top().second->getEid() << "\n";
// }


/*** R
# testVertexEdge()
# testGraph()
# testSubgraph()
# testMST()
# testMSF()
# testHeap()
set.seed(1234)
testRST()
#
# ### benchmark connected component performance
# library(igraph)
# library(microbenchmark)
#
# n = 500
# graph = make_tree(n, 3, "undirected")
# edge_list = as_edgelist(graph)
#
# benchmarkComponent_r <- function(edge_list) {
#   graph = graph_from_edgelist(edge_list, directed = F)
#   subgraph = induced_subgraph(graph, 1:500)
#   subgraph = delete.edges(subgraph, E(subgraph)[3])
#   connect_comp = components(subgraph)
# }
#
# print(microbenchmark(benchmarkComponent(edge_list), benchmarkComponent_r(edge_list),
#                      times = 1000), unit = "relative")
*/
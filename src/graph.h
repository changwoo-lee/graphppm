/**** Undirected graph classes and functions ****/

#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include "graph_data_type.h"
#include "rng.h"


/** General Graph class **/

/*
 IMPORTANT ASSUMPTION:
 Vertex ids are consecutive from 0 to vcount()-1
*/

class Graph {

private:
  // vertices, sorted by vid
  v_vec vertices;

  // edge list, a.k.a, incidence matrix
  e_vec edges;

  // vertex adjacency list
  // where the i-th element stores neighbors of vertex with vid == i
  std::vector<v_vec> v_adj_list;

  // edge adjacency list
  // where the i-th element stores edges adjecent to the vertex with vid == i
  std::vector<e_vec> e_adj_list;

  // function to obtain v_adj_list from edge list
  void computeVAdjList();
  // function to obtain e_adj_list from edge list
  void computeEAdjList();

public:
  Graph() = default;
  Graph(v_vec vertices, e_vec edges, bool sort_v = true): vertices(vertices), edges(edges) {
    // sort vertices by vid
    if (sort_v)
      std::sort(this->vertices.begin(), this->vertices.end(),
                [](v_ptr v_1, v_ptr v_2) {return v_1->getVid() < v_2->getVid();});
  }
  ~Graph() = default;

  // function to get number of vertices
  unsigned int vcount() const {return vertices.size();}
  // function to get number of edges
  unsigned int ecount() const {return edges.size();}

  // function to get vertices
  const v_vec &getVertices() const {return vertices;}
  // function to get edges
  const e_vec &getEdges() const {return edges;}

  // function to get the i-th vertex
  const v_ptr &getVertex(unsigned int i) const {return vertices.at(i);}

  // function to get vertex adjacency list
  const std::vector<v_vec> &getVAdjList() {
    // compute v_adj_list if not existed
    if(v_adj_list.size() == 0)
      computeVAdjList();
    return v_adj_list;
  }
  // function to get edge adjacency list
  const std::vector<e_vec> &getEAdjList() {
    // compute e_adj_list if not existed
    if(e_adj_list.size() == 0)
      computeEAdjList();
    return e_adj_list;
  }

  // function to get between cluster edges
  // @param membership: membership[i] is the cluster membership of vertices[i]
  e_set getBtwEdges(const id_vec &membership) const;

  // function to sample a within-cluster edge
  const e_ptr sampleWithinEdge(const e_set &edges_btw, RNG &rn) const;

  // function to get minimum spanning tree/forest using Prim's algorithm
  // when rooted == true, returned a "rooted" MST whose edge list satisfies:
  // (1) edges[0]->endpoints.first is the pointer to the root
  // (2) edges[i]->endpoints.first is the parent, and edges[i]->endpoints.second is the child
  // @param mst: resulting minimum spanning tree
  // @param max_weight: maximum possible edge weights
  // @param rooted: whether to obtain a rooted version of MST
  // @param update_edges: whether to update edges in mst only (useful for MCMC)
  void getMST(Graph &mst, double max_weight = 1000.0, bool rooted = false, bool update_edges = false);

  // function to sample new spanning tree/forest subject to cluster memberships
  // @param st_graph: resulting minimum spanning tree
  // @param membership: membership[i] is the cluster membership of vertices[i]
  // @param update_edges: whether to update edges spanning tree only (useful for MCMC)
  void sampleSpanningTree(Graph &st_graph, const id_vec &membership, RNG &rn, bool update_edges = false);

};



/** Subgraph class **/
// For simplicity, spanning subtrees are only stored as vertex adjacency list
// Less flexible functionality than the general graph class
// But enough for representing a cluster
class Subgraph {

private:
  // vertices
  v_vec vertices;

  // vertex adjacency map
  // where each element stores neighbors of a vertex
  v_unordered_map<v_set> v_adj_map;

  // number of edges
  unsigned int n_edges;

public:
  Subgraph() = default;
  Subgraph(v_vec vertices, v_unordered_map<v_set> v_adj_map, unsigned int n_edges):
    vertices(vertices), v_adj_map(v_adj_map), n_edges(n_edges) {}
  // construct induced subgraph from a Graph object by vertex ids
  Subgraph(Graph &graph, const v_unordered_set &vertices_subgraph);
  Subgraph(const Subgraph &sg):
    vertices(sg.vertices), v_adj_map(sg.v_adj_map), n_edges(sg.n_edges) {}
  ~Subgraph() = default;

  // function to get number of vertices
  unsigned int vcount() const {return vertices.size();}
  // function to get number of edges
  unsigned int ecount() const {return n_edges;}

  // function to get vertex adjacency map
  const v_unordered_map<v_set> &getVAdjMap() const {return v_adj_map;}

  // function to get vertices
  const v_vec &getVertices() const {return vertices;}
  // function to get edges
  // to-do (O(n) is required)
  e_vec getEdges() const;

  // function to check if a vertex belongs to this subgraph
  bool hasVertex(const v_ptr &v) const {return (v_adj_map.count(v) > 0);}

  // function to get connected components after removing an edge
  v_unordered_map<unsigned int> getComponents(const e_ptr &edge_removed) const;

  // function to connected components after removing an edge for trees
  // not applicable for forests
  // @param component: a vertex pointer set to store the resulting connected component
  //                   containing edge_removed->getEndpoint(0)
  void getTreeComponents(const e_ptr &edge_removed, v_unordered_set &component) const;

  // function to split a subgraph into two
  // by reducing vertices to vertices_old
  // @return: a subgraph with vertices_new
  std::shared_ptr<Subgraph> split(const e_ptr &e_removed,
                                  v_vec &v_vec_old,
                                  v_vec &v_vec_new);

  // function to merge another subgraph sg into this subgraph
  void merge(std::shared_ptr<Subgraph> &sg, const e_ptr &e_merged);

  // function to print all vids
  std::string printVids() const {
    std::string result("Vids: ");
    for (const auto &v : vertices)
      result += std::to_string(v->getVid()) + std::string(" ");
    return result;
  }
};

#endif

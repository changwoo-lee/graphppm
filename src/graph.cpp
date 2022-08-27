/**** Undirected graph classes and functions ****/

#include <stack>
#include <algorithm>
#include <iterator>
#include "graph.h"


/** General graph class **/

// function to obtain v_adj_list from edge list
void Graph::computeVAdjList() {
  v_adj_list.clear();
  v_adj_list.resize(vcount());
  for (unsigned int i = 0; i < ecount(); ++i) {
    v_ptr v_1 = (edges[i]->getEndpoints()).first;
    v_ptr v_2 = (edges[i]->getEndpoints()).second;

    v_adj_list[v_1->getVid()].push_back(v_2);
    v_adj_list[v_2->getVid()].push_back(v_1);
  }
}

// function to obtain e_adj_list from edge list
void Graph::computeEAdjList() {
  e_adj_list.clear();
  e_adj_list.resize(vcount());
  for (unsigned int i = 0; i < ecount(); ++i) {
    e_ptr edge = edges[i];
    v_ptr v_1 = (edge->getEndpoints()).first;
    v_ptr v_2 = (edge->getEndpoints()).second;

    e_adj_list[v_1->getVid()].push_back(edge);
    e_adj_list[v_2->getVid()].push_back(edge);
  }
}

// function to get between cluster edges
// @param membership: membership[i] is the cluster membership of vertices[i]
e_set Graph::getBtwEdges(const id_vec &membership) const {
  e_set edges_btw;
  // edges_btw.reserve(ecount() / 10);
  for (const auto &e : edges) {
    std::pair<v_ptr, v_ptr> endpoints = e->getEndpoints();
    unsigned int vid_1 = (endpoints.first)->getVid();
    unsigned int vid_2 = (endpoints.second)->getVid();
    if (membership[vid_1] != membership[vid_2])
      edges_btw.insert(e);
  }
  return edges_btw;
}

// function to sample a within-cluster edge
// using rejection sampler
const e_ptr Graph::sampleWithinEdge(const e_set &edges_btw, RNG &rn) const {
  e_ptr e;
  unsigned int edge_index;
  do {
    edge_index = rn.rdunif(0, ecount() - 1);
    e = edges[edge_index];
    // Rcpp::Rcout <<  edge_index << " " << e->print() << "\n";
  } while (edges_btw.count(e) > 0);
  return e;
}


// function to get minimum spanning tree/forest using Prim's algorithm
// when rooted == true, returned a "rooted" MST whose edge list satisfies:
// (1) edges[0]->endpoints.first is the pointer to the root
// (2) edges[i]->endpoints.first is the parent, and edges[i]->endpoints.second is the child
// A more efficient implementation can be based on Fibonacci heap in the boost library
// we will see if we need this
// @param mst: resulting minimum spanning tree
// @param max_weight: maximum possible edge weights
// @param rooted: whether to obtain a rooted version of MST
// @param update_edges: whether to update edges in mst only (useful for MCMC)
void Graph::getMST(Graph &mst, double max_weight, bool rooted, bool update_edges) {
  // check if e_adj_list precomputed
  if (e_adj_list.size() == 0)
    this->computeEAdjList();

  unsigned int nv = vcount();
  e_vec &edges_mst = mst.edges;
  if (!edges_mst.empty())  edges_mst.clear();
  edges_mst.reserve(nv - 1);

  // vector to record if a vertex is included in MST
  std::vector<bool> included(nv, false);
  // priority queue to store candidate edges to add to MST
  EdgeHeap pq(nv);

  // build MST
  for (const auto &v : vertices) {
    if (included[v->getVid()])
      continue;

    // use v as a root
    included[v->getVid()] = true;
    for (const auto &e : e_adj_list[v->getVid()]) {
      v_ptr v_neighbor = e->getNeighbor(v);
      if (!included[v_neighbor->getVid()])
        pq.push(v_neighbor, e);
    }

    while (!pq.empty()) {
      e_ptr e = pq.top();
      pq.pop();

      // check if e is in MST
      std::pair<v_ptr, v_ptr> endpoints = e->getEndpoints();
      unsigned int vid_1 = (endpoints.first)->getVid();
      unsigned int vid_2 = (endpoints.second)->getVid();
      if(included[vid_1] && included[vid_2])
        continue;

      v_ptr u; // vertex to be added to MST
      if (included[vid_1])
        u = endpoints.second;
      else
        u = endpoints.first;
      unsigned int vid_u = u->getVid();
      included[vid_u] = true;

      // add e to MST
      if(rooted) {
        // to get rooted tree, rearrange vertex pair to pair(parent, child)
        e_ptr e_directed = std::make_shared<Edge>(e->getEid(), e->getNeighbor(u), u);
        edges_mst.push_back(e_directed);
      } else {
        edges_mst.push_back(e);
      }

      // add edges adjacent to u
      for (const auto &e_u : e_adj_list[vid_u]) {
        double weight = e_u->getWeight();
        v_ptr u_neighbor = e_u->getNeighbor(u);
        unsigned int vid_u_neighbor = u_neighbor->getVid();
        if (!included[vid_u_neighbor]) {
          if (pq.contains(u_neighbor)) {
            if (weight < pq.getWeight(u_neighbor))
              pq.update(u_neighbor, e_u);
          }
          else {
            pq.push(u_neighbor, e_u);
          }
        }
      }
    }
  }

  if (!update_edges)
    mst.vertices = this->vertices;
}

// void Graph::getMST(Graph &mst, double max_weight, bool rooted, bool update_edges) {
//   // check if e_adj_list precomputed
//   if (e_adj_list.size() == 0)
//     this->computeEAdjList();
//
//   unsigned int nv = vcount();
//   e_vec &edges_mst = mst.edges;
//   if (!edges_mst.empty())  edges_mst.clear();
//   edges_mst.reserve(nv - 1);
//
//   // vector to record if a vertex is included in MST
//   std::vector<bool> included(nv, false);
//   // vrctor to record to minimum cost to include a vertex in MST
//   std::vector<double> min_cost(nv, max_weight + 1.0);
//   // priority queue to store candidate edges to add to MST
//   e_priority_queue pq;
//
//   // build MST
//   for (const auto &v : vertices) {
//     if (included[v->getVid()])
//       continue;
//
//     // use v as a root
//     included[v->getVid()] = true;
//     for (const auto &e : e_adj_list[v->getVid()]) {
//       v_ptr v_neighbor = e->getNeighbor(v);
//       if (!included[v_neighbor->getVid()])
//         pq.push(e);
//     }
//
//     while (!pq.empty()) {
//       e_ptr e = pq.top();
//       pq.pop();
//
//       // check if e is in MST
//       std::pair<v_ptr, v_ptr> endpoints = e->getEndpoints();
//       unsigned int vid_1 = (endpoints.first)->getVid();
//       unsigned int vid_2 = (endpoints.second)->getVid();
//       if(included[vid_1] && included[vid_2])
//         continue;
//
//       v_ptr u; // vertex to be added to MST
//       if (included[vid_1])
//         u = endpoints.second;
//       else
//         u = endpoints.first;
//       unsigned int vid_u = u->getVid();
//       included[vid_u] = true;
//
//       // add e to MST
//       if(rooted) {
//         // to get rooted tree, rearrange vertex pair to pair(parent, child)
//         e_ptr e_directed = std::make_shared<Edge>(e->getEid(), e->getNeighbor(u), u);
//         edges_mst.push_back(e_directed);
//       } else {
//         edges_mst.push_back(e);
//       }
//
//       // add edges adjacent to u
//       for (const auto &e_u : e_adj_list[vid_u]) {
//         double weight = e_u->getWeight();
//         unsigned int vid_u_neighbor = (e_u->getNeighbor(u))->getVid();
//         if (!included[vid_u_neighbor] && weight < min_cost[vid_u_neighbor]) {
//           min_cost[vid_u_neighbor] = weight;
//           pq.push(e_u);
//         }
//       }
//     }
//   }
//
//   if (!update_edges)
//     mst.vertices = this->vertices;
// }

// function to sample new spanning tree/forest subject to cluster memberships
// @param st_graph: resulting minimum spanning tree
// @param membership: membership[i] is the cluster membership of vertices[i]
// @param update_edges: whether to update edges spanning tree only (useful for MCMC)
void Graph::sampleSpanningTree(Graph &st_graph, const id_vec &membership, RNG &rn, bool update_edges) {
  // draw edge weights
  for (auto &e : edges) {
    std::pair<v_ptr, v_ptr> endpoints = e->getEndpoints();
    unsigned int vid_1 = (endpoints.first)->getVid();
    unsigned int vid_2 = (endpoints.second)->getVid();
    if (membership[vid_1] != membership[vid_2]) {
      // between cluster edge
      e->setWeight(rn.runif(0.5, 1.0));
    } else {
      // within cluster edge
      e->setWeight(rn.runif(0.0, 0.5));
    }
  }
  this->getMST(st_graph, 1000.0, false, update_edges);
  st_graph.computeVAdjList();
}


/** Subgraph class **/

// construct induced subgraph from a Graph object by vertex ids
Subgraph::Subgraph(Graph &graph, const v_unordered_set &vertices_subgraph) {
  v_vec vertices = graph.getVertices();
  std::vector<v_vec> v_adj_list = graph.getVAdjList();
  v_unordered_map<v_set> v_adj_map_subgraph;
  unsigned int ne_subgraph = 0;

  this->vertices.reserve(vertices_subgraph.size());
  // v_adj_map_subgraph.reserve(2 * vertices_subgraph.size());

  // O(E) version
  for (const auto &v : vertices_subgraph) {
    this->vertices.push_back(v);
    v_adj_map_subgraph[v] = v_set();

    for (const auto &v_neighbor : v_adj_list[v->getVid()]) {
      if ( vertices_subgraph.count(v_neighbor) == 1 ) {
        v_adj_map_subgraph[v].insert(v_neighbor);
        ne_subgraph += 1;
      }
    }
  }

  // // binary search version
  // for (const auto &v : vertices_subgraph) {
  //   this->vertices.push_back(v);
  //   v_adj_map_subgraph[v] = v_set();
  //
  //   // search for v in vertices
  //   std::size_t index;
  //   auto it = std::lower_bound(vertices.begin(), vertices.end(), v,
  //                              [](v_ptr v_1, v_ptr v_2) {return v_1->getVid() < v_2->getVid();});
  //   if (it == vertices.end() || *it != v) {
  //     throw std::invalid_argument("At least one vertex not in the original graph");
  //   } else {
  //     index = std::distance(vertices.begin(), it);
  //   }
  //
  //   for (const auto &v_neighbor : v_adj_list[index]) {
  //     if ( vertices_subgraph.count(v_neighbor) == 1 ) {
  //       v_adj_map_subgraph[v].insert(v_neighbor);
  //       ne_subgraph += 1;
  //     }
  //   }
  // }

  this->v_adj_map = v_adj_map_subgraph;
  this->n_edges = static_cast<unsigned int>(ne_subgraph / 2);
}

// function to get connected components after removing an edge
v_unordered_map<unsigned int> Subgraph::getComponents(const e_ptr &edge_removed) const {
  // get endpoints of the removed edge
  v_ptr v_1 = edge_removed->getEndpoints().first;
  v_ptr v_2 = edge_removed->getEndpoints().second;

  unsigned int nv = vcount();

  // connected component memberships
  v_unordered_map<unsigned int> components;
  // components.reserve(nv);

  // find connected components
  unsigned int component_cnt = 0;
  std::stack<v_ptr> v_stack;
  for (unsigned int i = 0; i < nv; ++i) {
    // check if visited
    if ( components.count(vertices[i]) > 0 ) continue;

    // DFS for the current component
    v_stack.push(vertices[i]);
    while (!v_stack.empty()) {
      v_ptr v = v_stack.top();
      v_stack.pop();
      // check if visited
      if ( components.count(v) == 0 ) {
        // visit v
        components[v] = component_cnt;
        // Rcpp::Rcout << v->getVid() << " ";
        for (const auto &v_neighbor : v_adj_map.at(v)) {
          // check if edge (v, v_neighbor) is removed
          if ( ((v == v_1) && (v_neighbor == v_2)) || ((v == v_2) && (v_neighbor == v_1)) )
            continue;
          v_stack.push(v_neighbor);
        }
      }
    }

    // go for next component
    component_cnt += 1;
  }

  return components;
}

// function to connected components after removing an edge for trees
// not applicable for forests
// @param component: a vertex pointer set to store the resulting connected component
//                   containing edge_removed->getEndpoint(0)
void Subgraph::getTreeComponents(const e_ptr &edge_removed, v_unordered_set &component) const {
  // get endpoints of the removed edge
  v_ptr v_1 = edge_removed->getEndpoint(0);
  v_ptr v_2 = edge_removed->getEndpoint(1);

  if (!hasVertex(v_1))
    throw std::invalid_argument("edge_removed->getEndpoint(0) not in this subgraph");

  unsigned int nv = vcount();
  if (!component.empty())  component.clear();
  component.reserve(nv);

  // find connected components
  std::stack<v_ptr> v_stack;
  v_stack.push(v_1);

  // DFS for the first component containing v_1
  while (!v_stack.empty()) {
    v_ptr v = v_stack.top();
    v_stack.pop();
    // check if visited
    if ( component.count(v) == 0 ) {
      // visit v
      component.insert(v);
      for (const auto &v_neighbor : v_adj_map.at(v)) {
        // check if edge (v, v_neighbor) is removed
        if ( ((v == v_1) && (v_neighbor == v_2)) || ((v == v_2) && (v_neighbor == v_1)) )
          continue;
        v_stack.push(v_neighbor);
      }
    }
  }
}

// function to split a subgraph into two
// by reducing vertices to vertices_old
// @return: a subgraph with vertices_new
std::shared_ptr<Subgraph> Subgraph::split(const e_ptr &e_removed,
                                          v_vec &v_vec_old,
                                          v_vec &v_vec_new) {
  v_ptr v_1 = e_removed->getEndpoint(0), v_2 = e_removed->getEndpoint(1);

  // remove e_removed from vertex adjacency map
  v_adj_map[v_1].erase(v_2);
  v_adj_map[v_2].erase(v_1);

  // remove v_vec_new from vertices
  vertices = v_vec_old;

  // build new subgraph
  v_vec vertices_new = v_vec_new;
  unsigned int ne_new = 0;
  v_unordered_map<v_set> v_adj_map_new;
  for (const auto &v : v_vec_new) {
    v_adj_map_new[v] = v_adj_map[v];
    v_adj_map.erase(v);
    ne_new += v_adj_map_new[v].size();
  }

  // update number of edges
  ne_new /= 2;
  this->n_edges -= ne_new - 1;

  return std::make_shared<Subgraph>(std::move(vertices_new), std::move(v_adj_map_new), ne_new);
}

// function to merge another subgraph sg into this subgraph
void Subgraph::merge(std::shared_ptr<Subgraph> &sg, const e_ptr &e_merged) {
  v_ptr v_1 = e_merged->getEndpoint(0), v_2 = e_merged->getEndpoint(1);

  unsigned int nv_new = this->vcount() + sg->vcount();
  this->vertices.reserve(nv_new);
  // this->v_adj_map.reserve(2 * nv_new);

  // update number of edges
  this->n_edges += sg->ecount() + 1;

  // merge vertices
  v_vec vertices_sg = sg->getVertices();
  this->vertices.insert(this->vertices.end(),
                        std::make_move_iterator(vertices_sg.begin()),
                        std::make_move_iterator(vertices_sg.end()));

  // merge vertex adjacency maps
  v_unordered_map<v_set> v_adj_map_sg = sg->getVAdjMap();
  this->v_adj_map.insert(std::make_move_iterator(v_adj_map_sg.begin()),
                         std::make_move_iterator(v_adj_map_sg.end()));

  // add e_merged to vertex adjacency maps
  this->v_adj_map[v_1].insert(v_2);
  this->v_adj_map[v_2].insert(v_1);
}

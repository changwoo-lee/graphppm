/**** Base Random Spanning Tree Partition Class ****/

#ifndef PARTITION_H
#define PARTITION_H

#include "graph.h"

// subgraph pointer
using subgraph_ptr = std::shared_ptr<Subgraph>;
// double vector
using double_vec = std::vector<double>;

/** Base partition class **/
class Partition {

protected:
  // subgraphs of spanning tree as clusters
  std::vector<subgraph_ptr> subgraphs;

  // membership[i]: cluster membership of vertices[i]
  id_vec membership;

  // spanning tree
  Graph spanning_tree;

  // between cluster edges in spanning tree
  e_set edges_btw_st;

  // number of clusters
  unsigned int n_clusters;

  // flag to indicate if the number of clusters has been changed
  bool n_clusters_changed = true;

  /* temporary variables for split proposals */

  // when a cluster c is split into c1 and c2, where c2 is a new cluster (of index k+1)
  // vertices_old_ stores vertices for c1 (AFTER splitting), vertices_new_ stores vertices for c2
  v_vec vertices_old_split_;
  v_vec vertices_new_split_;

  // cid_split_ stores cluster id of c
  unsigned int cid_split_;

  // _edge_removed stores the edge to be removed
  e_ptr e_removed_;


  /* temporary variables for merge proposals */

  // when cluster c2 is merged into c1
  // vertices_old_ stores vertices for c2, vertices_new_ stores in the merged cluster
  // v_unordered_set vertices_old_merge_;
  // v_unordered_set vertices_new_merge_;

  // cid_old_merge_ stores cluster id of c2 (BEFORE merging)
  unsigned int cid_old_merge_;
  // cid_new_merge_ stores cluster id of c1 (BEFORE merging)
  unsigned int cid_new_merge_;

  // _edge_merged stores the edge for merging
  e_ptr e_merged_;

  /* friend classes */

  friend class PartitionsBSCC;

public:
  Partition(std::vector<subgraph_ptr> subgraphs, id_vec membership,
            Graph spanning_tree, e_set edges_btw_st):
    subgraphs(subgraphs), membership(membership), spanning_tree(spanning_tree),
    edges_btw_st(edges_btw_st), n_clusters(subgraphs.size()) {}

  // @param membership has to be consecutive starting from 0
  Partition(const id_vec &membership, const Graph &spanning_tree);

  ~Partition() = default;

  // function to get subgraphs
  const std::vector<subgraph_ptr> &getSubgraphs() const {return subgraphs;}

  // function to get the i-th subgraphs
  const subgraph_ptr &getSubgraph(unsigned int i) const {return subgraphs.at(i);}

  // function to get membership
  const id_vec &getMembership() const {return membership;}

  // function to get membership of vertices[i]
  unsigned int getMembership(unsigned int i) const {return membership.at(i);}

  // function to get spanning tree
  const Graph &getSpanningTree() const {return spanning_tree;}

  // function to get between cluster edges in spanning tree
  const e_set &getBtwTreeEdges() const {return edges_btw_st;}

  // function to get number of clusters
  unsigned int ccount() const {return n_clusters;}

  // function to get size of the i-th cluster
  unsigned int csize(unsigned int i) const {
    if (i >= n_clusters)
      throw std::invalid_argument("Index out of bound");
    return subgraphs[i]->vcount();
  }

  // function to get vertices in the i-th cluster
  const v_vec &getVertices(unsigned int i) const {
    if (i >= n_clusters)
      throw std::invalid_argument("Index out of bound");
    return subgraphs[i]->getVertices();
  }

  // function to split a cluster
  // @param change: whether this is in a change move
  void splitCluster(RNG &rn, bool change = false);

  // function to merge two clusters
  void mergeCluster(RNG &rn);

  // function to update partition after split
  // @param change: whether this is in a change move
  void updateSplit(bool change = false);

  // function to update partition after merge
  // @param change: whether this is in a change move
  void updateMerge(bool change = false);

  // function to sample spanning tree
  void sampleSpanningTree(Graph &graph, RNG &rn);

  // functions to print split/merge results (for debug/testing purposes)
  std::string printSplit() const;
  std::string printMerge() const;
};

#endif

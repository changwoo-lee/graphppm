/**** Base Random Spanning Tree Partition Class ****/

#include <cmath>
#include <algorithm>
#include "partition.h"


/** Base partition class **/

// constructor
// @param membership has to be consecutive starting from 0
Partition::Partition(const id_vec &membership, const Graph &spanning_tree) {
  this->membership = membership;
  this->spanning_tree = spanning_tree;

  const v_vec &vertices = spanning_tree.getVertices();
  n_clusters = *std::max_element(membership.begin(), membership.end()) + 1;
  std::vector<v_unordered_set> v_set_vec(n_clusters);
  for (unsigned int i = 0; i < membership.size(); ++i)
    v_set_vec[membership.at(i)].insert( vertices.at(i) );

  // build subgraphs
  subgraphs.clear();
  for (unsigned int i = 0; i < n_clusters; ++i)
    subgraphs.push_back( std::make_shared<Subgraph>(this->spanning_tree, v_set_vec[i]) );

  // get between-cluster edges
  edges_btw_st.clear();
  for (const auto &e : spanning_tree.getEdges()) {
    unsigned int vid_1 = (e->getEndpoint(0))->getVid();
    unsigned int vid_2 = (e->getEndpoint(1))->getVid();
    if (membership[vid_1] != membership[vid_2])
      edges_btw_st.insert(e);
  }
}

// function to split a cluster
void Partition::splitCluster(RNG &rn, bool change) {
  // reset the flag to indicate if the number of clusters has been changed
  n_clusters_changed = false;

  // sample a within edge to remove
  if (!change) {
    e_removed_ = spanning_tree.sampleWithinEdge(edges_btw_st, rn);
  } else {
    // this is in a change step
    e_set edges_btw_st_new(edges_btw_st);
    edges_btw_st_new.erase(e_merged_);
    e_removed_ = spanning_tree.sampleWithinEdge(edges_btw_st_new, rn);
  }

  // get the cluster to split
  v_ptr endpoint = e_removed_->getEndpoint(0);
  cid_split_ = membership[ endpoint->getVid() ];
  // Rcpp::Rcout<< "cid_split_ = " << cid_split_ << "\n";

  v_unordered_set component_0;
  vertices_old_split_.clear();
  vertices_new_split_.clear();

  // get connected components after removing e_removed
  subgraphs[cid_split_]->getTreeComponents(e_removed_, component_0);
  // cluster sizes
  unsigned int csize_1 = component_0.size();
  unsigned int csize_2 = subgraphs[cid_split_]->vcount() - csize_1;

  if ((!change) ||
      (change && cid_split_ != cid_new_merge_ && cid_split_ != cid_old_merge_)) {
    /* Two cases:
     * 1) This is not in a change move; OR
     * 2) This is in a change move and we're not splitting the merged cluster
     */

    // set the smaller cluster containing endpoint as new
    if (csize_1 < csize_2) {
      // component_0 as new
      for (const auto &v : subgraphs[cid_split_]->getVertices()) {
        if (component_0.count(v) == 0)
          vertices_old_split_.push_back(v);
        else
          vertices_new_split_.push_back(v);
      }
    } else {
      // component_0 as old
      for (const auto &v : subgraphs[cid_split_]->getVertices()) {
        if (component_0.count(v) == 0)
          vertices_new_split_.push_back(v);
        else
          vertices_old_split_.push_back(v);
      }
    }

  } else {
    /* this is in a change step and we're splitting the merged cluster */

    // cluster id of the merged one not being split
    unsigned int cid_no_split = (cid_split_ == cid_new_merge_) ? cid_old_merge_ : cid_new_merge_;
    // cluster membership of one of the endpoints of e_merged_
    unsigned int membership_endpoint_merge = membership[ (e_merged_->getEndpoint(0))->getVid() ];
    // endpoint of e_merged_ that is in the cluster being split
    v_ptr endpoint_merge_split = (membership_endpoint_merge == cid_split_) ?
                                  e_merged_->getEndpoint(0) :
                                  e_merged_->getEndpoint(1);
    // check if cid_no_split is not in component_0
    bool cid_no_split_in_comp_0 = (e_merged_ != e_removed_) && (component_0.count(endpoint_merge_split) > 0);

    if (cid_no_split_in_comp_0)
      // cid_no_split is in component_0
      csize_1 += subgraphs[cid_no_split]->vcount();
    else
      // cid_no_split is not in component_0
      csize_2 += subgraphs[cid_no_split]->vcount();

    // set the smaller cluster containing endpoint as new
    if (csize_1 < csize_2) {
      // component_0 as new
      for (const auto &v : subgraphs[cid_split_]->getVertices()) {
        if (component_0.count(v) == 0)
          vertices_old_split_.push_back(v);
        else
          vertices_new_split_.push_back(v);
      }
      // add cid_no_split
      if (cid_no_split_in_comp_0) {
        // cid_no_split is in component_0
        for (const auto &v : subgraphs[cid_no_split]->getVertices())
          vertices_new_split_.push_back(v);
      } else {
        // cid_no_split is not in component_0
        for (const auto &v : subgraphs[cid_no_split]->getVertices())
          vertices_old_split_.push_back(v);
      }
    } else {
      // component_0 as old
      for (const auto &v : subgraphs[cid_split_]->getVertices()) {
        if (component_0.count(v) == 0)
          vertices_new_split_.push_back(v);
        else
          vertices_old_split_.push_back(v);
      }
      // add cid_no_split
      if (cid_no_split_in_comp_0) {
        // cid_no_split is in component_0
        for (const auto &v : subgraphs[cid_no_split]->getVertices())
          vertices_old_split_.push_back(v);
      } else {
        // cid_no_split is not in component_0
        for (const auto &v : subgraphs[cid_no_split]->getVertices())
          vertices_new_split_.push_back(v);
      }
    }
  }

  if (change) {
    // set cid_split_ to be the cluster id AFTER merging
    if (cid_split_ == cid_old_merge_)
      cid_split_ = cid_new_merge_;
    if (cid_split_ > cid_old_merge_)
      cid_split_ -= 1;
  }
}

// function to merge two clusters
void Partition::mergeCluster(RNG &rn) {
  // reset the flag to indicate if the number of clusters has been changed
  n_clusters_changed = false;

  // sample a between cluster edge
  e_vec edges_btw_st_vec;
  edges_btw_st_vec.reserve(n_clusters - 1);
  for (const auto &e_btw : edges_btw_st)
    edges_btw_st_vec.push_back(e_btw);
  e_merged_ = edges_btw_st_vec[ rn.rdunif(0, n_clusters - 2) ];
  // Rcpp::Rcout <<  e_merged_->print() << "\n";

  // get clusters to merge; merge smaller one into larger one
  v_ptr endpoint_1 = e_merged_->getEndpoint(0);
  v_ptr endpoint_2 = e_merged_->getEndpoint(1);
  unsigned int membership_1 = membership[ endpoint_1->getVid() ];
  unsigned int membership_2 = membership[ endpoint_2->getVid() ];
  unsigned int csize_1 = subgraphs[membership_1]->vcount();
  unsigned int csize_2 = subgraphs[membership_2]->vcount();
  if (csize_1 < csize_2) {
    cid_new_merge_ = membership_2;
    cid_old_merge_ = membership_1;
  } else {
    cid_new_merge_ = membership_1;
    cid_old_merge_ = membership_2;
  }

  // // cluster sizes
  // unsigned int csize_1 = subgraphs[cid_new_merge_]->vcount();
  // unsigned int csize_2 = subgraphs[cid_old_merge_]->vcount();
  //
  // merge c2 to c1
  // vertices_old_merge_.clear();
  // vertices_new_merge_.clear();
  // vertices_old_merge_.reserve(csize_2);
  // vertices_new_merge_.reserve(csize_1 + csize_2);
  // for (const auto &v : subgraphs[cid_new_merge_]->getVertices())
  //   vertices_new_merge_.insert(v);
  // for (const auto &v : subgraphs[cid_old_merge_]->getVertices()) {
  //   vertices_old_merge_.insert(v);
  //   vertices_new_merge_.insert(v);
  // }
}



// function to update partition after split
// @param change: whether this is in a change move
void Partition::updateSplit(bool change) {
  // update subgraphs
  subgraphs.push_back( subgraphs[cid_split_]->split(e_removed_, vertices_old_split_, vertices_new_split_) );
  // update number of clusters
  n_clusters += 1;
  // update cluster membership
  for (const auto &v : vertices_new_split_)
    membership[v->getVid()] = n_clusters - 1;
  // update between cluster edge
  edges_btw_st.insert(e_removed_);

  if (!change) n_clusters_changed = true;
}

// function to update partition after merge
// @param change: whether this is in a change move
void Partition::updateMerge(bool change) {
  // update subgraphs
  subgraphs[cid_new_merge_]->merge(subgraphs[cid_old_merge_], e_merged_);
  subgraphs[cid_old_merge_].reset();
  subgraphs.erase(subgraphs.begin() + cid_old_merge_);
  // update number of clusters
  n_clusters -= 1;
  // update cluster membership
  for (unsigned int i = 0; i < membership.size(); ++i) {
    if (membership[i] == cid_old_merge_)
      membership[i] = cid_new_merge_;
    if (membership[i] > cid_old_merge_)
      membership[i] -= 1;
  }
  // update between cluster edge
  edges_btw_st.erase(e_merged_);

  if (!change) n_clusters_changed = true;
}

// function to sample spanning tree
void Partition::sampleSpanningTree(Graph &graph, RNG &rn) {
  // sample a new spanning tree
  graph.sampleSpanningTree(spanning_tree, membership, rn, true);

  // update between cluster edges in spanning tree
  edges_btw_st.clear();
  for (const auto &e : spanning_tree.getEdges()) {
    unsigned int vid_1 = (e->getEndpoint(0))->getVid();
    unsigned int vid_2 = (e->getEndpoint(1))->getVid();
    if (membership[vid_1] != membership[vid_2])
      edges_btw_st.insert(e);
  }

  // update subgraphs
  v_unordered_set vertices_subgraph;
  // vertices_subgraph.reserve(spanning_tree.vcount());
  for (unsigned int k = 0; k < n_clusters; ++k) {
    vertices_subgraph.clear();
    for (const auto &v : subgraphs[k]->getVertices())
      vertices_subgraph.insert(v);
    subgraphs[k] = std::make_shared<Subgraph>(spanning_tree, vertices_subgraph);
  }

  n_clusters_changed = false;
}

// functions to print split/merge results (for debug/testing purposes)
std::string Partition::printSplit() const {
  std::string result = std::string("removed eid = ") + e_removed_->print() +
    std::string(" in Cluster ") + std::to_string(cid_split_) + std::string("\n");

  result += std::string("vids in c1: \n");
  for (const auto &v : vertices_old_split_)
    result += std::to_string(v->getVid()) + std::string(" ");

  result += std::string("\n") + std::string("vids in c2: \n");
  for (const auto &v : vertices_new_split_)
    result += std::to_string(v->getVid()) + std::string(" ");

  return result;
}

std::string Partition::printMerge() const {
  std::string result = std::string("merge eid = ") + std::to_string(e_merged_->getEid()) + std::string("\n");
  result += std::string("merging Cluster ") + std::to_string(cid_old_merge_) +
    std::string(" into Cluster ") + std::to_string(cid_new_merge_) + std::string("\n");

  result += std::string("vids in c2: \n");
  for (const auto &v : subgraphs[cid_old_merge_]->getVertices())
    result += std::to_string(v->getVid()) + std::string(" ");

  result += std::string("\n") + std::string("vids in c1 after merging: \n");
  for (const auto &v : subgraphs[cid_new_merge_]->getVertices())
    result += std::to_string(v->getVid()) + std::string(" ");
  for (const auto &v : subgraphs[cid_old_merge_]->getVertices())
    result += std::to_string(v->getVid()) + std::string(" ");

  return result;
}

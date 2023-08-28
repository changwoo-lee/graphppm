/**** Graph Related Data Types ****/

#ifndef GRAPH_DATA_TYPE_H
#define GRAPH_DATA_TYPE_H

#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include "vertex.h"
#include "edge.h"
#include <stdexcept>
// edge pointer
using e_ptr = std::shared_ptr<Edge>;
// vector of vertices
using v_vec = std::vector<v_ptr>;
// vector of edges
using e_vec = std::vector<e_ptr>;
// vector of vertex/edge/cluster ids
using id_vec = std::vector<unsigned int>;
// set of vertex/edge ids
using id_set = std::unordered_set<unsigned int>;

/** Ordered Sets for Vertex/Edge Pointers **/

// function to compare two edge pointers
struct CompareEPtr {
  bool operator() (const e_ptr &e_1, const e_ptr &e_2) const {
    return (e_1->getEid() > e_2->getEid());
  }
};

// function to compare two vertex pointers
struct CompareVPtr {
  bool operator() (const v_ptr &v_1, const v_ptr &v_2) const {
    return (v_1->getVid() > v_2->getVid());
  }
};

// set of edge pointers
using e_set = std::set<e_ptr, CompareEPtr>;
// set of vertex pointers
using v_set = std::set<v_ptr, CompareVPtr>;


/** Hash Maps Using Vertex/Edge Pointers as Keys **/

// hash function for vertex pointers
// in case the default one is giving good enough performance
struct VertexHashFunc {
  std::size_t operator()(const v_ptr &v) const {
    return static_cast<std::size_t>( v->getVid() );
  }
};

// hash function for edge pointers
struct EdgeHashFunc {
  std::size_t operator()(const e_ptr &e) const {
    return static_cast<std::size_t>( e->getEid() );
  }
};

// hash map using vertex pointers as keys
template <class T>
using v_unordered_map = std::unordered_map<v_ptr, T, VertexHashFunc>;

// unordered set of vertex pointers
using v_unordered_set = std::unordered_set<v_ptr, VertexHashFunc>;

// hash map using edge pointers as keys
template <class T>
using e_unordered_map = std::unordered_map<e_ptr, T, EdgeHashFunc>;

// unordered set of edge pointers
using e_unordered_set = std::unordered_set<e_ptr, EdgeHashFunc>;


/** Priority queue/min heap for weighed edge pointers **/

// function to compare two weighted edge pointers by weights
struct CompareEWeights {
  bool operator() (const e_ptr &e_1, const e_ptr &e_2) const {
    return (e_1->getWeight() > e_2->getWeight());
  }
};

using e_priority_queue = std::priority_queue< e_ptr, e_vec, CompareEWeights >;

/*
 * For the propose of Prim's algorithm, the heap stores e_ptr and
 * is indexed by vid, where vid is the id of an endpoint
 * of the edge e_ptr points to, which is the vertex to be added to MST
 * By design, the maximum heap size is the number of vertices in the graph,
 * and we can update the edge reaching vid in O(|V|) time
 */

class EdgeHeap {

private:
  // edge vector where edges[i] is the edge reaching vertex i
  e_vec edges;

  // min-heap of edge indices: edge[ heap[0] ] is the edge with minimum weight
  id_vec heap;

  // heap indices: heap[ heap_idx[i] ] is the index in edges reaching vertex i
  // note: heap_idx[heap[i]] == i
  id_vec heap_idx;

  // // size of heap
  // unsigned int heap_size = 0;

  // heap capacity, typically |V|
  unsigned int capacity = 0;

  // function to build a heap rooted at heap[i]
  // assuming that the binary trees rooted at heap[left(i)] and heap[right(i)] are heaps
  void heapify(unsigned int i);

  /* Useful functions */

  // function to obtain index of parent
  inline static unsigned int parent(unsigned int i) {return (i - 1) / 2;}

  // function to obtain index of left child
  inline static unsigned int left(unsigned int i) {return 2 * i + 1;}

  // function to obtain index of right child
  inline static unsigned int right(unsigned int i) {return 2 * i + 2;}

  // function to compare weights of edges[heap[i]] and edges[heap[j]]
  inline bool smaller(unsigned int i, unsigned int j) {
    return ( edges[heap[i]]->getWeight() < edges[heap[j]]->getWeight() );
  }

public:
  EdgeHeap() = default;

  EdgeHeap(unsigned int capacity): capacity(capacity) {
    edges.resize(capacity);
    heap.reserve(capacity);
    heap_idx.resize(capacity, capacity + 1);
  }

  EdgeHeap(e_vec edges): edges(edges), capacity(edges.size()) {
    unsigned int n = edges.size();
    heap.reserve(n);
    heap_idx.reserve(n);
    for (unsigned int i = 0; i < n; ++i) {
      heap.push_back(i);
      heap_idx.push_back(i);
    }
    // build heap
    for (unsigned int i = 1; i <= n / 2; ++i)
      heapify(n / 2 - i);
  }

  ~EdgeHeap() = default;

  // function to check if heap is empty
  bool empty() const {return heap.size() == 0;}

  // function to get heap size
  unsigned int size() const {return heap.size();}

  // function to check if the heap contains an edge reaching to a vertex v
  bool contains(const v_ptr &v) const {
    return (v->getVid() < heap_idx.size()) && (heap_idx[v->getVid()] != capacity + 1);
  }

  // function to obtain the edge with minimum edge weight
  const e_ptr &top() {return edges[heap[0]];}

  // function to remove the edge with minimum edge weight
  void pop();

  // function to add a new edge reaching vertex <to>
  void push(const v_ptr &to, const e_ptr &edge);

  // function to change the edge reaching vertex <to> to an edge with smaller weight
  void update(const v_ptr &to, const e_ptr &edge);

  // function to get weight of the edge reaching v
  double getWeight(const v_ptr v) const {
    if (!contains(v))
      throw std::invalid_argument("vertex has no associated edge in heap");
    return edges[v->getVid()]->getWeight();
  }

  // function to print heap
  std::string printHeap() const {
    std::string info("");
    for (unsigned int i = 0; i < heap.size(); ++i)
      info += edges[heap[i]]->print() + ", weight = " +
        std::to_string(edges[heap[i]]->getWeight()) + "\n";
    return info;
  }

  // function to print heap index
  std::string printHeapIndex() const {
    std::string info("");
    for (unsigned int i = 0; i < heap.size(); ++i)
      info += edges[heap[i]]->print() + ", weight = " +
        std::to_string(edges[heap[i]]->getWeight()) +
        ", heap idx = " + std::to_string(heap_idx[heap[i]]) + "\n";
    return info;
  }
};


#endif

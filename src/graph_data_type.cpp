/**** Graph Related Data Types ****/

#include <utility>
#include "graph_data_type.h"

/** Priority queue/min heap for weighed edge pointers **/

// function to build a heap rooted at heap[i]
// assuming that the binary trees rooted at heap[left(i)] and heap[right(i)] are heaps
void EdgeHeap::heapify(unsigned int i) {
  unsigned int l = left(i), r = right(i);
  // index with smallest edge weight among heap[i], heap[l], heap[r]
  unsigned int min_idx = i;

  if (l < heap.size() && smaller(l, i))
    min_idx = l;
  if (r < heap.size() && smaller(r, min_idx))
    min_idx = r;

  if (min_idx != i) {
    // swap heap[i] and heap[min_idx]
    heap_idx[heap[i]] = min_idx;
    heap_idx[heap[min_idx]] = i;
    std::swap(heap[i], heap[min_idx]);
    this->heapify(min_idx);
  }
}

// function to remove the edge with minimum edge weight
void EdgeHeap::pop() {
  if (empty())
    throw std::runtime_error("cannot pop an empty heap");

  unsigned int last_idx = heap.size() - 1;
  heap_idx[heap[0]] = capacity + 1;
  heap_idx[heap[last_idx]] = 0;
  std::swap(heap[0], heap[last_idx]);
  heap.pop_back();

  // maintain heap property
  heapify(0);
}

// function to add a new edge
void EdgeHeap::push(const v_ptr &to, const e_ptr &edge) {
  unsigned int vid = to->getVid();
  if (vid >= capacity)
    throw std::invalid_argument("not enough heap capacity");

  // check if there is already an edge reaching to
  if (!contains(to)) {
    edges[vid] = edge;
    heap_idx[vid] = heap.size();
    heap.push_back(vid);

    // maintain heap property
    unsigned int idx = heap.size() - 1;
    unsigned int parent_idx = parent(idx);
    while (idx > 0 && smaller(idx, parent_idx)) {
      // swap i with its parent
      heap_idx[heap[idx]] = parent_idx;
      heap_idx[heap[parent_idx]] = idx;
      std::swap(heap[idx], heap[parent_idx]);
      idx = parent_idx;
      parent_idx = parent(idx);
    }
  } else {
    throw std::invalid_argument("there is an edge associated to this vertex in heap");
  }
}

// function to change the edge reaching vertex <to> to an edge with smaller weight
void EdgeHeap::update(const v_ptr &to, const e_ptr &edge) {
  if (!contains(to))
    throw std::invalid_argument("vertex has no associated edge in heap");

  unsigned int vid = to->getVid();
  if (edges[vid]->getWeight() >= edge->getWeight()) {
    edges[vid] = edge;
    // maintain heap property
    unsigned int idx = heap_idx[vid];
    unsigned int parent_idx = parent(idx);
    while (idx > 0 && smaller(idx, parent_idx)) {
      // swap i with its parent
      heap_idx[heap[idx]] = parent_idx;
      heap_idx[heap[parent_idx]] = idx;
      std::swap(heap[idx], heap[parent_idx]);
      idx = parent_idx;
      parent_idx = parent(idx);
    }
  } else {
    throw std::invalid_argument(edge->print() + " has larger weight");
  }
}


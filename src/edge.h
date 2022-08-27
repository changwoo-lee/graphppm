/**** Edge Class ****/

#ifndef EDGE_H
#define EDGE_H

#include <utility>
#include <memory>
#include <string>
#include "vertex.h"
#include <stdexcept>
// vertex pointer
using v_ptr = std::shared_ptr<Vertex>;


class Edge {

private:
  // edge id
  unsigned int eid;
  // endpoints
  std::pair<const v_ptr, const v_ptr> endpoints;
  // edge weight
  double weight = 0.0;

public:
  Edge(const unsigned int id, const std::pair<const v_ptr, const v_ptr> &endpoints):
    eid(id), endpoints(endpoints) {}
  Edge(const unsigned int id, const v_ptr endpoint_1, const v_ptr endpoint_2):
    eid(id), endpoints(endpoint_1, endpoint_2) {}
  Edge(const unsigned int id, const std::pair<const v_ptr, const v_ptr> &endpoints, double weight):
    eid(id), endpoints(endpoints), weight(weight) {}
  Edge(const unsigned int id, const v_ptr endpoint_1, const v_ptr endpoint_2, double weight):
    eid(id), endpoints(endpoint_1, endpoint_2), weight(weight) {}
  ~Edge() = default;

  // function to get eid
  const unsigned int &getEid() const {return eid;}

  // function to get endpoints
  const std::pair<const v_ptr, const v_ptr> &getEndpoints() {return endpoints;}

  // function to get one endpoint
  const v_ptr getEndpoint(unsigned int idx) {
    if (idx == 0)
      return endpoints.first;
    else if (idx == 1)
      return endpoints.second;
    else
      throw std::invalid_argument("idx is either 0 or 1");
  }

  // function to get neighbor of a vertex connected by this edge
  const v_ptr getNeighbor(const v_ptr &v) const {
    v_ptr u1 = endpoints.first;
    v_ptr u2 = endpoints.second;
    if(v == u1) return u2;
    if(v == u2) return u1;
    throw std::invalid_argument("vertex is not connected by the edge");
  }

  // function to get edge weight
  double getWeight() const {return weight;}

  // function to set edge weight
  void setWeight(double weight) {this->weight = weight;}

  // function to print edge
  // preferably this should be done by overloading operator <<
  std::string print() const {
    std::string info = std::string("Edge id = ") + std::to_string(eid) + std::string(": ") +
      std::to_string((endpoints.first)->getVid()) + std::string(" --- ") +
      std::to_string((endpoints.second)->getVid());
    return info;
  }
};


#endif

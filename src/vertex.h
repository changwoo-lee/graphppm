/**** Vertex Class ****/

#ifndef VERTEX_H
#define VERTEX_H


class Vertex {

private:
  // vertex id
  unsigned int vid;

public:
  explicit Vertex(const unsigned int id) : vid(id) {}
  ~Vertex() = default;

  // function to get vertex id
  unsigned int getVid() const {return vid;}

  bool operator==(const Vertex &v) const {return this->vid == v.vid;}
  bool operator< (const Vertex &v) const {return this->vid < v.vid;}
};

#endif

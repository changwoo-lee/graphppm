#include <Rcpp.h>
using namespace Rcpp;

#include "vertex.h"
#include "edge.h"
#include "graph_data_type.h"
#include "graph.h"
#include "rng.h"

// [[Rcpp::export]]
Rcpp::IntegerMatrix runcppWilson(int n,
                                 Rcpp::IntegerMatrix graph_edge_list,
                                 int rootvertex,
                                 Rcpp::IntegerVector ordering,
                                 int maxwalk) {
  /* convert R objects to C++ */
  v_vec vertices;
  e_vec edges;
  vertices.reserve(n);
  edges.reserve(graph_edge_list.nrow());
  for (unsigned int i = 0; i < n; ++i)
    vertices.push_back(std::make_shared<Vertex>(i));
  for (int i = 0; i < graph_edge_list.nrow(); ++i) {
    unsigned int vid_1 = graph_edge_list(i, 0) - 1; // change to c++ indexing
    unsigned int vid_2 = graph_edge_list(i, 1) - 1; // change to c++ indexing
    edges.push_back(std::make_shared<Edge>(i, vertices[vid_1], vertices[vid_2]));
  }
  Graph g(vertices, edges);
  
  rootvertex = rootvertex - 1; // change to c++ indexing
  
  /* convert R objects to C++ */
  //Rcout << "vcount = " << g.vcount() << "\n";
  //Rcout << "ecount = " << g.ecount() << "\n";
  
  //Rcout << "Vertex ids: \n";
  //for (const auto v : g.getVertices())
  //  Rcout << v->getVid() << "\n";
  
  const std::vector<v_vec> v_adj_list = g.getVAdjList();
  //Rcout << "v_adj_list size = " << v_adj_list.size() << "\n";
  
  // print adjacency list
  //Rcout << "Vertex adjacency list: \n";
  //for(unsigned int i = 0; i < v_adj_list.size(); ++i) {
  //  Rcout << g.getVertices()[i]->getVid() << " ->" << "\n";
  //  Rcout << "   ";
  //  for(auto v : v_adj_list[i])
  //    Rcout << v->getVid() << ", ";
  //  Rcout << '\n';
  //}
  
  // random number generator
  RNG rn;
  
  int neighborsize;
  int idx;
  int randomwalk;
  int startrow;
  // for loop-erase
  std::vector<int> ivec; 
  int i_j, inext;
  std::vector<int> trace_new;
  
  
  // initialize tree_vertex
  std::vector<int> tree_vertex;
  tree_vertex.push_back(rootvertex);
  //Rcout << "tree_vertex initial:"<< tree_vertex[0] << "\n";
  // initialize tree_edge
  Rcpp::IntegerMatrix tree_edge(n-1, 2);
  
  std::vector<int> trace;
  
  for(int i=0; i<n; i++){
    if(tree_vertex.size()==n) break;
    
    trace = {ordering[i]-1}; // change to c++ ordering; starting from 0
    if(trace[0]==rootvertex) continue;
    //Rcout << "trace initial:"<< trace[0] << "\n";
    
    
    // random walk
    // check break
    for(int iwalk = 0; iwalk < maxwalk; iwalk++ ){
      // check if hits tree_vertex
      if( std::find(tree_vertex.begin(), tree_vertex.end(), trace.back()) != tree_vertex.end() )
      {
        break; //trace.back() is an element of tree_vertex
      }
      //trace.back() : access last element
      neighborsize = (&v_adj_list[trace.back()])->size();
      idx = rn.rdunif(0, neighborsize -1);
      randomwalk = v_adj_list[trace.back()][idx]->getVid();
      trace.push_back(randomwalk);
      //Rcout << randomwalk << ": randomwalk" << "\n";
      if(iwalk == maxwalk - 1) Rcout << "maximum random walk iteration reached!!\n" ;
    }
    //for(auto t : trace) Rcout << t << ",";
    //Rcout << "before loop erase \n";
    
    // loop_erase the trace
    if(trace.size()>2){
      //Rcout << "loop erase start \n";
      ivec = {0};
      trace_new = {trace[0]};
      for(int i = 0; i < trace.size(); i++ ){
        i_j = ivec.back() ;
        //Rcout << "find" << trace[i_j] << "\n";
        //https://stackoverflow.com/questions/24997910/get-index-in-vector-from-reverse-iterator
        auto rit = std::find(trace.rbegin(), trace.rend(), trace[i_j]);
        //if (rit != trace.rend()) {
        inext = std::distance(begin(trace), rit.base()); //- 1 + 1;
        //}else{
        //  Rcout << "error\n";
        //}
        //Rcout << "inext:"<<inext<<"\n";
        ivec.push_back(inext);
        //for(auto t : ivec) Rcout << t << ",";
        //Rcout << "ivec \n";
        trace_new.push_back(trace[inext]);
        if(trace[inext]==trace.back()) break;
      }
    }else{
      trace_new = trace;
    }
    //for(auto t : trace_new) Rcout << t << ", ";
    //Rcout << "after loop erase \n";
    
    if(trace_new.size()>1){
      startrow = tree_vertex.size()-1;
      for(int j = 0; j < trace_new.size()-1 ; j++){
        tree_edge(startrow + j,0) = trace_new[j] + 1; // return back to R indexing
        tree_edge(startrow + j,1) = trace_new[j+1] + 1; // return back to R indexing
        tree_vertex.push_back(trace_new[j]);
      }
    }
    
  }
  
  return tree_edge;
}


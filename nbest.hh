/*** written by Thomas Schoenemann as an employee of the University of Pisa, Italy, October 2011 ***/

#include "vector.hh"
#include <vector>

class DAGEdge;

class DAGNode {
public:

  DAGNode(double forward);

  void add_edge(DAGEdge* edge);

  double forward_potential_; //minimal distance to start node

  Storage1D<DAGEdge*> incoming_edge_;

  uint nVisits_;
};


class DAGEdge {
public:

  DAGEdge(DAGNode* from, DAGNode* to, double weight, const Storage1D<uint>& output);

  DAGNode* from_;
  DAGNode* to_;
  double weight_;
  Storage1D<uint> output_;
};


class DAG {
public:

  DAG(uint nNodes, uint nEdges);

  ~DAG();

  uint add_node(double forward_potential);

  void add_edge(uint from, uint to, double weight, const Storage1D<uint>& output);

  DAGNode* node(uint num);

  //note: this assumes that accurate forward potentials have been set
  void nbest(uint n, uint start_node, uint end_node,
	     std::vector<Storage1D<uint> >& sequence, std::vector<double>* score);

protected:

  Storage1D<DAGNode*> node_;
  Storage1D<DAGEdge*> edge_;

  uint nUsedNodes_;
  uint nUsedEdges_;
};

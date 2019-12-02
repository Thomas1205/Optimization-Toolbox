/*** written by Thomas Schoenemann as an employee of the University of Pisa, Italy, October 2011 ***/

#include "nbest.hh"

#include <queue>

struct NBestState {

  NBestState(DAGNode* node, double cost, NBestState* prev, const Storage1D<uint>& output);

  DAGNode* node_;
  double cost_;

  NBestState* prev_;
  Storage1D<uint> output_;
};

NBestState::NBestState(DAGNode* node, double cost, NBestState* prev, const Storage1D<uint>& output) :
  node_(node), cost_(cost), prev_(prev), output_(output) {}

struct QueueEntry {

  QueueEntry(NBestState* s);

  NBestState* s_;
};

QueueEntry::QueueEntry(NBestState* s) : s_(s) {}

bool operator<(const QueueEntry& s1, const QueueEntry& s2)
{

  double p1 = s1.s_->cost_ + s1.s_->node_->forward_potential_;
  double p2 = s2.s_->cost_ + s2.s_->node_->forward_potential_;

  return (p1 > p2);
}

DAGNode::DAGNode(double forward) : forward_potential_(forward), nVisits_(0) {}

void DAGNode::add_edge(DAGEdge* edge)
{

  assert(edge->to_ == this);

  uint size = incoming_edge_.size();
  incoming_edge_.resize(size+1);
  incoming_edge_[size] = edge;
}

/************/

DAGEdge::DAGEdge(DAGNode* from, DAGNode* to, double weight, const Storage1D<uint>& output) :
  from_(from), to_(to), weight_(weight), output_(output)
{
  to_->add_edge(this);
}


/************/

DAG::DAG(uint nNodes, uint nEdges) : nUsedNodes_(0), nUsedEdges_(0)
{
  node_.resize(nNodes,0);
  edge_.resize(nEdges,0);
}

DAG::~DAG()
{

  for (uint e=0; e < nUsedEdges_; e++) {
    if (edge_[e] != 0)
      delete edge_[e];
  }

  for (uint n=0; n < nUsedNodes_; n++) {
    if (node_[n] != 0)
      delete node_[n];
  }
}

void DAG::clear()
{

  for (uint e=0; e < nUsedEdges_; e++) {
    if (edge_[e] != 0)
      delete edge_[e];
  }

  for (uint n=0; n < nUsedNodes_; n++) {
    if (node_[n] != 0)
      delete node_[n];
  }

  nUsedNodes_ = 0;
  nUsedEdges_ = 0;
}


uint DAG::add_node(double forward_potential)
{

  if (nUsedNodes_ == node_.size()) {
    node_.resize((node_.size()+1) * 1.8,0);
  }

  node_[nUsedNodes_] = new DAGNode(forward_potential);

  nUsedNodes_++;

  return (nUsedNodes_-1);
}

void DAG::add_edge(uint from, uint to, double weight, const Storage1D<uint>& output)
{

  if(nUsedEdges_ == edge_.size()) {
    edge_.resize((edge_.size()+1)*1.8);
  }

  edge_[nUsedEdges_] = new DAGEdge(node_[from],node_[to],weight,output);

  nUsedEdges_++;
}

DAGNode* DAG::node(uint num)
{
  return node_[num];
}


//note: this assumes that accurate forward potentials have been set
void DAG::nbest(uint N, uint start_node, uint end_node,
                std::vector<Storage1D<uint> >& sequence, std::vector<double>* score)
{

  //std::cerr << "A" << std::endl;

  for (uint k=0; k < nUsedNodes_; k++)
    node_[k]->nVisits_ = 0;

  //std::cerr << "B" << std::endl;

  std::vector<NBestState* > states;

  Storage1D<uint> empty_seq;

  //start-state
  states.push_back( new NBestState(node_[end_node],0.0,0,empty_seq)  );

  std::priority_queue<QueueEntry> queue;

  queue.push(QueueEntry(states[0]));

  double last_found = -1e300;

  while (!queue.empty()) {

    //std::cerr << "next iter" << std::endl;

    QueueEntry cur_entry = queue.top();
    queue.pop();

    NBestState* cur_state = cur_entry.s_;
    double cur_score = cur_state->cost_;

    DAGNode* cur_node = cur_state->node_;

    cur_node->nVisits_++;

    if (cur_node == node_[start_node]) {
      //new sequence found -> backtrace

      //std::cerr << "new sequence found" << std::endl;

      assert(cur_score >= last_found - 0.01);

      last_found = cur_score;

      std::vector<uint> cur_sequence;

      while (cur_state->prev_ != 0) {

        const Storage1D<uint>& out = cur_state->output_;

        if (out.size() != 0) {
          for (uint j=0; j < out.size(); j++)
            cur_sequence.push_back(out[j]);
        }

        cur_state = cur_state->prev_;
      }

      //reverse sequence
      sequence.push_back(Storage1D<uint>(cur_sequence.size()));
      for (uint k=0; k < cur_sequence.size(); k++)
        sequence.back()[k] = cur_sequence[k];

      if (score != 0)
        score->push_back(cur_score);

      if (cur_node->nVisits_ >= N)
        break;
    }

    if (cur_node->nVisits_ <= N) {
      for (uint k=0; k < cur_node->incoming_edge_.size(); k++) {

        DAGNode* from = cur_node->incoming_edge_[k]->from_;
        double edge_weight = cur_node->incoming_edge_[k]->weight_;

        states.push_back( new NBestState(from,cur_score + edge_weight, cur_state,
                                         cur_node->incoming_edge_[k]->output_) );

        queue.push(QueueEntry(states.back()));
      }
    }
  }

  //cleanup
  for (uint k=0; k < states.size(); k++) {
    delete states[k];
  }
}


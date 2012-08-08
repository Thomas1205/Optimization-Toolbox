/******* written by Thomas Schoenemann as a private person without employment, July 2011 *****/

#include "factorDualOpt.hh"
#include <map>

DualVariableNode::DualVariableNode(const Math1D::Vector<float>& cost) : cost_(cost) {
  dual_var_.resize(cost.size(),0);
}
  
void DualVariableNode::add_factor(DualFactorNode* node) {

  uint nPrevFactors = neighboring_factor_.size();
  neighboring_factor_.resize(nPrevFactors+1,node);
  dual_var_.resize(nLabels(),nPrevFactors+1,0.0); 
}

const double* DualVariableNode::get_dual_var_start() const {

  return dual_var_.direct_access();
}

double* DualVariableNode::get_dual_var_start() {

  return dual_var_.direct_access();
}

  
double* DualVariableNode::get_dual_vars(const DualFactorNode* node) {
  
  double* ptr = dual_var_.direct_access();

  const uint label_size = dual_var_.xDim();
  const uint nFactors = neighboring_factor_.size();

  bool found = false;
  for (uint i=0; i < nFactors; i++) {
    if (neighboring_factor_[i] == node) {
      found = true;
      break;
    }
    else
      ptr += label_size;
  }

  if (!found) {
    INTERNAL_ERROR << "node not found" << std::endl;
    exit(1);
  }

  return ptr;
}

const double* DualVariableNode::get_dual_vars(const DualFactorNode* node) const {
  
  const double* ptr = dual_var_.direct_access();

  const uint label_size = dual_var_.xDim();
  const uint nFactors = neighboring_factor_.size();

  bool found = false;
  for (uint i=0; i < nFactors; i++) {
    if (neighboring_factor_[i] == node) {
      found = true;
      break;
    }
    else
      ptr += label_size;
  }

  if (!found) {
    INTERNAL_ERROR << "node not found" << std::endl;
    exit(1);
  }

  return ptr;
}


  
uint DualVariableNode::nLabels() const {
  return dual_var_.xDim();
}

void DualVariableNode::init_dual_vars() {
  dual_var_.set_constant(0.0);
}

const Storage1D<DualFactorNode*>& DualVariableNode::neighboring_factor() const {
  return neighboring_factor_;
}


/**********************************/

void DualVariableNode::compute_message(const DualFactorNode* node, Math1D::Vector<double>& msg) const {

  const uint label_size = cost_.size();
  const uint nFactors = neighboring_factor_.size();

  msg.resize(label_size);
  
  for (uint l=0; l < label_size; l++)
    msg[l] = cost_[l];

  for (uint k=0; k < nFactors; k++) {

    if (neighboring_factor_[k] != node) {

      for (uint l=0; l < label_size; l++)
        msg[l] += dual_var_(l,k);
    }
  }
}

double DualVariableNode::cost(uint label) const {
  return cost_[label];
}

const Math1D::Vector<float>& DualVariableNode::cost() const {
  return cost_;
}

double DualVariableNode::dual_value(uint& arg_min) const {

  Math1D::Vector<double> msg;
  this->compute_message(0,msg);

  double min_val = 1e300;
  arg_min = MAX_UINT;

  for (uint k=0; k < msg.size(); k++) {
    if (msg[k] < min_val) {
      min_val = msg[k];
      arg_min = k;
    }
  }

  return min_val;
}

/**********************************/

DualFactorNode::DualFactorNode(const Storage1D<DualVariableNode*>& participating_vars)
  : participating_var_(participating_vars) {

  for (uint i=0; i < participating_var_.size(); i++)
    participating_var_[i]->add_factor(this);
}
  
const Storage1D<DualVariableNode*>& DualFactorNode::participating_nodes() const {
  return participating_var_;
}

/*virtual*/ double DualFactorNode::dual_value() const {

  Math1D::Vector<uint> labels;
  return compute_minimizer(labels);
}

/**********************************/

BinaryDualFactorNodeBase::BinaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars) :
  DualFactorNode(participating_vars) {
}

void BinaryDualFactorNodeBase::update_duals(const Math2D::Matrix<float>& cost, DualBCDMode mode) {

  //std::cerr << "bin" << std::endl;

  NamedStorage1D<Math1D::Vector<double> > msg(2, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();

  //important for MSD-mode
  msg[0].resize(nLabels1);
  msg[1].resize(nLabels2);

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());

  for (uint v=0; v < 2; v++) {
    if (mode == DUAL_BCD_MODE_MPLP)
      participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const double* msg0 = msg[0].direct_access();
  const double* msg1 = msg[1].direct_access();

  double* dp0 = dual_ptr[0];
  double* dp1 = dual_ptr[1];

  //for var 1
  if (mode == DUAL_BCD_MODE_MSD)
    participating_var_[0]->compute_message(this, msg[0]);

  for (uint l1=0; l1 < nLabels1; l1++) {

    double min_val = 1e300;

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      double hyp = cost(l1,l2);
      if (mode == DUAL_BCD_MODE_MPLP)
        hyp += msg0[l1] + msg1[l2];
      else
        hyp -= dp1[l2];

      if (hyp < min_val)
        min_val = hyp;
    }

    dp0[l1] = min_val;
  }

  for (uint l=0; l < nLabels1; l++) {
    if (mode == DUAL_BCD_MODE_MPLP) {
      dp0[l] *= 0.5;
      dp0[l] -= msg[0][l];
    }
    else {
      dp0[l] = 0.5 * (dp0[l] - msg0[l]);
    }
  }
  
  //for var 2
  if (mode == DUAL_BCD_MODE_MSD)
    participating_var_[1]->compute_message(this, msg[1]);

  for (uint l2=0; l2 < nLabels2; l2++) {

    double min_val = 1e300;

    for (uint l1=0; l1 < nLabels1; l1++) {

      double hyp = cost(l1,l2);
      if (mode == DUAL_BCD_MODE_MPLP) 
        hyp += msg[0][l1] + msg[1][l2];
      else
        hyp -= dp0[l1];

      if (hyp < min_val)
        min_val = hyp;      
    }

    dp1[l2] = min_val;
  }

  for (uint l=0; l < nLabels2; l++) {
    if (mode == DUAL_BCD_MODE_MPLP) {
      dp1[l] *= 0.5;
      dp1[l] -= msg[1][l];
    }
    else {
      dp1[l] = 0.5 * (dp1[l] - msg1[l]);
    }
  }
}
  
double BinaryDualFactorNodeBase::dual_value(const Math2D::Matrix<float>& cost) const {

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  uint nLabels1 = cost.xDim();
  uint nLabels2 = cost.yDim();
  
  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double min_cost = 1e300;
  
  for (uint l1=0; l1 < nLabels1; l1++) {
    for (uint l2=0; l2 < nLabels2; l2++) {
      double cur_cost = cost(l1,l2) - dual_ptr[0][l1] - dual_ptr[1][l2];
      if (cur_cost < min_cost)
        min_cost = cur_cost;
    }
  }

  return min_cost;
}

double BinaryDualFactorNodeBase::compute_minimizer(const Math2D::Matrix<float>& cost, Math1D::Vector<uint>& min_labels) const {

  min_labels.resize(2);

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  uint nLabels1 = cost.xDim();
  uint nLabels2 = cost.yDim();
  
  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double min_cost = 1e300;
  
  for (uint l1=0; l1 < nLabels1; l1++) {
    for (uint l2=0; l2 < nLabels2; l2++) {
      double cur_cost = cost(l1,l2) - dual_ptr[0][l1] - dual_ptr[1][l2];
      if (cur_cost < min_cost) {
        min_cost = cur_cost;
        min_labels[0] = l1;
        min_labels[1] = l2;
      }
    }
  }

  return min_cost;
}

/**********************************/

PottsDualFactorNode::PottsDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda) :
  DualFactorNode(participating_vars), lambda_(lambda) {
}

/*virtual*/ double PottsDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return (labels[0] == labels[1]) ? 0.0 : lambda_;
}

/*virtual*/ void PottsDualFactorNode::update_duals(DualBCDMode mode) {

  NamedStorage1D<Math1D::Vector<double> > msg(2, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();

  for (uint v=0; v < 2; v++) {
    if (mode == DUAL_BCD_MODE_MPLP)
      participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double* dp0 = dual_ptr[0];
  double* dp1 = dual_ptr[1];

  if (mode == DUAL_BCD_MODE_MPLP) {

    Math1D::Vector<double> msg_min(2); //note: already includes lambda

    for (uint k=0; k < 2; k++) {
      
      const Math1D::Vector<double>& cur_msg = msg[k];
      
      const uint cur_size = participating_var_[k]->nLabels();
      
      double cur_min = 1e300;
      for (uint l=0; l < cur_size; l++) {
        if (cur_msg[l] < cur_min)
          cur_min = cur_msg[l];
      }
      
      msg_min[k] = cur_min + lambda_;
    }      

    //NOTE: this code assumes that nLabels1 == nLabels2

    //Message 1
    for (uint l=0; l < nLabels1; l++) {
	
      dp0[l] = 0.5 * (std::min(msg_min[1], msg[1][l]) - msg[0][l]);
    }
    
    //Message 2
    for (uint k=0; k < nLabels2; k++) {

      dp1[k] = 0.5 * (std::min(msg_min[0], msg[0][k]) - msg[1][k]);
    }  
  }
  else {
    // MSD-mode
    //NOTE: this code assumes that nLabels1 == nLabels2

    //variable 1
    participating_var_[0]->compute_message(this, msg[0]);

    double min2 = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++) {
      if (-dp1[l2] < min2)
        min2 = -dp1[l2];
    }

    for (uint l1=0; l1 < nLabels1; l1++) {
      dp0[l1] = 0.5 * ( std::min(dp1[l1],min2+lambda_) - msg[0][l1]);
    }

    //variable 2
    participating_var_[1]->compute_message(this, msg[1]);

    double min1 = 1e300;
    for (uint l1=0; l1 < nLabels1; l1++) {
      if (-dp0[l1] < min1)
        min1 = -dp0[l1];
    }
    
    for (uint l2=0; l2 < nLabels2; l2++) {
      dp1[l2] = 0.5 * ( std::min(-dp0[l2],min1+lambda_) - msg[1][l2]);
    }
  }
}

/*virtual*/ double PottsDualFactorNode::dual_value() const {

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }
  
  Math1D::Vector<double> neg_dual_min(2);

  for (uint k=0; k < 2; k++) {

    double cur_min = 1e300;
    for (uint l=0; l < participating_var_[k]->nLabels(); l++) {
      double hyp = -dual_ptr[k][l];

      if (hyp < cur_min)
        cur_min = hyp;
    }
    
    neg_dual_min[k] = cur_min;
  }

  double min_cost = neg_dual_min.sum() + lambda_;

  for (uint l=0; l < std::min(nLabels1,nLabels2); l++) {

    double hyp =  -dual_ptr[0][l] - dual_ptr[1][l];

    if (hyp < min_cost)
      min_cost = hyp;
  }
  
  return min_cost;
}

/*virtual*/ double PottsDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  min_labels.resize(2);

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }
  
  Math1D::Vector<double> neg_dual_min(2);

  for (uint k=0; k < 2; k++) {

    double cur_min = 1e300;
    for (uint l=0; l < participating_var_[k]->nLabels(); l++) {
      double hyp = -dual_ptr[k][l];

      if (hyp < cur_min) {
        cur_min = hyp;
        min_labels[k] = l;
      }
    }
    
    neg_dual_min[k] = cur_min;
  }

  double min_cost = neg_dual_min.sum() + lambda_;

  for (uint l=0; l < std::min(nLabels1,nLabels2); l++) {

    double hyp =  -dual_ptr[0][l] - dual_ptr[1][l];

    if (hyp < min_cost) {
      min_cost = hyp;
      
      for (uint k=0; k < 2; k++)
        min_labels[k] = l;
    }
  }
  
  return min_cost;
}


/**********************************/

BinaryDualFactorNode::BinaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                           const Math2D::Matrix<float>& cost) :
  BinaryDualFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ void BinaryDualFactorNode::update_duals(DualBCDMode mode) {
  BinaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double BinaryDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1]);
}

/*virtual*/ double BinaryDualFactorNode::dual_value() const {
  return BinaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double BinaryDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return BinaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}

/**********************************/

BinaryDualRefFactorNode::BinaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                 const Math2D::Matrix<float>& cost) :
  BinaryDualFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ void BinaryDualRefFactorNode::update_duals(DualBCDMode mode) {
  BinaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double BinaryDualRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1]);
}

/*virtual*/ double BinaryDualRefFactorNode::dual_value() const {
  return BinaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double BinaryDualRefFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return BinaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}


/**********************************/

GenericDualFactorNode::GenericDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                             const VarDimStorage<float>& cost) :
  DualFactorNode(participating_vars), cost_(cost)
{
}

/*virtual*/ void GenericDualFactorNode::update_duals(DualBCDMode mode) {

  assert(mode == DUAL_BCD_MODE_MPLP);

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  for (uint v=0; v < nVars; v++) {
    participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    nLabels[v] = participating_var_[v]->nLabels();

    for (uint l=0; l < nLabels[v]; l++) 
      dual_ptr[v][l] = 1e300; 
  }

  Math1D::Vector<size_t> labeling(nVars,0);

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      cost += msg[v][labeling[v]];
    }

    for (uint v=0; v < nVars; v++) {

      if (cost < dual_ptr[v][labeling[v]])
        dual_ptr[v][labeling[v]] = cost; 
    }

    //increase labeling
    uint l;
    for (l=0; l < nVars; l++) {

      labeling[l] = (labeling[l] + 1) % nLabels[l];
      if (labeling[l] != 0)
        break;
    }

    if (l == nVars) //all zero after increase => cycle completed
      break;
  }

  for (uint v=0; v < nVars; v++) {
    for (uint l=0; l < nLabels[v]; l++) {
      dual_ptr[v][l] /= double(nVars);
      dual_ptr[v][l] -= msg[v][l];
    }
  }
}

/*virtual*/ double GenericDualFactorNode::dual_value() const {

  uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    nLabels[v] = participating_var_[v]->nLabels();
  }

  Math1D::Vector<size_t> labeling(nVars,0);

  double min_val = 1e300;

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      cost -= dual_ptr[v][labeling[v]];
    }

    if (cost < min_val)
      min_val = cost;

    //increase labeling
    uint l;
    for (l=0; l < nVars; l++) {

      labeling[l] = (labeling[l] + 1) % nLabels[l];
      if (labeling[l] != 0)
        break;
    }

    if (l == nVars) //all zero after increase => cycle completed
      break;
  }

  return min_val;
}

/*virtual*/ double GenericDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  min_labels.resize(nVars);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    nLabels[v] = participating_var_[v]->nLabels();
  }

  Math1D::Vector<size_t> labeling(nVars,0);

  double min_val = 1e300;

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      cost -= dual_ptr[v][labeling[v]];
    }

    if (cost < min_val) {
      min_val = cost;
      for (uint v=0; v < nVars; v++)
        min_labels[v] = labeling[v];
    }

    //increase labeling
    uint l;
    for (l=0; l < nVars; l++) {

      labeling[l] = (labeling[l] + 1) % nLabels[l];
      if (labeling[l] != 0)
        break;
    }

    if (l == nVars) //all zero after increase => cycle completed
      break;
  }

  return min_val;
}

/*virtual*/ double GenericDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  Math1D::Vector<size_t> size_t_labels(labels.size());

  for (uint k=0; k < labels.size(); k++)
    size_t_labels[k] = labels[k];

  return cost_(size_t_labels);
}

/**********************************/

TernaryDualFactorNodeBase::TernaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars)
  : DualFactorNode(participating_vars) {}

void TernaryDualFactorNodeBase::update_duals(const Math3D::Tensor<float>& cost, DualBCDMode mode) {

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();
  const uint nLabels3 = cost.zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());

  NamedStorage1D<Math1D::Vector<double> > msg(3, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {

    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);

    if (mode == DUAL_BCD_MODE_MPLP) {
      participating_var_[v]->compute_message(this, msg[v]);

      for (uint l=0; l < participating_var_[v]->nLabels() ; l++) {
        dual_ptr[v][l] = 1e300;
      }
    }
    else
      msg[v].resize(participating_var_[v]->nLabels());
  }

  if (mode == DUAL_BCD_MODE_MPLP) {
    for (uint l1=0; l1 < nLabels1; l1++) {

      const double inter_val =  msg[0][l1];

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        const double part_sum = inter_val + msg[1][l2];
	
        for (uint l3=0; l3 < nLabels3; l3++) {
	  
          const double hyp = cost(l1,l2,l3) + part_sum + msg[2][l3];
	  
          if (hyp < dual_ptr[0][l1])
            dual_ptr[0][l1] = hyp;
          if (hyp < dual_ptr[1][l2])
            dual_ptr[1][l2] = hyp;
          if (hyp < dual_ptr[2][l3])
            dual_ptr[2][l3] = hyp;
        }
      }
    }

    for (uint k=0; k < 3; k++) {
      for (uint l=0; l < participating_var_[k]->nLabels() ; l++) {
        dual_ptr[k][l] *= 1.0 / 3.0;
        dual_ptr[k][l] -= msg[k][l];
      }
    }
  }
  else {
    //MSD-Mode

    //var 1
    participating_var_[0]->compute_message(this, msg[0]);

    for (uint l1=0; l1 < nLabels1; l1++) {

      double min_cost = 1e300;

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        const double inter_val = dual_ptr[1][l2];

        for (uint l3=0; l3 < nLabels3; l3++) {

          double hyp = cost(l1,l2,l3) - inter_val - dual_ptr[2][l3];
	  
          if (hyp < min_cost)
            min_cost = hyp;
        }
      }
      
      dual_ptr[0][l1] = 0.5 * (min_cost - msg[0][l1]);
    }
    
    //var 2
    participating_var_[1]->compute_message(this, msg[1]);

    for (uint l2=0; l2 < nLabels2; l2++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double inter_val = dual_ptr[0][l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          double hyp = cost(l1,l2,l3) - inter_val - dual_ptr[2][l3];
	  
          if (hyp < min_cost)
            min_cost = hyp;
        }
      }

      dual_ptr[1][l2] = 0.5 * (min_cost - msg[1][l2]);
    }

    //var 3
    participating_var_[2]->compute_message(this, msg[2]);

    for (uint l3=0; l3 < nLabels3; l3++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double inter_val = dual_ptr[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          double hyp = cost(l1,l2,l3) - inter_val - dual_ptr[1][l2];
	  
          if (hyp < min_cost)
            min_cost = hyp;
        }
      }

      dual_ptr[2][l3] = 0.5 * (min_cost - msg[2][l3]);
    }
  }
}

double TernaryDualFactorNodeBase::dual_value(const Math3D::Tensor<float>& cost) const {

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();
  const uint nLabels3 = cost.zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3=0; l3 < nLabels3; l3++) {

        const double hyp = cost(l1,l2,l3) - dual_ptr[0][l1] - dual_ptr[1][l2] - dual_ptr[2][l3];

        if (hyp < min_val)
          min_val = hyp;
      }
    }
  }  

  return min_val;
}

double TernaryDualFactorNodeBase::compute_minimizer(const Math3D::Tensor<float>& cost, Math1D::Vector<uint>& min_labels) const {

  min_labels.resize(3);

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();
  const uint nLabels3 = cost.zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3=0; l3 < nLabels3; l3++) {

        const double hyp = cost(l1,l2,l3) - dual_ptr[0][l1] - dual_ptr[1][l2] - dual_ptr[2][l3];

        if (hyp < min_val) {
          min_val = hyp;
          min_labels[0] = l1;
          min_labels[1] = l2;
          min_labels[2] = l3;
        }
      }
    }
  }  

  return min_val;
}

/**********************************/

TernaryDualFactorNode::TernaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                             const Math3D::Tensor<float>& cost) :
  TernaryDualFactorNodeBase(participating_vars), cost_(cost) {}

/*virtual*/ void TernaryDualFactorNode::update_duals(DualBCDMode mode) {
  TernaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double TernaryDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1],labels[2]);
}

/*virtual*/ double TernaryDualFactorNode::dual_value() const {
  return TernaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double TernaryDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return TernaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}

/**********************************/

TernaryDualRefFactorNode::TernaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                   const Math3D::Tensor<float>& cost) :
  TernaryDualFactorNodeBase(participating_vars), cost_(cost) {}

/*virtual*/ void TernaryDualRefFactorNode::update_duals(DualBCDMode mode) {
  TernaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double TernaryDualRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1],labels[2]);
}

/*virtual*/ double TernaryDualRefFactorNode::dual_value() const {
  return TernaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double TernaryDualRefFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return TernaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}

/**********************************/

SecondDiffDualFactorNode::SecondDiffDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda) :
  DualFactorNode(participating_vars), lambda_(lambda) {
  assert(participating_vars.size() == 3);
}  

/*virtual*/
void SecondDiffDualFactorNode::update_duals(DualBCDMode mode) {

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();
  const uint nLabels3 = participating_var_[2]->nLabels();

  NamedStorage1D<Math1D::Vector<double> > msg(3, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {

    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);

    if (mode == DUAL_BCD_MODE_MPLP) {
      participating_var_[v]->compute_message(this, msg[v]);

      for (uint l=0; l < participating_var_[v]->nLabels() ; l++) {
        dual_ptr[v][l] = 1e300;
      }
    }
    else
      msg[v].resize(participating_var_[v]->nLabels());
  }

  if (mode == DUAL_BCD_MODE_MPLP) {

    Math1D::Vector<double> min_msg(3);
    for (uint v=0; v < 3; v++) {
      min_msg[v] = msg[v].min();
    }

    for (int l1=0; l1 < int(nLabels1); l1++) {

      double best = min_msg[1] + min_msg[2] + 3*lambda_;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = msg[1][l2] + msg[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = msg[1][l2] + msg[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      dual_ptr[0][l1] = best + msg[0][l1];
    }

    for (int l2=0; l2 < int(nLabels2); l2++) {

      double best = min_msg[0] + min_msg[2] + 3*lambda_;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = msg[0][l1] + msg[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = msg[0][l1] + msg[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      dual_ptr[1][l2] = best + msg[1][l2];
    }

    for (int l3=0; l3 < int(nLabels3); l3++) {

      double best = min_msg[0] + min_msg[1] + 3*lambda_;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {
        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = msg[0][l1] + msg[1][l2];
          }
          else if (abs(so_diff) <= 1) {
            hyp = msg[0][l1] + msg[1][l2] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      dual_ptr[2][l3] = best + msg[2][l3];
    }

    for (uint k=0; k < 3; k++) {
      for (uint l=0; l < participating_var_[k]->nLabels() ; l++) {
        dual_ptr[k][l] *= 1.0 / 3.0;
        dual_ptr[k][l] -= msg[k][l];
      }
    }
  }
  else {
    //MSD-Mode

    //var 1
    participating_var_[0]->compute_message(this, msg[0]);

    double base_cost = 3*lambda_;

    double mincost2 = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++) {

      if (-dual_ptr[1][l2] < mincost2)
        mincost2 = -dual_ptr[1][l2];
    }

    base_cost += mincost2;

    double mincost3 = 1e300;
    for (uint l3=0; l3 < nLabels3; l3++) {

      if (-dual_ptr[2][l3] < mincost3)
        mincost3 = -dual_ptr[2][l3];
    }

    base_cost += mincost3;

    for (int l1=0; l1 < int(nLabels1); l1++) {

      double min_cost = base_cost;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -dual_ptr[1][l2] - dual_ptr[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -dual_ptr[1][l2] - dual_ptr[2][l3] + lambda_;
          }

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }      
      
      dual_ptr[0][l1] = 0.5 * (min_cost - msg[0][l1]);
    }

    //var 2
    participating_var_[1]->compute_message(this, msg[1]);

    base_cost = 3*lambda_ + mincost3;

    double mincost1 = 1e300;
    for (uint l1=0; l1 < nLabels1; l1++) {

      if (-dual_ptr[0][l1] < mincost1)
        mincost1 = -dual_ptr[0][l1];
    }

    base_cost += mincost1;
    
    for (int l2=0; l2 < int(nLabels2); l2++) {

      double min_cost = base_cost;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -dual_ptr[0][l1] - dual_ptr[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -dual_ptr[0][l1] - dual_ptr[2][l3] + lambda_;
          }

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }

      dual_ptr[1][l2] = 0.5 * (min_cost - msg[1][l2]);      
    }    
    
    //var 3
    participating_var_[2]->compute_message(this, msg[2]);

    base_cost = 3*lambda_ + mincost1;

    mincost2 = 1e300;

    for (uint l2=0; l2 < nLabels2; l2++) {

      if (-dual_ptr[1][l2] < mincost2)
        mincost2 = -dual_ptr[1][l2];
    }

    base_cost += mincost2;

    for (int l3=0; l3 < int(nLabels3); l3++) {

      double min_cost = base_cost;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {
        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -dual_ptr[0][l1] - dual_ptr[1][l2];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -dual_ptr[0][l1] - dual_ptr[1][l2] + lambda_;
          }

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }

      dual_ptr[2][l3] = 0.5 * (min_cost - msg[2][l3]);
    }
  }

}

/*virtual*/ 
double SecondDiffDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  int diff1 = labels[1] - labels[0];
  int diff2 = labels[2] - labels[1];
  
  int so_diff = diff2 - diff1;
  
  if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
    return 0.0; //no cost
  else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
    return lambda_;

  return 3*lambda_;
}

/*virtual*/ 
double SecondDiffDualFactorNode::dual_value() const {

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double best = 3*lambda_;
  for (uint v=0; v < 3; v++) {
    double cur_min = 1e300;

    for (uint l=0; l < participating_var_[v]->nLabels(); l++) {
      if (-dual_ptr[v][l] < cur_min)
        cur_min = -dual_ptr[v][l];
    }

    best += cur_min;
  }

  int nLabels1 = participating_var_[0]->nLabels();
  int nLabels2 = participating_var_[1]->nLabels();
  int nLabels3 = participating_var_[2]->nLabels();

  for (int l1=0; l1 < nLabels1; l1++) {

    double best_inter = 1e300;

    for (int l2=std::max(0,l1-1); l2 < std::min(l1+1,nLabels2-1); l2++) {
      for (int l3=std::max(0,l2-1); l3 < std::min(l2+1,nLabels3-1); l3++) {
      
        assert(abs(l2-l1) <= 1);
        assert(abs(l3-l2) <= 1);

        int so_diff = l3 - 2*l2 + l1;
	
        double hyp = 1e300;
	
        if (so_diff == 0) {
          hyp = -dual_ptr[1][l2] - dual_ptr[2][l3];
        }
        else if (abs(so_diff) <= 1) {
          hyp = -dual_ptr[1][l2] - dual_ptr[2][l3] + lambda_;
        }

        if (hyp < best_inter)
          best_inter = hyp;
      }
    }

    if (best_inter - dual_ptr[0][l1] < best)
      best = best_inter - dual_ptr[0][l1];
  }

  return best;
}

/*virtual*/
double SecondDiffDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  min_labels.resize(3);

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double best = 3*lambda_;
  for (uint v=0; v < 3; v++) {
    double cur_min = 1e300;

    for (uint l=0; l < participating_var_[v]->nLabels(); l++) {
      if (-dual_ptr[v][l] < cur_min) {
        cur_min = -dual_ptr[v][l];
        min_labels[v] = l;
      }
    }

    best += cur_min;
  }

  int nLabels1 = participating_var_[0]->nLabels();
  int nLabels2 = participating_var_[1]->nLabels();
  int nLabels3 = participating_var_[2]->nLabels();

  for (int l1=0; l1 < nLabels1; l1++) {

    double best_inter = 1e300;

    int opt_l2 = MAX_UINT;
    int opt_l3 = MAX_UINT;


    for (int l2=std::max(0,l1-1); l2 < std::min(l1+1,nLabels2-1); l2++) {
      for (int l3=std::max(0,l2-1); l3 < std::min(l2+1,nLabels3-1); l3++) {
      
        assert(abs(l2-l1) <= 1);
        assert(abs(l3-l2) <= 1);

        int so_diff = l3 - 2*l2 + l1;
	
        double hyp = 1e300;
	
        if (so_diff == 0) {
          hyp = -dual_ptr[1][l2] - dual_ptr[2][l3];
        }
        else if (abs(so_diff) <= 1) {
          hyp = -dual_ptr[1][l2] - dual_ptr[2][l3] + lambda_;
        }

        if (hyp < best_inter) {
          best_inter = hyp;
          opt_l2 = l2;
          opt_l3 = l3;
        }
      }
    }

    if (best_inter - dual_ptr[0][l1] < best) {
      best = best_inter - dual_ptr[0][l1];
      min_labels[0] = l1;
      min_labels[1] = opt_l2;
      min_labels[2] = opt_l3;
    }
  }

  return best;
}


/**********************************/

FourthOrderDualFactorNodeBase::FourthOrderDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars) :
  DualFactorNode(participating_vars) {
}

void FourthOrderDualFactorNodeBase::update_duals(const Storage1D<Math3D::Tensor<float> >& cost, DualBCDMode mode) {

  const uint nLabels1 = cost.size();
  const uint nLabels2 = cost[0].xDim();
  const uint nLabels3 = cost[0].yDim();
  const uint nLabels4 = cost[0].zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());
  assert(nLabels4 == participating_var_[3]->nLabels());

  NamedStorage1D<Math1D::Vector<double> > msg(4, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(4, MAKENAME(dual_ptr));

  for (uint v=0; v < 4; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);

    if (mode == DUAL_BCD_MODE_MPLP) {
      participating_var_[v]->compute_message(this, msg[v]);
      for (uint l=0; l < participating_var_[v]->nLabels() ; l++) {
        dual_ptr[v][l] = 1e300;
      }
    }
    else
      msg[v].resize(participating_var_[v]->nLabels());
  }

  if (mode == DUAL_BCD_MODE_MPLP) {
    for (uint l1=0; l1 < nLabels1; l1++) {

      double* cur_dp0 = dual_ptr[0] + l1;

      const Math3D::Tensor<float>& cur_cost = cost[l1];

      const double inter1 = msg[0][l1];

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        double* cur_dp1 = dual_ptr[1] + l2;
	
        const double part_sum1 = inter1 + msg[1][l2];
	
        for (uint l3=0; l3 < nLabels3; l3++) {
	  
          const double part_sum2 = part_sum1 + msg[2][l3];
	  
          for (uint l4=0; l4 < nLabels4; l4++) {
	    
            const double hyp = cur_cost(l2,l3,l4) + part_sum2 + msg[3][l4];
	    
            if (hyp < (*cur_dp0))
              *cur_dp0 = hyp;
            if (hyp < (*cur_dp1))
              *cur_dp1 = hyp;
            if (hyp < dual_ptr[2][l3])
              dual_ptr[2][l3] = hyp;
            if (hyp < dual_ptr[3][l4])
              dual_ptr[3][l4] = hyp;
          }
        }
      }
    }

    for (uint k=0; k < 4; k++) {

      double* cur_ptr = dual_ptr[k];

      for (uint l=0; l < participating_var_[k]->nLabels() ; l++) {
        cur_ptr[l] *= 0.25;
        cur_ptr[l] -= msg[k][l];
      }
    }
  }
  else {
    
    //var 1
    participating_var_[0]->compute_message(this, msg[0]);

    for (uint l1=0; l1 < nLabels1; l1++) {

      double min_cost = 1e300;

      const Math3D::Tensor<float>& cur_cost = cost[l1];

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        const double inter1 = dual_ptr[1][l2];

        for (uint l3=0; l3 < nLabels3; l3++) {

          const double inter2 = inter1 + dual_ptr[2][l3];

          for (uint l4=0; l4 < nLabels4; l4++) {
	    
            double hyp = cur_cost(l2,l3,l4) - inter2 - dual_ptr[3][l4];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }
      
      dual_ptr[0][l1] = 0.5 * (min_cost - msg[0][l1]);
    }

    //var 2
    participating_var_[1]->compute_message(this, msg[1]);

    for (uint l2=0; l2 < nLabels2; l2++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = dual_ptr[0][l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          const double inter2 = inter1 + dual_ptr[2][l3];

          for (uint l4=0; l4 < nLabels4; l4++) {

            double hyp = cur_cost(l2,l3,l4) - inter2 - dual_ptr[3][l4];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }

      dual_ptr[1][l2] = 0.5 * (min_cost - msg[1][l2]);
    }

    //var 3
    participating_var_[2]->compute_message(this, msg[2]);
    
    for (uint l3=0; l3 < nLabels3; l3++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = dual_ptr[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          const double inter2 = inter1 + dual_ptr[1][l2];

          for (uint l4=0; l4 < nLabels4; l4++) {

            double hyp = cur_cost(l2,l3,l4) - inter2 - dual_ptr[3][l4];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }

      dual_ptr[2][l3] = 0.5 * (min_cost - msg[2][l3]);
    }

    //var 4
    participating_var_[3]->compute_message(this, msg[3]);
    
    for (uint l4=0; l4 < nLabels4; l4++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = dual_ptr[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          const double inter2 = inter1 + dual_ptr[1][l2];

          for (uint l3=0; l3 < nLabels3; l3++) {

            double hyp = cur_cost(l2,l3,l4) - inter2 - dual_ptr[2][l3];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }

      dual_ptr[3][l4] = 0.5 * (min_cost - msg[3][l4]);
    }
  }
}

double FourthOrderDualFactorNodeBase::dual_value(const Storage1D<Math3D::Tensor<float> >& cost) const {

  NamedStorage1D<const double*> dual_ptr(4, MAKENAME(dual_ptr));

  for (uint v=0; v < 4; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.size();
  const uint nLabels2 = cost[0].xDim();
  const uint nLabels3 = cost[0].yDim();
  const uint nLabels4 = cost[0].zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());
  assert(nLabels4 == participating_var_[3]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    const Math3D::Tensor<float>& cur_cost = cost[l1];

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3=0; l3 < nLabels3; l3++) {
	
        for (uint l4=0; l4 < nLabels4; l4++) {
	  
          const double hyp = cur_cost(l2,l3,l4) - dual_ptr[0][l1] - dual_ptr[1][l2] - dual_ptr[2][l3] - dual_ptr[3][l4];
	  
          if (hyp < min_val)
            min_val = hyp;
        }
      }
    }
  }  

  return min_val;
}

double FourthOrderDualFactorNodeBase::compute_minimizer(const Storage1D<Math3D::Tensor<float> >& cost, 
                                                        Math1D::Vector<uint>& min_labels) const {

  min_labels.resize(4);

  NamedStorage1D<const double*> dual_ptr(4, MAKENAME(dual_ptr));

  for (uint v=0; v < 4; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.size();
  const uint nLabels2 = cost[0].xDim();
  const uint nLabels3 = cost[0].yDim();
  const uint nLabels4 = cost[0].zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());
  assert(nLabels4 == participating_var_[3]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    const Math3D::Tensor<float>& cur_cost = cost[l1];

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3=0; l3 < nLabels3; l3++) {
	
        for (uint l4=0; l4 < nLabels4; l4++) {
	  
          const double hyp = cur_cost(l2,l3,l4) - dual_ptr[0][l1] - dual_ptr[1][l2] - dual_ptr[2][l3] - dual_ptr[3][l4];
	  
          if (hyp < min_val) {
            min_val = hyp;

            min_labels[0] = l1;
            min_labels[1] = l2;
            min_labels[2] = l3;
            min_labels[3] = l4;
          }
        }
      }
    }
  }  

  return min_val;
}


/**********************************/

FourthOrderDualFactorNode::FourthOrderDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                     const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderDualFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ void FourthOrderDualFactorNode::update_duals(DualBCDMode mode) {
  FourthOrderDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double FourthOrderDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_[labels[0]](labels[1],labels[2],labels[3]);
}

/*virtual*/ double FourthOrderDualFactorNode::dual_value() const {
  return FourthOrderDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double FourthOrderDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return FourthOrderDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}


/**********************************/

FourthOrderDualRefFactorNode::FourthOrderDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                           const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderDualFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ void FourthOrderDualRefFactorNode::update_duals(DualBCDMode mode) {
  FourthOrderDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double FourthOrderDualRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_[labels[0]](labels[1],labels[2],labels[3]);
}

/*virtual*/ double FourthOrderDualRefFactorNode::dual_value() const {
  return FourthOrderDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double FourthOrderDualRefFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return FourthOrderDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}


/**********************************/

OneOfNDualFactorNode::OneOfNDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars) :
  DualFactorNode(participating_vars) {}

/*virtual*/ double OneOfNDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return (labels.sum() == 1) ? 0.0 : 1e75;
}

/*virtual*/ void OneOfNDualFactorNode::update_duals(DualBCDMode mode) {

  assert(mode == DUAL_BCD_MODE_MPLP);

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Math1D::Vector<double> rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    rel_msg[k] = msg[k][1] - msg[k][0];
    offs += msg[k][0];
  }

  double best = 1e15;
  uint arg_best = MAX_UINT;
  double second_best = 1e15;
  
  for (uint k=0; k < nVars; k++) {

    if (rel_msg[k] < best) {

      second_best = best;

      best = rel_msg[k];
      arg_best = k;
    }
    else if (rel_msg[k] < second_best) {

      second_best = rel_msg[k];
    }
  }

  best += offs;
  best /= nVars;
  second_best += offs;
  second_best /= nVars;

  for (uint k=0; k < nVars; k++) {

    //msg0
    if (k == arg_best)
      dual_ptr[k][0] = second_best - msg[k][0];
    else
      dual_ptr[k][0] = best - msg[k][0];
    
    //msg 1
    dual_ptr[k][1] = ((offs + rel_msg[k]) / nVars) - msg[k][1];
  }
}

/*virtual*/ double OneOfNDualFactorNode::dual_value() const {

  uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  double dual_zero_sum = 0.0;

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    dual_zero_sum -= dual_ptr[v][0];
  }

  double min_val = 1e300;

  for (uint v=0; v < nVars; v++) {

    double hyp = dual_zero_sum + dual_ptr[v][0] - dual_ptr[v][1];

    if (hyp < min_val)
      min_val = hyp;
  }
  
  return min_val;
}

/*virtual*/ double OneOfNDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  min_labels.resize(nVars);
  min_labels.set_constant(0);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  double dual_zero_sum = 0.0;

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    dual_zero_sum -= dual_ptr[v][0];
  }

  double min_val = 1e300;
  uint arg_min = MAX_UINT;

  for (uint v=0; v < nVars; v++) {

    double hyp = dual_zero_sum + dual_ptr[v][0] - dual_ptr[v][1];

    if (hyp < min_val) {
      min_val = hyp;
      arg_min = v;
    }
  }
  
  min_labels[arg_min] = 1;
  return min_val;
}

/**********************************/

CardinalityDualFactorNode::CardinalityDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                                                     const Math1D::Vector<float>& card_cost) :
  DualFactorNode(participating_vars), card_cost_(card_cost) {
}

/*virtual*/ double CardinalityDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  
  uint sum = labels.sum();
  return card_cost_[sum];
}

/*virtual*/ void CardinalityDualFactorNode::update_duals(DualBCDMode mode) {

  assert(mode == DUAL_BCD_MODE_MPLP);

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }
  

  //TEMPORARY
  // for (uint idx=0; idx < nVars; idx++) {

  //   Math1D::NamedVector<double> rel_param(nVars-1,MAKENAME(rel_param));
    
  //   double offs = 0.0;
    
  //   uint next = 0;
  //   for (uint k=0; k < nVars; k++) {
      
  //     if (k != idx) {
  // 	rel_param[next] = (msg[k][1]) - (msg[k][0]);
  // 	offs += msg[k][0];
	
  // 	next++;
  //     }
  //   }

  //   std::sort(rel_param.direct_access(), rel_param.direct_access() + nVars-1);
  
  //   double cum_sum = 0.0;
    
  //   Math1D::Vector<double> message(2,1e300);

  //   for (uint c=0; c < nVars; c++) {
      
  //     double hyp0 = cum_sum + card_cost_[c];
  //     if (hyp0 < message[0])
  // 	message[0] = hyp0;
      
  //     double hyp1 = cum_sum + card_cost_[c+1];
  //     if (hyp1 < message[1])
  // 	message[1] = hyp1;
      
  //     if (c+1 < nVars) 
  // 	cum_sum += rel_param[c];
  //   }

  //   for (uint l=0; l < 2; l++)
  //     message[l] += offs + msg[idx][l];

  //   for (uint l=0; l < 2; l++) {

  //     dual_ptr[idx][l] = (message[l] / double(nVars)) - msg[idx][l]; 
  //   }
  // }
  //END_TEMPORARY


  Math1D::Vector<std::pair<double,uint> > rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    rel_msg[k].first = msg[k][1] - msg[k][0];
    rel_msg[k].second = k;
    offs += msg[k][0];
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);
  
  Math1D::Vector<uint> order(nVars);
  for (uint k=0; k < nVars; k++) {
    order[rel_msg[k].second] = k + 1;
  }

  Math1D::Vector<double> cum_sum(nVars+1);
  cum_sum[0] = rel_msg[0].first;
  for (uint k=1; k < nVars; k++) {
    assert(rel_msg[k].first >= rel_msg[k-1].first);
    cum_sum[k] = cum_sum[k-1] + rel_msg[k].first;
  }

  Math1D::Vector<double> cum_best(nVars + 1);
  cum_best[0] = card_cost_[0]; 
  
  for (uint k=1; k <= nVars; k++) {

    double hyp = card_cost_[k] + cum_sum[k-1];
    cum_best[k] = std::min(hyp, cum_best[k-1]);
  }


  Math1D::Vector<double> cum_best_m1(nVars+1);
  cum_best_m1[0] = 1e300;
  cum_best_m1[1] = card_cost_[1];

  for (uint k=2; k <= nVars; k++) {
    cum_best_m1[k] = std::min(cum_best_m1[k-1], card_cost_[k] + cum_sum[k-2]  ); 
  }

  Math1D::Vector<double> rev_cum_best(nVars + 1);
  rev_cum_best[nVars] = card_cost_[nVars] + cum_sum[nVars-1];
  for (int k=nVars-1; k >= 1; k--) {

    double hyp = card_cost_[k] + cum_sum[k-1];
    rev_cum_best[k] = std::min(hyp, rev_cum_best[k+1]);
  }
  rev_cum_best[0] = 1e300;

  Math1D::Vector<double> rev_cum_best_p1(nVars+1);
  rev_cum_best_p1[nVars] = 1e300;
  for (int k=nVars-1; k >= 0; k--) {
    rev_cum_best_p1[k] = std::min(rev_cum_best_p1[k+1], card_cost_[k] + cum_sum[k]);
  }

  for (uint k=0; k < nVars; k++) {

    const uint cur_order = order[k];

    double cur_rel_msg = msg[k][1] - msg[k][0];

    const double val0 = std::min(cum_best[cur_order-1], rev_cum_best_p1[cur_order] - cur_rel_msg);
    dual_ptr[k][0] = val0 + offs;

    const double val1 = std::min(rev_cum_best[cur_order], cum_best_m1[cur_order-1] + cur_rel_msg );
    dual_ptr[k][1] = val1 + offs;

    //DEBUG
#if 0
    //a) check msg0
    double best0 = card_cost_[0];
    uint arg_best0 = 0;
    for (uint l=1; l < nVars; l++) {

      double hyp = card_cost_[l];
      if (cur_order > l)
        hyp += cum_sum[l-1];
      else {
        hyp += cum_sum[l] - cur_rel_msg;
      }

      if (hyp < best0) {
        best0 = hyp;
        arg_best0 = l;
      }
    }
    if ( ! (fabs(best0-val0) < 1e-2)) {
      std::cerr << "best0: " << best0 << ", msg0: " << val0 << ", arg best: "  << arg_best0 
                << ", order: " << cur_order << ", factor size: " << nVars << std::endl;

      std::cerr << "cum_best value: " << cum_best[cur_order-1] << std::endl;

      std::cerr << "p1 value: " << (rev_cum_best_p1[cur_order] - cur_rel_msg) << std::endl;

      std::cerr << "card_cost[arg0]: " << card_cost_[arg_best0] << std::endl;
      std::cerr << "cum_sum[arg0]: " << cum_sum[arg_best0] << std::endl;
      std::cerr << "rel msg: " << cur_rel_msg << std::endl;
    }
    assert(fabs(best0-val0) < 1e-2);

    //b) check msg1
    double best1 = 1e300;
    uint arg_best1 = 0;
    for (uint l=1; l <= nVars; l++) {

      double hyp = card_cost_[l];
      if (cur_order <= l)
        hyp += cum_sum[l-1] - cur_rel_msg;
      else {
        if (l >= 2)
          hyp += cum_sum[l-2];
      }

      if (hyp < best1) {
        best1 = hyp;
        arg_best1 = l;
      }
    }
    
    if ( ! (fabs(best1-val1) < 1e-2)) {

      std::cerr << "best1: " << best1 << ", msg1: " << val1 << ", arg_best1: " << arg_best1 
                << ", factor size: " << nVars << std::endl;
      std::cerr << "rev_cum_best: " << rev_cum_best[cur_order] << std::endl;
      std::cerr << "cum_best_m1 value: " << cum_best_m1[cur_order-1] << std::endl;


      std::cerr << "cur_order: " << cur_order << std::endl;
      std::cerr << "cur_rel_msg: " << cur_rel_msg << std::endl;
    }

    assert(fabs(best1-val1) < 1e-2);
#endif
    //END_DEBUG

    for (uint i=0; i < 2; i++) {
      dual_ptr[k][i] /= double(nVars);
      dual_ptr[k][i] -= msg[k][i];
    }
  }
}

/*virtual*/ double CardinalityDualFactorNode::dual_value() const {

  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Math1D::Vector<double> rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    const double val0 = - dual_ptr[k][0];
    const double val1 = - dual_ptr[k][1];

    rel_msg[k] = val1 - val0;
    offs += val0;
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);

  double min_val = card_cost_[0];

  double cum_sum = 0.0;
  for (uint c=1; c <= nVars; c++) {
    cum_sum += rel_msg[c-1];

    const double hyp = cum_sum + card_cost_[c];

    if (hyp < min_val) {
      min_val = hyp;
    }
  }

  return min_val + offs;
}

/*virtual*/ double CardinalityDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  min_labels.resize(nVars);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Storage1D<std::pair<double,uint> > rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    const double val0 = - dual_ptr[k][0];
    const double val1 = - dual_ptr[k][1];

    rel_msg[k] = std::make_pair(val1 - val0,k);
    offs += val0;
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);

  double min_val = card_cost_[0];
  uint arg_min = 0;
  
  double cum_sum = 0.0;
  for (uint c=1; c <= nVars; c++) {
    cum_sum += rel_msg[c-1].first;
    
    const double hyp = cum_sum + card_cost_[c];

    if (hyp < min_val) {
      min_val = hyp;
      arg_min = c;
    }
  }

  min_labels.set_constant(0);
  for (uint k=0; k < arg_min; k++)
    min_labels[rel_msg[k].second] = 1;

  return min_val + offs;
}

/**********************************/

BILPConstraintDualFactorNode::BILPConstraintDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                                                           const Storage1D<bool>& positive, int rhs_lower, int rhs_upper) :
  DualFactorNode(participating_vars), positive_(positive), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

  assert(rhs_lower_ <= rhs_upper_);

  int nPositive = 0;
  int nNegative = 0;

  for (uint k=0; k < positive_.size(); k++) {
    if (positive_[k])
      nPositive++;
    else
      nNegative++;
  }

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  if (rhs_upper_ + zero_offset < 0 || rhs_lower + zero_offset >= range) {
    INTERNAL_ERROR << "constraint is unsatisfiable. Exiting..." << std::endl;
    exit(1);
  }

  //adjust the bounds to the actually possible range of values
  if (rhs_lower_ + zero_offset < 0) {
    rhs_lower_ -= (rhs_lower_ + zero_offset);
  }
  if (rhs_upper_ + zero_offset >= range) {
    rhs_upper_ -= (rhs_upper_ + zero_offset - range +1);
  }
}

/*virtual*/ double BILPConstraintDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  int sum = 0;

  for (uint k=0; k < participating_var_.size(); k++) {

    const int label = labels[k];
    if (positive_[k])
      sum += label;
    else
      sum -= label;
  }

  return (sum == 0) ? 0.0 : 1e15;
}

/*virtual*/ void BILPConstraintDualFactorNode::update_duals(DualBCDMode mode) {

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {

    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);

    if (mode == DUAL_BCD_MODE_MPLP)
      participating_var_[v]->compute_message(this, msg[v]);
    else
      msg[v].resize(participating_var_[v]->nLabels());
  }

  int nPositive = 0;
  int nNegative = 0;

  for (uint k=0; k < nVars; k++) {
    if (positive_[k])
      nPositive++;
    else
      nNegative++;
  }

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  //some intermediate values can never result in a satisfied constraint 
  // => prune the range to the useful values
  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  Math3D::Tensor<double> forward(range,2,nVars);
  Math2D::Matrix<double> forward_light(range,nVars);
    
  Math2D::Matrix<double> backward_light(range,nVars);

  const uint last_var = nVars-1;

  if (mode == DUAL_BCD_MODE_MPLP) {

    /**** forward ****/

    //init
    for (int sum=0; sum < range; sum++) {
      
      forward_light(sum,0) = 1e100;
      for (int l=0; l < 2; l++) {
        forward(sum,l,0) = 1e100;
      }
    }

    forward(zero_offset,0,0) = msg[0][0];
    forward_light(zero_offset,0) = msg[0][0];
    const int init_mul = (positive_[0]) ? 1 : -1;

    if (int(zero_offset)+init_mul >= 0
	&& int(zero_offset)+init_mul < range ) {
      forward(zero_offset+init_mul,1,0) = msg[0][1];
      forward_light(zero_offset+init_mul,0) = msg[0][1];
    }

    //proceed
    for (uint v=1; v < nVars; v++) {
      
      for (int sum=0; sum < range; sum++) {

        for (int l=0; l < 2; l++) {
	
          double best_prev = 1e75;
	
          int move = l;
          if (positive_[v]) //since we are tracing backward here
            move *= -1;
	  
          const int dest = sum + move;
          if (dest >= 0 && dest < range) {
	    
            best_prev = forward_light(dest,v-1);
          }
	  
          forward(sum,l,v) = best_prev + msg[v][l];
        }
        forward_light(sum,v) = std::min(forward(sum,0,v), forward(sum,1,v));
      }
    }

    /**** backward ****/

    //init
    for (int sum=0; sum < range; sum++) 
      backward_light(sum,last_var) = 1e100;
    
    backward_light(zero_offset,last_var) = msg[last_var][0];
    const int end_mul = (positive_[last_var]) ? 1 : -1;
    if (int(zero_offset)+end_mul >= 0
	&& int(zero_offset)+end_mul < range ) {
      backward_light(zero_offset + end_mul,last_var) = msg[last_var][1];
    }
    
    //proceed
    for (int v=last_var-1; v >= 1; v--) {
      
      for (int sum=0; sum < range; sum++) {
	
        double best_prev = 1e75;
	
        for (int l=0; l < 2; l++) {
	  
          int move = l;
          if (positive_[v]) //since we are tracing backward here
            move *= -1;
	  
          const int dest = sum + move;
          double hyp = 1e75;
          if (dest >= 0 && dest < range) {
            hyp = backward_light(dest,v+1) + msg[v][l];
          }

          if (hyp < best_prev)
            best_prev = hyp;
        }
	
        backward_light(sum,v) = best_prev;
      }
    }
    
    /*** derive messages ***/
    
    for (uint k=0; k < nVars; k++) {

      for (uint l=0; l < 2; l++) {

        double min_msg = 1e300;

        for (int s=0; s < (int) range; s++) {

          double hyp = forward(s,l,k);

          if (k+1 < positive_.size()) {

            double best_bwd = 1e300;
	    
            const int diff = (s - zero_offset);
	    
            for (int r=rhs_lower_; r <= rhs_upper_; r++) {
              const int other = r + zero_offset - diff; 
	      
              if (other >= 0 && other < (int) range) {

                best_bwd = std::min(best_bwd,backward_light(other,k+1));
              }
            }

            hyp += best_bwd;
          }
          else {
            if (s < int(rhs_lower_ + zero_offset) || s > int(rhs_upper_ + zero_offset)) 
              hyp = 1e300;
          }
	  
          if (hyp < min_msg)
            min_msg = hyp;
        }

        assert(!isnan(min_msg));
	
        dual_ptr[k][l] = (min_msg) / nVars - msg[k][l];
      }
    }
  }
  else {
    //min sum diffusion mode

    //a) compute backward completely (and just once)
    //init
    for (int sum=0; sum < range; sum++) 
      backward_light(sum,last_var) = 1e100;
    
    backward_light(zero_offset,last_var) = -dual_ptr[last_var][0];
    const int end_mul = (positive_[last_var]) ? 1 : -1;
    backward_light(zero_offset + end_mul,last_var) = -dual_ptr[last_var][1];
    
    //proceed
    for (int v=last_var-1; v >= 0; v--) {
      
      for (int sum=0; sum < range; sum++) {
	
        double best_prev = 1e75;
	
        for (int l=0; l < 2; l++) {
	  
          int move = l;
          if (positive_[v]) //since we are tracing backward here
            move *= -1;
	  
          const int dest = sum + move;
          double hyp = 1e75;
          if (dest >= 0 && dest < range) {
            hyp = backward_light(dest,v+1) - dual_ptr[v][l];
          }

          if (hyp < best_prev)
            best_prev = hyp;
        }
	
        backward_light(sum,v) = best_prev;
      }
    }

    //b) compute forward incrementally and derive messages
    
    for (uint v=0; v < nVars; v++) {
      
      //correct?
      participating_var_[v]->compute_message(this, msg[v]);

      if (v == 0) {
        for (int sum=0; sum < range; sum++) {
      
          forward_light(sum,0) = 1e100;
          for (int l=0; l < 2; l++) {
            forward(sum,l,0) = 1e100;
          }
        }
	
        forward(zero_offset,0,0) = 0.0;
        forward_light(zero_offset,0) = 0.0;
        const int init_mul = (positive_[0]) ? 1 : -1;
        forward(zero_offset+init_mul,1,0) = 0.0;
        forward_light(zero_offset+init_mul,0) = 0.0;
      }
      else {

        for (int sum=0; sum < range; sum++) {

          for (int l=0; l < 2; l++) {
	    
            double best_prev = 1e75;
	    
            int move = l;
            if (positive_[v]) //since we are tracing backward here
              move *= -1;
	    
            const int dest = sum + move;
            if (dest >= 0 && dest < range) {
	      
              best_prev = forward_light(dest,v-1);
            }
	    
            forward(sum,l,v) = best_prev;
          }
          //note: forward_light is set below
        }
      }

      //now compute new duals
      for (uint l=0; l < 2; l++) {

        double min_msg = 1e300;

        for (int s=0; s < (int) range; s++) {

          double hyp = forward(s,l,v);

          if (v+1 < positive_.size()) {

            double best_bwd = 1e300;
	    
            const int diff = (s - zero_offset);
	    
            for (int r=rhs_lower_; r <= rhs_upper_; r++) {
              const int other = r + zero_offset - diff; 
	      
              if (other >= 0 && other < (int) range) {

                best_bwd = std::min(best_bwd,backward_light(other,v+1));
              }
            }

            hyp += best_bwd;
          }
          else {
            if (s < int(rhs_lower_ + zero_offset) || s > int(rhs_upper_ + zero_offset)) 
              hyp = 1e300;
          }
	  
          if (hyp < min_msg)
            min_msg = hyp;
        }
	
        assert(!isnan(min_msg));


        dual_ptr[v][l] = 0.5 * (min_msg - msg[v][l]);
      }


      //correct the freshly computed forward term
      for (int sum=0; sum < range; sum++) {

        for (uint l=0; l < 2; l++)
          forward(sum,l,v) -= dual_ptr[v][l];
        forward_light(sum,v) = std::min(forward(sum,0,v), forward(sum,1,v));
      }
    }

  }
}

/*virtual*/ double BILPConstraintDualFactorNode::dual_value() const {


  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  int nPositive = 0;
  int nNegative = 0;

  for (uint k=0; k < nVars; k++) {
    if (positive_[k])
      nPositive++;
    else
      nNegative++;
  }

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  //some intermediate values can never result in a satisfied constraint 
  // => prune the range to the useful values
  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  /**** forward ****/

  Math2D::Matrix<double> forward_light(range,participating_var_.size(),1e100);

  //init
  forward_light(zero_offset,0) = - dual_ptr[0][0];
  int init_mul = (positive_[0]) ? 1 : -1;
  forward_light(zero_offset+init_mul,0) = - dual_ptr[0][1];

  //proceed
  for (uint v=1; v < nVars; v++) {

    for (int sum=0; sum < range; sum++) {

      double best_val = 1e75;

      for (int l=0; l < 2; l++) {
	
        double best_prev = 1e75;
	
        int move = l;
        if (positive_[v]) //since we are tracing backward here
          move *= -1;

        const int dest = sum + move;
        if (dest >= 0 && dest < range) {
	  
          best_prev = forward_light(dest,v-1);
        }

        const double hyp = best_prev - dual_ptr[v][l];
        if (hyp < best_val)
          best_val = hyp;
      }
      forward_light(sum,v) = best_val;
    }
  }

  double min_val = 1e300;

  for (int k=rhs_lower_; k <= rhs_upper_; k++) {

    double min_msg = forward_light(zero_offset + k,nVars-1); 
    
    if (min_msg < min_val)
      min_val = min_msg;
  }

  return min_val;
}

/*virtual*/ double BILPConstraintDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  //std::cerr << "cm" << std::endl;

  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  int nPositive = 0;
  int nNegative = 0;

  for (uint k=0; k < nVars; k++) {
    if (positive_[k])
      nPositive++;
    else
      nNegative++;
  }

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  //some intermediate values can never result in a satisfied constraint 
  // => prune the range to the useful values
  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  /**** forward ****/

  Math3D::NamedTensor<double> forward(range,2,participating_var_.size(),1e100,MAKENAME(forward));
  Math3D::NamedTensor<uint> trace(range,2,participating_var_.size(),MAX_UINT,MAKENAME(trace));

  //init
  forward(zero_offset,0,0) = - dual_ptr[0][0];
  int init_mul = (positive_[0]) ? 1 : -1;
  forward(zero_offset+init_mul,1,0) = - dual_ptr[0][1];

  //proceed
  for (uint v=1; v < nVars; v++) {

    for (int sum=0; sum < range; sum++) {

      for (int l=0; l < 2; l++) {
	
        double best_prev = 1e75;
	
        int move = l;
        if (positive_[v]) //since we are tracing backward here
          move *= -1;

        const int dest = sum + move;
        if (dest >= 0 && dest < range) {
          for (int l_prev = 0; l_prev < 2; l_prev ++) {
            const double hyp = forward(dest,l_prev,v-1);
            if (hyp < best_prev) {
              best_prev = hyp;
              trace(sum,l,v) = l_prev;
            }
          }
        }

        forward(sum,l,v) = best_prev - dual_ptr[v][l];
      }
    }
  }

  double min_val = 1e300;

  uint best_l = MAX_UINT;
  uint best_k = MAX_UINT;

  for (uint l=0; l < 2; l++) {

    for (int k=rhs_lower_; k <= rhs_upper_; k++) {

      double min_msg = forward(zero_offset + k,l,nVars-1); 
    
      if (min_msg < min_val) {
        min_val = min_msg;
        best_k = k + zero_offset;
        best_l = l;
      }
    }
  }

  //std::cerr << "best_k: " << best_k << std::endl;
  //std::cerr << "best_l: " << best_l << std::endl;

  // traceback
  min_labels[nVars-1] = best_l;

  for (int v=nVars-2; v >= 0; v--) {

    uint prev_k = best_k;

    if (positive_[v+1])
      best_k -= best_l;
    else
      best_k += best_l;

    best_l = trace(prev_k,best_l,v+1);

    assert(best_l != MAX_UINT);

    min_labels[v] = best_l;
  }

  //std::cerr << "finished" << std::endl;

  return min_val;
}

/**********************************/

FactorDualOpt::FactorDualOpt(uint nVars, uint nFactors) : nUsedVars_(0), nUsedFactors_(0), optimize_called_(false) {

  var_.resize(nVars,0);
  factor_.resize(nFactors,0);

  first_shared_var_ = MAX_UINT;
  first_shared_factor_ = MAX_UINT;
}

//delete all allocated owned var. and factor nodes
FactorDualOpt::~FactorDualOpt() {
  for (uint i=0; i < std::min<uint>(nUsedVars_,first_shared_var_); i++)
    delete var_[i];

  for (uint i=0; i < std::min<uint>(nUsedFactors_,first_shared_factor_); i++)
    delete factor_[i];
}

/***** Setting up a factor graph ******/

uint FactorDualOpt::add_factor(DualFactorNode* node, bool owned) {
    
  if (owned && first_shared_factor_ != MAX_UINT) {
    INTERNAL_ERROR << " cannot create owned factor nodes after external factor nodes have been passed in";
    exit(1);
  }

  if (!owned && first_shared_factor_ == MAX_UINT)
    first_shared_factor_ = nUsedFactors_;
  
  nUsedFactors_++;

  uint prev_size = factor_.size();
  if (nUsedFactors_ > prev_size)
    factor_.resize(size_t(1.2*prev_size)+1,0);

  factor_[nUsedFactors_-1] = node;

  return nUsedFactors_-1;
}

//return var. number
uint FactorDualOpt::add_var(const Math1D::Vector<float>& unary_cost) {

  if (first_shared_var_ != MAX_UINT) {
    INTERNAL_ERROR << " cannot create generic var-nodes after external var-nodes have been passed in";
    exit(1);
  }

  DualVariableNode* ptr = new DualVariableNode(unary_cost);

  uint nPrevVars = nUsedVars_;

  nUsedVars_++;

  uint prev_size = var_.size();
  if (nUsedVars_ > prev_size) {
    var_.resize(size_t(1.2*prev_size)+1,0);
  }
  var_[nPrevVars] = ptr;

  return nPrevVars;
}

//NOTE: after calling this routine, owned vars can no longer be created
uint FactorDualOpt::pass_in_var_node(DualVariableNode* var) {

  uint nPrevVars = nUsedVars_;
  
  nUsedVars_++;
  
  uint prev_size = var_.size();
  if (first_shared_var_ == MAX_UINT)
    first_shared_var_ = prev_size;
  if (nUsedVars_ > prev_size) {
    var_.resize(size_t(1.2*prev_size)+1,0);
  }
  var_[nPrevVars] = var;
  
  return prev_size;
}

uint FactorDualOpt::add_generic_factor(const Math1D::Vector<uint> var, VarDimStorage<float>& cost) {

  assert(var.size() == cost.nDims());

  Storage1D<DualVariableNode*> var_nodes(var.size());
  for (uint k=0; k < var.size(); k++)
    var_nodes[k] = var_[var[k]];

  GenericDualFactorNode* ptr = new GenericDualFactorNode(var_nodes,cost);
  return add_factor(ptr);
}

uint FactorDualOpt::add_generic_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref) {

  assert(var1 < var_.size());
  assert(var2 < var_.size());

  Storage1D<DualVariableNode*> var_nodes(2);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];

  if (ref) {
    BinaryDualRefFactorNode* ptr = new BinaryDualRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    BinaryDualFactorNode* ptr = new BinaryDualFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}
  
uint FactorDualOpt::add_potts_factor(uint var1, uint var2, double lambda) {

  assert(var1 < var_.size());
  assert(var2 < var_.size());

  Storage1D<DualVariableNode*> var_nodes(2);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];

  PottsDualFactorNode* ptr = new PottsDualFactorNode(var_nodes, lambda);
  return add_factor(ptr);
}

//return factor number
uint FactorDualOpt::add_generic_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {


  assert(var1 < var_.size());
  assert(var2 < var_.size());
  assert(var3 < var_.size());

  Storage1D<DualVariableNode*> var_nodes(3);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];

  if (ref) {
    TernaryDualRefFactorNode* ptr = new TernaryDualRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    TernaryDualFactorNode* ptr = new TernaryDualFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}

uint FactorDualOpt::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  assert(var1 < var_.size());
  assert(var2 < var_.size());
  assert(var3 < var_.size());

  Storage1D<DualVariableNode*> var_nodes(3);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];

  SecondDiffDualFactorNode* ptr = new SecondDiffDualFactorNode(var_nodes,lambda);
  return add_factor(ptr);
}


//return factor number
uint FactorDualOpt::add_generic_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                                    const Storage1D<Math3D::Tensor<float> >& cost, bool ref) {

  assert(var1 < var_.size());
  assert(var2 < var_.size());
  assert(var3 < var_.size());
  assert(var4 < var_.size());

  Storage1D<DualVariableNode*> var_nodes(4);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];
  var_nodes[3] = var_[var4];

  if (ref) {
    FourthOrderDualRefFactorNode* ptr = new FourthOrderDualRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    FourthOrderDualFactorNode* ptr = new FourthOrderDualFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }

}


//all participating variables must be binary
uint FactorDualOpt::add_one_of_N_factor(const Math1D::Vector<uint>& var) {

  Storage1D<DualVariableNode*> var_nodes(var.size());

  for (uint k=0; k < var.size(); k++) {

    assert(var[k] < var_.size());
    var_nodes[k] = var_[var[k]];

    if (var_[var[k]]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of 1-of-N nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
  }

  OneOfNDualFactorNode* ptr = new OneOfNDualFactorNode(var_nodes);

  return add_factor(ptr);
}

//all participating variables must be binary
uint FactorDualOpt::add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost) {

  Storage1D<DualVariableNode*> var_nodes(var.size());

  for (uint k=0; k < var.size(); k++) {

    assert(var[k] < var_.size());
    var_nodes[k] = var_[var[k]];

    if (var_[var[k]]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of Cardinality nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
  }

  CardinalityDualFactorNode* ptr = new CardinalityDualFactorNode(var_nodes, card_cost);

  return add_factor(ptr);
}

//all participating variables must be binary
uint FactorDualOpt::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                                          int rhs_lower, int rhs_upper) {

  uint nUseful = 0;

  //check for variables whose value is essentially fixed due to the cost vector
  for (uint k=0; k < var.size(); k++) {
    
    const Math1D::Vector<float>& cur_cost = var_[var[k]]->cost();

    if (fabs(cur_cost[0] - cur_cost[1]) < 1e10)
      nUseful++;
    else {
      assert(cur_cost[0] < cur_cost[1]); 
      //otherwise need to adjust rhs_lower and upper (currently not implemented)
    }
  }

  if (nUseful != 0) {

    assert(nUseful >= 2);
    
    Storage1D<DualVariableNode*> vars(nUseful);
    Storage1D<bool> reduced_positive(nUseful);
  
    uint next = 0;
    
    for (uint k=0; k < var.size(); k++) {
    
      if (fabs(var_[var[k]]->cost()[0] - var_[var[k]]->cost()[1]) < 1e10) {
        vars[next] = var_[var[k]];
        reduced_positive[next] = positive[k];
        next++;
      }
    }

    assert(next == nUseful);
    assert(nUsedFactors_ < factor_.size());

    BILPConstraintDualFactorNode* ptr = new BILPConstraintDualFactorNode(vars, reduced_positive, rhs_lower, rhs_upper);

    return add_factor(ptr);
  }
  else {

    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;
    return MAX_UINT;
  }
}

//NOTE: after calling this routine, owned factors can no longer be created
uint FactorDualOpt::pass_in_factor_node(DualFactorNode* factor) {

  return add_factor(factor,false);
}

/**** run inference ***/

double FactorDualOpt::dual_bcd(uint nIter, DualBCDMode mode, bool init, bool quiet) {


  if (optimize_called_) {
    INTERNAL_ERROR << "you can only call dual_bcd() or subgradient_opt() once!. Exiting... " << std::endl;
    return 1e300;
  }
  optimize_called_ = true;

  labeling_.resize(nUsedVars_,0);

  uint arg_min;

  std::cerr.precision(10);

  if (init) {
    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->init_dual_vars();
    }
  }

  size_t effort_per_iter = 0;
  for (uint f=0; f < nUsedFactors_; f++) {
    effort_per_iter += factor_[f]->participating_nodes().size() * factor_[f]->participating_nodes().size();
  }

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "******** iteration " << iter << " ************" << std::endl;

    for (uint f=0; f < nUsedFactors_; f++) {
      factor_[f]->update_duals(mode);
    }
    
    double lower_bound = 0.0;

    if (!quiet && (iter % 1) == 0) {
      for (uint v=0; v < nUsedVars_; v++) {
        lower_bound += var_[v]->dual_value(arg_min);
        labeling_[v] = arg_min;
      }
      
      //std::cerr << "lb part for vars: " << lower_bound << std::endl;

      for (uint f=0; f < nUsedFactors_; f++)
        lower_bound += factor_[f]->dual_value();
      
      std::cerr << "lower bound: " << lower_bound << std::endl;
      
      std::cerr << "primal energy: " << labeling_energy() << std::endl;
    }
  }

  size_t message_effort = 0;
  for (uint f=0; f < nUsedFactors_; f++) {
    message_effort += factor_[f]->participating_nodes().size() * factor_[f]->participating_nodes().size();
  }
  message_effort *= nIter;
  std::cerr << "message efffort " << message_effort << std::endl;

  double lower_bound = 0.0;

  for (uint v=0; v < nUsedVars_; v++) {
    lower_bound += var_[v]->dual_value(arg_min);
    labeling_[v] = arg_min;
  }
  
  for (uint f=0; f < nUsedFactors_; f++)
    lower_bound += factor_[f]->dual_value();

  std::cerr << "lower bound: " << lower_bound << std::endl;

  return lower_bound;
}

double FactorDualOpt::subgradient_opt(uint nIter, double start_step_size) {

  if (optimize_called_) {
    INTERNAL_ERROR << "you can only call dual_bcd() or subgradient_opt() once!. Exiting... " << std::endl;
    return 1e300;
  }
  optimize_called_ = true;

  double best_dual = -1e300;

  Math1D::Vector<uint> var_label(nUsedVars_);

  Storage1D<Math1D::Vector<uint> > factor_label(nUsedFactors_);
  
  for (uint f=0; f < nUsedFactors_; f++) {
    factor_label[f].resize(factor_[f]->participating_nodes().size());
  }

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "iteration #" << iter << std::endl;

    double step_size = start_step_size / iter;

    double cur_bound = 0.0;

    for (uint v=0; v < nUsedVars_; v++) {
      uint cur_label;
      cur_bound += var_[v]->dual_value(cur_label);
      var_label[v] = cur_label;
    }
    
    std::cerr << "A" << std::endl;

    for (uint f=0; f < nUsedFactors_; f++) {
      cur_bound += factor_[f]->compute_minimizer(factor_label[f]);
    }

    if (cur_bound > best_dual)
      best_dual = cur_bound;

    std::cerr << "cur bound: " << cur_bound << ", best ever: " << best_dual << std::endl;

    for (uint v=0; v < nUsedVars_; v++) {
      const uint cur_label = var_label[v];

      const uint nCurFactors = var_[v]->neighboring_factor().size();
      const uint nLabels = var_[v]->nLabels();

      double* ptr = var_[v]->get_dual_var_start();

      for (uint k=0; k < nCurFactors; k++) {
        ptr[cur_label] += step_size;
        ptr += nLabels;
      }
    }

    for (uint f=0; f < nUsedFactors_; f++) {

      const Storage1D<DualVariableNode*>& cur_nodes = factor_[f]->participating_nodes();

      for (uint k=0; k < cur_nodes.size(); k++) {

        uint cur_label = factor_label[f][k];

        double* ptr = cur_nodes[k]->get_dual_vars(factor_[f]);
        ptr[cur_label] -= step_size;
      }
    }
  }

  labeling_ = var_label;

  return best_dual;
}


/**** get the state of the solver ****/

double FactorDualOpt::labeling_energy() {

  double energy = 0.0;

  std::map<DualVariableNode*,uint> label;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    energy += var_[k]->cost(labeling_[k]);
  }

  //std::cerr << "unary cost: " << energy << std::endl;

  //  std::cerr << "sum of labeling: " << labeling_.sum() << std::endl;
  //std::cerr << "nFactors: " << nUsedFactors_ << std::endl;

  for (uint k=0; k < nUsedFactors_; k++) {

    //std::cerr << "k: " << k << std::endl;

    const Storage1D<DualVariableNode*>& nodes = factor_[k]->participating_nodes();

    Math1D::Vector<uint> sub_labeling(nodes.size());
    for (uint i=0; i < nodes.size(); i++) {
      assert(label.find(nodes[i]) != label.end());
      sub_labeling[i] = label[nodes[i]];
    }

    energy += factor_[k]->cost(sub_labeling);
  }
  
  return energy;
}

const Math1D::Vector<uint>& FactorDualOpt::labeling() {

  return labeling_;
}


void FactorDualOpt::icm(uint nIter) {

  std::map<DualVariableNode*,uint> label;
  std::map<DualVariableNode*,uint> num;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    num[var_[k]] = k;
  }

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "ICM iter " << iter << std::endl;

    for (uint v=0; v < nUsedVars_; v++) {

      const Storage1D<DualFactorNode*>& factors = var_[v]->neighboring_factor();

      Storage1D<Math1D::Vector<uint> > sublabeling(factors.size());
      
      Math1D::Vector<uint> var_index(factors.size(),MAX_UINT);

      for (uint f=0; f < factors.size(); f++) {

        const Storage1D<DualVariableNode*>& factor_var = factors[f]->participating_nodes();

        sublabeling[f].resize(factor_var.size());

        for (uint k=0; k < factor_var.size(); k++) {

          sublabeling[f][k] = label[factor_var[k]];
          if (factor_var[k] == var_[v])
            var_index[f] = k;
        }
      }

      /****** now search for the conditional mode *****/
      double mode_val = 1e300;
      uint mode = MAX_UINT;

      for (uint l=0; l < var_[v]->nLabels(); l++) {

        double hyp = var_[v]->cost(l);
        for (uint f=0; f < factors.size(); f++) {
          sublabeling[f][var_index[f]] = l;
          hyp += factors[f]->cost(sublabeling[f]);
        }

        if (hyp < mode_val) {
          mode_val = hyp;
          mode = l;
        }
      }

      labeling_[v] = mode;
      label[var_[v]] = mode;
    }
  }

}
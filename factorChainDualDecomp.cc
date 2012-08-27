/******* written by Thomas Schoenemann as an employee of the University of Pisa, Italy, 2011 *****/

#include "factorChainDualDecomp.hh"

#include <map>
#include <vector>
#include <set>
#include "stl_out.hh"

//#define PRIMAL_DUAL_STEPSIZE

ChainDDVar::ChainDDVar(const Math1D::Vector<float>& cost) : cost_(cost) {}
  
void ChainDDVar::add_factor(ChainDDFactor* factor) {

  uint size = neighboring_factor_.size();

  neighboring_factor_.resize(size+1);
  neighboring_factor_[size] = factor;
}

void ChainDDVar::add_cost(const Math1D::Vector<float>& add_cost) {

  if (add_cost.size() != cost_.size()) {
    INTERNAL_ERROR << "cannot add cost due to incompatible vector sizes: " << cost_.size() << " and " << add_cost.size() << std::endl;
    exit(1);
  }

  float nC = nChains();

  for (uint i=0; i < cost_.size(); i++)
    cost_[i] += add_cost[i] / nC;
}

const Math1D::Vector<float>& ChainDDVar::cost() const {

  return cost_;
}

uint ChainDDVar::nLabels() const {

  return cost_.size();
}

const Storage1D<ChainDDFactor*>& ChainDDVar::neighboring_factor() const {

  return neighboring_factor_;
}

uint ChainDDVar::nChains() const {

  uint nChains = 0;

  for (uint k=0; k < neighboring_factor_.size(); k++) {
    if (neighboring_factor_[k]->prev_var() != this)
      nChains++;
  }

  return nChains;
}

void ChainDDVar::set_up_chains() {

  uint nChains = 0;

  for (uint k=0; k < neighboring_factor_.size(); k++) {
    if (neighboring_factor_[k]->prev_var() != this)
      nChains++;
  }

  cost_ *= 1.0 / nChains;
}

double ChainDDVar::dual_value(uint& arg_min) {

  Math1D::Vector<double> sum(cost_.size());
  for (uint l=0; l < cost_.size(); l++)
    sum[l] = cost_[l];

  for (uint f=0; f < neighboring_factor_.size(); f++) {

    sum += neighboring_factor_[f]->get_duals(this);
  }

  double best = 1e300;

  for (uint l=0; l < cost_.size(); l++) {

    if (sum[l] < best) {
      best = sum[l];
      arg_min = l;
    }
  }

  return best;
}

/********************************************/

ChainDDFactor::ChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars) : 
  prev_var_(0), next_var_(0), prev_factor_(0), next_factor_(0), involved_var_(involved_vars) {

  dual_var_.resize(involved_var_.size());

  for (uint v=0; v < involved_var_.size(); v++) {
    involved_var_[v]->add_factor(this);
    dual_var_[v].resize(involved_var_[v]->nLabels(),0.0);
  }
}


Math1D::Vector<double>& ChainDDFactor::get_duals(ChainDDVar* var) {

  for (uint k=0; k < involved_var_.size(); k++) {

    if (involved_var_[k] == var)
      return dual_var_[k];
  }

  assert(false);
  return dual_var_[0];
}

const Storage1D<ChainDDVar*>& ChainDDFactor::involved_vars() {
  return involved_var_;
}

ChainDDVar* ChainDDFactor::prev_var() const {
  return prev_var_;
}

ChainDDVar* ChainDDFactor::next_var() const {
  return next_var_;
}

ChainDDFactor* ChainDDFactor::prev_factor() const {
  return prev_factor_;
}

ChainDDFactor* ChainDDFactor::next_factor() const {
  return next_factor_;
}

void ChainDDFactor::set_prev_var(ChainDDVar* var) {
  prev_var_ = var;
}

void ChainDDFactor::set_next_var(ChainDDVar* var) {
  next_var_ = var;
}

void ChainDDFactor::set_prev_factor(ChainDDFactor* factor) {
  prev_factor_ = factor;
}

void ChainDDFactor::set_next_factor(ChainDDFactor* factor) {
  next_factor_ = factor;
}

/********************************************/

BinaryChainDDFactorBase::BinaryChainDDFactorBase(const Storage1D<ChainDDVar*>& involved_vars)
  : ChainDDFactor(involved_vars) {

  assert(involved_vars.size() == 2);
}

double BinaryChainDDFactorBase::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                Math2D::Matrix<uint>& trace, const Math2D::Matrix<float>& cost) {

  assert(out_var != in_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }

  forward.resize(involved_var_[idx]->nLabels());
  trace.resize(2,forward.size());

  if (idx == 0) {

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      double best = 1e300;
      uint arg_best = MAX_UINT;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double hyp = cost(l1,l2) - param[1][l2];

        if (hyp < best) {
          best = hyp;
          arg_best = l2;
        }
      }

      forward[l1] = best - param[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = arg_best;
    }

  }
  else {
    assert(idx == 1);

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      double best = 1e300;
      uint arg_best = MAX_UINT;

      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double hyp = cost(l1,l2) - param[0][l1];

        if (hyp < best) {
          best = hyp;
          arg_best = l1;
        }
      }

      forward[l2] = best - param[1][l2];
      trace(0,l2) = arg_best;
      trace(1,l2) = l2;
    }

  }

  return 0.0; //currently not removing a constant offset
}

/******/

BinaryChainDDFactor::BinaryChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math2D::Matrix<float>& cost) :
  BinaryChainDDFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ 
double BinaryChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {
  return cost_(labeling[0],labeling[1]);
}


/*virtual */
double BinaryChainDDFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                            const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                            Math2D::Matrix<uint>& trace) {


  return BinaryChainDDFactorBase::compute_forward(in_var,out_var,prev_forward,forward,trace,cost_);
}

/*************/

BinaryChainDDRefFactor::BinaryChainDDRefFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math2D::Matrix<float>& cost) :
  BinaryChainDDFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ 
double BinaryChainDDRefFactor::cost(const Math1D::Vector<uint>& labeling) {
  return cost_(labeling[0],labeling[1]);
}


/*virtual */
double BinaryChainDDRefFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                               const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                               Math2D::Matrix<uint>& trace) {

  return BinaryChainDDFactorBase::compute_forward(in_var,out_var,prev_forward,forward,trace,cost_);
}

/********************************************/

TernaryChainDDFactorBase::TernaryChainDDFactorBase(const Storage1D<ChainDDVar*>& involved_vars)
  : ChainDDFactor(involved_vars) {

  assert(involved_vars.size() == 3);
}

double TernaryChainDDFactorBase::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                 const Math1D::Vector<double>& prev_forward, const Math3D::Tensor<float>& cost,
                                                 Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) {

  assert(out_var != in_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }

  if (idx == 0) {

    forward.resize(nLabels1);
    trace.resize(3,nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      uint argbest2 = MAX_UINT;
      uint argbest3 = MAX_UINT;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double inter1 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - inter1 - param[2][l3];
	  
          if (hyp < best) {
            best = hyp;
            argbest2 = l2;
            argbest3 = l3;
          }
        }
      }
      
      forward[l1] = best - param[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = argbest2;
      trace(2,l1) = argbest3;
    }
  }
  else if (idx == 1) {

    forward.resize(nLabels2);
    trace.resize(3,nLabels2);

    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest3 = MAX_UINT;

      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        double inter1 = param[0][l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - inter1 - param[2][l3];
	  
          if (hyp < best) {
            best = hyp;
            argbest1 = l1;
            argbest3 = l3;
          }
        }
      }
      
      forward[l2] = best - param[1][l2];
      trace(0,l2) = argbest1;
      trace(1,l2) = l2;
      trace(2,l2) = argbest3;
    }    
  }
  else {
    assert(out_var == involved_var_[2]);

    forward.resize(nLabels3);
    trace.resize(3,nLabels3);

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest2 = MAX_UINT;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        const double inter1 = param[0][l1];

        for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
          double hyp = cost(l1,l2,l3) - inter1 - param[1][l2];
	  
          if (hyp < best) {
            best = hyp;
            argbest1 = l1;
            argbest2 = l2;
          }
        }
      }
      
      forward[l3] = best - param[2][l3];
      trace(0,l3) = argbest1;
      trace(1,l3) = argbest2;
      trace(2,l3) = l3;
    }
  }

  return 0.0; //presently not subtracting an offset
}

/*********/

TernaryChainDDFactor::TernaryChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math3D::Tensor<float>& cost)
  : TernaryChainDDFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ double TernaryChainDDFactor::compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                                         const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                         Math2D::Matrix<uint>& trace) {

  return TernaryChainDDFactorBase::compute_forward(incoming,outgoing,prev_forward,cost_,forward,trace);
}

/*virtual*/ 
double TernaryChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {
  return cost_(labeling[0],labeling[1],labeling[2]);
}

/*********/

TernaryChainDDRefFactor::TernaryChainDDRefFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math3D::Tensor<float>& cost)
  : TernaryChainDDFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ double TernaryChainDDRefFactor::compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                                            const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                            Math2D::Matrix<uint>& trace) {

  return TernaryChainDDFactorBase::compute_forward(incoming,outgoing,prev_forward,cost_,forward,trace);
}

/*virtual*/ 
double TernaryChainDDRefFactor::cost(const Math1D::Vector<uint>& labeling) {
  return cost_(labeling[0],labeling[1],labeling[2]);
}

SecondDiffChainDDFactor::SecondDiffChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, float lambda)
  : ChainDDFactor(involved_vars), lambda_(lambda) {

  assert(involved_vars.size() == 3);
}

/*virtual*/ 
double SecondDiffChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {
  
  int diff1 = int(labeling[1]) - int(labeling[0]);
  int diff2 = int(labeling[2]) - int(labeling[1]);
  
  int so_diff = diff2 - diff1;
  
  if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
    return 0.0; //no cost
  else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
    return lambda_;

  return 3*lambda_;
}


/*virtual*/
double SecondDiffChainDDFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                Math2D::Matrix<uint>& trace) {

  assert(out_var != in_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }

  if (idx == 0) {

    forward.resize(nLabels1);
    trace.resize(3,nLabels1);

    uint default_base2 = MAX_UINT;
    double best2 = 1e300;

    for (uint l=0; l < nLabels2; l++) {

      if (-param[1][l] < best2) {
        best2 = -param[1][l];
        default_base2 = l;
      }
    }

    uint default_base3 = MAX_UINT;
    double best3 = 1e300;

    for (uint l=0; l < nLabels3; l++) {

      if (-param[2][l] < best3) {
        best3 = -param[2][l];
        default_base3 = l;
      }
    }

    double base_cost = best2 + best3 + 3*lambda_;

    for (int l1=0; l1 < int(nLabels1); l1++) {
      
      double best = base_cost;

      uint best2 = default_base2;
      uint best3 = default_base3;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[1][l2] - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[1][l2] - param[2][l3] + lambda_;
          }

          if (hyp < best) {
            best = hyp;
            best2 = l2;
            best3 = l3;
          }
        }
      }

      forward[l1] = best - param[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = best2;
      trace(2,l1) = best3;
    }
  }
  else if (idx == 1) {  

    forward.resize(nLabels2);
    trace.resize(3,nLabels2);
    

    uint default_base1 = MAX_UINT;
    double best1 = 1e300;

    for (uint l=0; l < nLabels1; l++) {

      if (-param[0][l] < best1) {
        best1 = -param[0][l];
        default_base1 = l;
      }
    }

    uint default_base3 = MAX_UINT;
    double best3 = 1e300;

    for (uint l=0; l < nLabels3; l++) {

      if (-param[2][l] < best3) {
        best3 = -param[2][l];
        default_base3 = l;
      }
    }

    double base_cost = best1 + best3 + 3*lambda_;

    for (int l2=0; l2 < int(nLabels2); l2++) {
      
      double best = base_cost;
      uint best1 = default_base1;
      uint best3 = default_base3;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[0][l1] - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[0][l1] - param[2][l3] + lambda_;
          }

          if (hyp < best) {
            best = hyp;
            best1 = l1;
            best3 = l3;
          }
        }
      }

      forward[l2] = best - param[1][l2];
      trace(0,l2) = best1;
      trace(1,l2) = l2;
      trace(2,l2) = best3;
    }    
  }
  else {

    forward.resize(nLabels3);
    trace.resize(3,nLabels3);
    
    uint default_base1 = MAX_UINT;
    double best1 = 1e300;

    for (uint l=0; l < nLabels1; l++) {

      if (-param[0][l] < best1) {
        best1 = -param[0][l];
        default_base1 = l;
      }
    }

    uint default_base2 = MAX_UINT;
    double best2 = 1e300;

    for (uint l=0; l < nLabels2; l++) {

      if (-param[1][l] < best2) {
        best2 = -param[1][l];
        default_base2 = l;
      }
    }

    double base_cost = best1 + best2 + 3*lambda_;

    for (int l3=0; l3 < int(nLabels3); l3++) {
      
      double best = base_cost;
      uint best1 = default_base1;
      uint best2 = default_base2;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {
        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[0][l1] - param[1][l2];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[0][l1] - param[1][l2] + lambda_;
          }

          if (hyp < best) {
            best = hyp;
            best1 = l1;
            best2 = l2;
          }
        }
      }

      forward[l3] = best - param[2][l3];
      trace(0,l3) = best1;
      trace(1,l3) = best2;
      trace(2,l3) = l3;
    }
  }

  return 0.0;
}


/********************************************/

FourthOrderChainDDFactorBase::FourthOrderChainDDFactorBase(const Storage1D<ChainDDVar*>& involved_vars) 
  : ChainDDFactor(involved_vars) {
  assert(involved_vars.size() == 4);
}

double FourthOrderChainDDFactorBase::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                     const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                     Math2D::Matrix<uint>& trace, const Storage1D<Math3D::Tensor<float> >& cost) {

  assert(out_var != in_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }


  if (idx == 0) {

    forward.resize(nLabels1);
    trace.resize(4,nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      uint argbest2 = MAX_UINT;
      uint argbest3 = MAX_UINT;
      uint argbest4 = MAX_UINT;

      const Math3D::Tensor<float>& cur_cost = cost[l1];

      for (uint l2 = 0; l2 < nLabels2; l2++) {

       const double inter1 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double inter2 = inter1 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) - inter2 - param[3][l4];
	  
            if (hyp < best) {
              best = hyp;
              argbest2 = l2;
              argbest3 = l3;
              argbest4 = l4;
            }
          }
        }
      }
      
      forward[l1] = best - param[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = argbest2;
      trace(2,l1) = argbest3;
      trace(3,l1) = argbest4;
    }
  }
  else if (idx == 1) {

    forward.resize(nLabels2);
    trace.resize(4,nLabels2);
    
    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest3 = MAX_UINT;
      uint argbest4 = MAX_UINT;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];
	
        const double inter1 = param[0][l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double inter2 = inter1 + param[2][l3];
          
          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) - inter2 - param[3][l4];
	  
            if (hyp < best) {
              best = hyp;
              argbest1 = l1;
              argbest3 = l3;
              argbest4 = l4;
            }
          }
        }
      }
      
      forward[l2] = best - param[1][l2];
      trace(0,l2) = argbest1;
      trace(1,l2) = l2;
      trace(2,l2) = argbest3;
      trace(3,l2) = argbest4;
    }    
  }
  else if (idx == 2) {
    
    forward.resize(nLabels2);
    trace.resize(4,nLabels2);

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest2 = MAX_UINT;
      uint argbest4 = MAX_UINT;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
        
        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double inter2 = inter1 + param[1][l2];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) - inter2 - param[3][l4];
	  
            if (hyp < best) {
              best = hyp;
              argbest1 = l1;
              argbest2 = l2;
              argbest4 = l4;
            }
          }
        }
      }
      
      forward[l3] = best - param[2][l3];
      trace(0,l3) = argbest1;
      trace(1,l3) = argbest2;
      trace(2,l3) = l3;
      trace(3,l3) = argbest4;
    }
  }
  else {
    
    assert(idx == 3);

    forward.resize(nLabels3);
    trace.resize(4,nLabels3);

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest2 = MAX_UINT;
      uint argbest3 = MAX_UINT;

      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = param[0][l1];

        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double inter2 = inter1 + param[1][l2];

          for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
            double hyp = cur_cost(l2,l3,l4) - inter2 - param[2][l3];
	  
            if (hyp < best) {
              best = hyp;
              argbest1 = l1;
              argbest2 = l2;
              argbest3 = l3;
            }
          }
        }
      }
      
      forward[l4] = best - param[3][l4];
      trace(0,l4) = argbest1;
      trace(1,l4) = argbest2;
      trace(2,l4) = argbest3;
      trace(3,l4) = l4;
    }    
  }
  
  return 0.0; //presently not subtracting an offset
}

/*****/

FourthOrderChainDDFactor::FourthOrderChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, 
                                                   const Storage1D<Math3D::Tensor<float> >& cost) 
  : FourthOrderChainDDFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ 
double FourthOrderChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {

  return cost_[labeling[0]](labeling[1],labeling[2],labeling[3]);
}

/*virtual*/ double FourthOrderChainDDFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                             const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                             Math2D::Matrix<uint>& trace) {


  return FourthOrderChainDDFactorBase::compute_forward(in_var,out_var,prev_forward,forward,trace,cost_);
}

FourthOrderChainDDRefFactor::FourthOrderChainDDRefFactor(const Storage1D<ChainDDVar*>& involved_vars, 
                                                         const Storage1D<Math3D::Tensor<float> >& cost) 
  : FourthOrderChainDDFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ 
double FourthOrderChainDDRefFactor::cost(const Math1D::Vector<uint>& labeling) {

  return cost_[labeling[0]](labeling[1],labeling[2],labeling[3]);
}

/*virtual*/ double FourthOrderChainDDRefFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                                const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                                Math2D::Matrix<uint>& trace) {

  return FourthOrderChainDDFactorBase::compute_forward(in_var,out_var,prev_forward,forward,trace,cost_);
}

/****/

OneOfNChainDDFactor::OneOfNChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars) : 
  ChainDDFactor(involved_vars) {
  
  for (uint v=0; v < involved_vars.size(); v++)
    assert(involved_vars[v]->nLabels() == 2);
}

/*virtual*/ 
double OneOfNChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {

  if (labeling.sum() == 1)
    return 0.0;
  return 1e30;
}

/*virtual*/ 
double OneOfNChainDDFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                            const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                            Math2D::Matrix<uint>& trace) {


  uint nVars = involved_var_.size();

  assert(out_var != in_var);

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }

  forward.resize(involved_var_[idx]->nLabels());
  trace.resize(nVars,involved_var_[idx]->nLabels());
  trace.set_constant(0);

  double best_gain = 1e300;
  uint argmin = MAX_UINT;

  double sum = 0.0;

  for (uint i=0; i < nVars; i++) {

    if (involved_var_[i] == out_var) {
      forward[0] = -param[i][0];
      forward[1] = -param[i][1]; 
      trace(i,1) = 1;
    }
    else {
      double hyp = -param[i][1] + param[i][0];

      if (hyp < best_gain) {
        best_gain = hyp;
        argmin = i;
      }
      
      sum -= param[i][0];
    }
  }

  trace(argmin,0) = 1;
  forward[0] += sum + param[argmin][0] - param[argmin][1];
  forward[1] += sum;

  return 0.0; //presently not subtracting an offset
}


CardinalityChainDDFactor::CardinalityChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, 
                                                   const Math1D::Vector<float>& cost) :
  ChainDDFactor(involved_vars), cost_(cost) {

  for (uint v=0; v < involved_vars.size(); v++)
    assert(involved_vars[v]->nLabels() == 2);
}

/*virtual*/ 
double CardinalityChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {

  return cost_[labeling.sum()];
}

/*virtual*/
double CardinalityChainDDFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                                 Math2D::Matrix<uint>& trace) {

  uint nVars = involved_var_.size();

  assert(nVars >= 2);

  assert(out_var != in_var);

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }

  forward.resize(involved_var_[idx]->nLabels());
  trace.resize(nVars,involved_var_[idx]->nLabels());
  trace.set_constant(0);
  trace(idx,1) = 1;
    
  Storage1D<std::pair<double,uint> > rel_msg(nVars-1);

  double offs = 0.0;

  uint next = 0;
  for (uint k=0; k < nVars; k++) {

    if (k != idx) {
      const double val0 = - param[k][0];
      const double val1 = - param[k][1];

      rel_msg[next] = std::make_pair(val1 - val0,k);
      offs += -param[k][0];

      next++;
    }
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars-1);

  forward.set_constant(1e300);
  
  int best_c0 = 0;
  int best_c1 = 0;

  double cum_sum = 0.0;

  for (uint c=0; c < nVars; c++) {

    double hyp0 = cum_sum + cost_[c];
    if (hyp0 < forward[0]) {
      forward[0] = hyp0;
      best_c0 = c;
    }
    
    double hyp1 = cum_sum + cost_[c+1];
    if (hyp1 < forward[1]) {
      forward[1] = hyp1;
      best_c1 = c;
    }

    if (c+1 < nVars) 
      cum_sum += rel_msg[c].first;
  }

  for (uint l=0; l < 2; l++) {
    forward[l] += offs - param[idx][l];
  }

  for (int c = 0; c < best_c0-1; c++)
    trace(rel_msg[c].second,0) = 1;

  for (int c = 0; c < best_c1-1; c++)
    trace(rel_msg[c].second,1) = 1;

  return 0.0; //presently not subtracting an offset
}


BILPChainDDFactor::BILPChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Storage1D<bool>& positive,
                                     int rhs_lower, int rhs_upper) : 
  ChainDDFactor(involved_vars), positive_(positive), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

  for (uint v=0; v < involved_vars.size(); v++)
    assert(involved_vars[v]->nLabels() == 2);

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

  range_ = range;
  zero_offset_ = zero_offset;
}

/*virtual*/ 
double BILPChainDDFactor::cost(const Math1D::Vector<uint>& labeling) {

  short sum = 0;
  for (uint v=0; v < involved_var_.size(); v++) {
    if (positive_[v])
      sum += labeling[v];
    else
      sum -= labeling[v];
  }

  if (sum >= rhs_lower_ && sum <= rhs_upper_)
    return 0.0;
  return 1e30;
}

/*virtual*/ 
double BILPChainDDFactor::compute_forward(const ChainDDVar* in_var, const ChainDDVar* out_var,
                                          const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward_msg, 
                                          Math2D::Matrix<uint>& trace) {

  //based on [Potetz & Lee CVIU 2007]

  uint nVars = involved_var_.size();

  assert(out_var != in_var);

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] == in_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_forward[l];
      }
    }
    else {
      if (involved_var_[k] == out_var) {
        idx = k;
      }

      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= involved_var_[k]->cost()[l];
      }
    }
  }

  forward_msg.resize(involved_var_[idx]->nLabels());
  trace.resize(nVars,involved_var_[idx]->nLabels());

  /**** forward ****/

  Math3D::Tensor<double> forward(range_,2,idx+1);
  Math2D::Matrix<double> forward_light(range_,idx+1);
  Math2D::Matrix<uchar> forward_light_trace(range_,idx+1,255);
  

  //init
  for (int sum=0; sum < range_; sum++) {

    forward_light(sum,0) = 1e100;
    for (int l=0; l < 2; l++) {
      forward(sum,l,0) = 1e100;
    }
  }

  forward(zero_offset_,0,0) = -param[0][0];
  forward_light(zero_offset_,0) = -param[0][0];
  forward_light_trace(zero_offset_,0) = 0;
  const int init_mul = (positive_[0]) ? 1 : -1;
  if (int(zero_offset_)+init_mul >= 0
      && int(zero_offset_)+init_mul < range_) {
    forward(zero_offset_+init_mul,1,0) = -param[0][1];
    forward_light(zero_offset_+init_mul,0) = -param[0][1];
    forward_light_trace(zero_offset_+init_mul,0) = 1;
  }

  //proceed
  for (uint v=1; v <= idx; v++) {

    for (int sum=0; sum < range_; sum++) {

      for (int l=0; l < 2; l++) {
	
        double best_prev = 1e75;

        int move = l;
        if (positive_[v]) //since we are tracing backward here
          move *= -1;

        const int dest = sum + move;
        if (dest >= 0 && dest < range_) {

          best_prev = forward_light(dest,v-1);
        }

        forward(sum,l,v) = best_prev - param[v][l];
      }
      if (forward(sum,0,v) <= forward(sum,1,v)) {
        forward_light(sum,v) = forward(sum,0,v);
        forward_light_trace(sum,v) = 0;
      }
      else {
        forward_light(sum,v) = forward(sum,1,v);
        forward_light_trace(sum,v) = 1;
      }
    }
  }

  Math2D::Matrix<double> backward_light(range_,nVars,1e100);
  Math2D::Matrix<uchar> backward_light_trace(range_,nVars,255);

  /**** backward ****/

  const uint last_var = nVars-1;

  //init
  for (int sum=0; sum < range_; sum++) 
    backward_light(sum,last_var) = 1e100;

  backward_light(zero_offset_,last_var) = -param[last_var][0];
  assert(!isinf(backward_light(zero_offset_,last_var)));
  backward_light_trace(zero_offset_,last_var) = 0;
  const int end_mul = (positive_[last_var]) ? 1 : -1;
  backward_light(zero_offset_ + end_mul,last_var) = -param[last_var][1];
  backward_light_trace(zero_offset_ + end_mul,last_var) = 1;

  if (isinf(backward_light(zero_offset_ + end_mul,last_var))) {
    std::cerr << "last var is in var: " << (involved_var_[last_var] == in_var) << std::endl;
    std::cerr << "dual var: " << dual_var_[last_var][1] << std::endl;
  }

  assert(!isinf(backward_light(zero_offset_ + end_mul,last_var)));
	 

  //proceed
  for (int v=last_var-1; v > int(idx); v--) {

    for (int sum=0; sum < range_; sum++) {
      
      double best_prev = 1e75;

      for (int l=0; l < 2; l++) {

      	int move = l;
      	if (positive_[v]) //since we are tracing backward here
      	  move *= -1;

      	const int dest = sum + move;
        double hyp = 1e75;
      	if (dest >= 0 && dest < range_) {
      	  hyp = backward_light(dest,v+1) - param[v][l];
          if (isinf(hyp)) {

            std::cerr << "bw: " << backward_light(dest,v+1) << std::endl;
            std::cerr << "param: " << param[v][l] << std::endl;
            std::cerr << "v: " << v << ", nVars: " << nVars << std::endl;
            std::cerr << "is in-var: " << (involved_var_[v] == in_var) << std::endl;
          }

          assert(!isinf(hyp));
      	}

        if (hyp < best_prev) {
          best_prev = hyp;
          backward_light_trace(sum,v) = l;
        }
      }

      backward_light(sum,v) = best_prev;
    }
  }


  for (uint l=0; l < 2; l++) {

    double min_msg = 1e300;
    int best_s = 0;
    int best_r = 0;
    
    for (int s=0; s < (int) range_; s++) {
      
      double hyp = forward(s,l,idx);

      assert(!isinf(hyp));

      uint br = 0;

      if (idx+1 < nVars) {

        double best_bwd = 1e300;

        const int diff = (s - zero_offset_);

        for (int r=rhs_lower_; r <= rhs_upper_; r++) {
          const int other = r + zero_offset_ - diff; 
	  
          if (other >= 0 && other < (int) range_) {

            if (backward_light(other,idx+1) < best_bwd) {
              best_bwd = backward_light(other,idx+1);
              br = other;
            }
          }
        }

        hyp += best_bwd;
        assert(!isinf(hyp));
      }
      else {
        if (s < int(rhs_lower_ + zero_offset_) || s > int(rhs_upper_ + zero_offset_)) 
          hyp = 1e300;
      }

      if (hyp < min_msg) {
        min_msg = hyp;
        best_s = s;
        best_r = br;
      }
    }

    forward_msg[l] = min_msg;
    trace(idx,l) = l;

    //std::cerr << "trace back" << std::endl;

    //a) trace backward to start
    uint cur_l = l;
    uint cur_sum = best_s;
    for (int k=idx; k >= 0; k--) {

      assert(cur_l < 2);
      
      trace(k,l) = cur_l;
      if (positive_[k])
        cur_sum -= cur_l;
      else
        cur_sum += cur_l;

      if (k > 0)
        cur_l = forward_light_trace(cur_sum,k-1);
    }

    //std::cerr << "trace for" << std::endl;

    //a) trace forward to end
    if (idx + 1 < nVars) {
      cur_sum = best_r;
      cur_l = backward_light_trace(cur_sum,idx+1);
      for (uint k=idx+1; k < nVars; k++) {

        if (!(cur_l < 2)) {
          std::cerr << "k: " << k << ", nVars: " << nVars << ", cur_l: " << cur_l << std::endl;
          std::cerr << "associated cost: " << forward_msg[l] << std::endl;
        }

        assert(cur_l < 2);

        trace(k,l) = cur_l;
        if (positive_[k])
          cur_sum -= cur_l;
        else
          cur_sum += cur_l;
	
        if (k+1 < nVars)
          cur_l = backward_light_trace(cur_sum,k+1);
      }
    }
  }

  return 0.0; //currently not subtracting an offest
}

/********************************************/

FactorChainDualDecomposition::FactorChainDualDecomposition(uint nVars, uint nFactors) : 
  nUsedVars_(0), nUsedFactors_(0), optimize_called_(false) {
  var_.resize(nVars);
  factor_.resize(nFactors);
}

FactorChainDualDecomposition::~FactorChainDualDecomposition() {

  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];
  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

void FactorChainDualDecomposition::add_var(const Math1D::Vector<float>& cost) {

  assert(!optimize_called_);

  if (nUsedVars_ == var_.size())
    var_.resize(uint(1.2*nUsedVars_)+4);

  assert(nUsedVars_ < var_.size());
  var_[nUsedVars_] = new ChainDDVar(cost);

  nUsedVars_++;
}

void FactorChainDualDecomposition::add_factor(ChainDDFactor* fac) {

  assert(!optimize_called_);

  if (nUsedFactors_ == factor_.size())
    factor_.resize(uint(1.2*nUsedFactors_)+4);

  factor_[nUsedFactors_] = fac;
  nUsedFactors_++;
}

void FactorChainDualDecomposition::add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref) {

  assert(var1 < nUsedVars_);
  assert(var2 < nUsedVars_);

  Storage1D<ChainDDVar*> vars(2);
  vars[0] = var_[var1];
  vars[1] = var_[var2];

  ChainDDFactor* newFac;

  if (ref)
    newFac = new BinaryChainDDRefFactor(vars,cost);
  else
    newFac = new BinaryChainDDFactor(vars,cost);

  add_factor(newFac);
}

void FactorChainDualDecomposition::add_ternary_factor(uint var1, uint var2, uint var3, 
                                                      const Math3D::Tensor<float>& cost, bool ref) {

  Storage1D<ChainDDVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  ChainDDFactor* new_fac = 0;

  if (!ref)
    new_fac = new TernaryChainDDFactor(vars,cost);
  else
    new_fac = new TernaryChainDDRefFactor(vars,cost);

  add_factor(new_fac);
}

void FactorChainDualDecomposition::add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                                           const Storage1D<Math3D::Tensor<float> >& cost, bool ref) {

  Storage1D<ChainDDVar*> vars(4);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];
  vars[3] = var_[var4];

  ChainDDFactor* new_fac;

  if (ref)
    new_fac = new FourthOrderChainDDRefFactor(vars,cost);
  else
    new_fac = new FourthOrderChainDDFactor(vars,cost);

  add_factor(new_fac);
}

void FactorChainDualDecomposition::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  Storage1D<ChainDDVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  add_factor(new SecondDiffChainDDFactor(vars,lambda));
}

void FactorChainDualDecomposition::add_one_of_n_factor(const Math1D::Vector<uint>& var) {

  Storage1D<ChainDDVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  add_factor(new OneOfNChainDDFactor(vars));
}

void FactorChainDualDecomposition::add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost) {

  if (var.size() == 1) {
    var_[var[0]]->add_cost(cost);
  }
  else {

    Storage1D<ChainDDVar*> vars(var.size());
  
    for (uint k=0; k < var.size(); k++)
      vars[k] = var_[var[k]];

    add_factor(new CardinalityChainDDFactor(vars,cost));
  }
}

void FactorChainDualDecomposition::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                                                         int rhs_lower, int rhs_upper) {

  uint nUseful = 0;
  for (uint k=0; k < var.size(); k++) {
    
    if (fabs(var_[var[k]]->cost()[0] - var_[var[k]]->cost()[1]) < 1e10)
      nUseful++;
    else {
      assert(var_[var[k]]->cost()[0] < var_[var[k]]->cost()[1]); 
      //otherwise need to adjust rhs_lower and upper (currently not implemented)
    }
  }

  if (nUseful != 0) {
    // if (nUseful < 2) {
    //   std::cerr << "only " << nUseful << " out of " << var.size() << " variables are actually not fixed" << std::endl;
    
    //   for (uint k=0; k < var.size(); k++)
    //     std::cerr << "cost: " << var_[var[k]]->cost() << std::endl;

    //   std::cerr << "var: " << var << std::endl;
    // }

    assert(nUseful >= 2);
    
    Storage1D<ChainDDVar*> vars(nUseful);
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

    add_factor(new BILPChainDDFactor(vars,reduced_positive,rhs_lower,rhs_upper));
  }
  else {
    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;
  }
}

const Math1D::Vector<uint>& FactorChainDualDecomposition::labeling() {
  return labeling_;
}

void FactorChainDualDecomposition::set_up_chains() {

  //use non-monotonic chains

  uint nChains = 0;
  uint nAtLeast5 = 0;
  uint nAtLeast10 = 0;
  uint nAtLeast25 = 0;

  for (uint f=0; f < nUsedFactors_; f++) {
    
    if (factor_[f]->prev_var() == 0 && factor_[f]->next_var() == 0) {

      uint length = 1;

      std::set<ChainDDVar*> current_vars;

      //std::cerr << "extend lower" << std::endl;

      bool extension_found = true;

      ChainDDFactor* cur_factor = factor_[f];

#if 1
      //extend lower end
      while (extension_found) {

        extension_found = false;

        const Storage1D<ChainDDVar*>& involved_vars = cur_factor->involved_vars();

        for (uint k=0; k < involved_vars.size(); k++) {
          current_vars.insert(involved_vars[k]);
        }
	
        for (double k=0; k < involved_vars.size(); k++) {
	
          ChainDDVar* var = involved_vars[k];

          if (var == cur_factor->next_var())
            continue;

          const Storage1D<ChainDDFactor*>& adjacent_factor = var->neighboring_factor();

          for (uint l=0; l < adjacent_factor.size(); l++) {
	    
            ChainDDFactor* hyp_factor = adjacent_factor[l];
            const Storage1D<ChainDDVar*>& hyp_involved_vars = hyp_factor->involved_vars();

            bool is_valid_extension = false;

            if (hyp_factor != cur_factor
                && hyp_factor->prev_var() == 0 && hyp_factor->next_var() == 0) {
	      
              is_valid_extension = true;
	      
              for (uint v=0; v < hyp_involved_vars.size(); v++) {
		
                if (hyp_involved_vars[v] != var && current_vars.find(hyp_involved_vars[v]) != current_vars.end())
                  is_valid_extension = false;
              }
	      
              if (is_valid_extension) {

                extension_found = true;
                cur_factor->set_prev_var(var);
                cur_factor->set_prev_factor(hyp_factor);

                hyp_factor->set_next_var(var);
                hyp_factor->set_next_factor(cur_factor);
		
                cur_factor = hyp_factor;

                length++;
		
                break;
              }
            }
          }
	  
          if (extension_found)
            break;
        }
      }
#endif

#if 1
      //extend upper end

      //std::cerr << "extend upper" << std::endl;

      cur_factor = factor_[f];
      extension_found = true;

      while (extension_found) {

        extension_found = false;

        const Storage1D<ChainDDVar*>& involved_vars = cur_factor->involved_vars();

        for (uint k=0; k < involved_vars.size(); k++) {
          current_vars.insert(involved_vars[k]);
        }

        for (double k=0; k < involved_vars.size(); k++) {
	
          ChainDDVar* var = involved_vars[k];

          if (var == cur_factor->prev_var())
            continue;
	  
          const Storage1D<ChainDDFactor*>& adjacent_factor = var->neighboring_factor();
	  
          for (uint l=0; l < adjacent_factor.size(); l++) {
	    
            ChainDDFactor* hyp_factor = adjacent_factor[l];

            bool is_valid_extension = false;

            if (hyp_factor != cur_factor
                && hyp_factor->prev_var() == 0 && hyp_factor->next_var() == 0) {

              const Storage1D<ChainDDVar*>& hyp_involved_vars = hyp_factor->involved_vars();
	      
              is_valid_extension = true;
	      
              for (uint v=0; v < hyp_involved_vars.size(); v++) {
		
                if (hyp_involved_vars[v] != var && current_vars.find(hyp_involved_vars[v]) != current_vars.end())
                  is_valid_extension = false;
              }
	      
              if (is_valid_extension) {

                extension_found = true;
                cur_factor->set_next_var(var);
                cur_factor->set_next_factor(hyp_factor);

                hyp_factor->set_prev_var(var);
                hyp_factor->set_prev_factor(cur_factor);
		
                cur_factor = hyp_factor;

                length++;

                break;
              }
            }	    

          }
          if (extension_found)
            break;
        }
      }
#endif
      
      nChains++;

      if (length >= 5)
        nAtLeast5++;
      if (length >= 10)
        nAtLeast10++;
      if (length >= 25)
        nAtLeast25++;


      // if (length > 15)
      //   std::cerr << "chain length " << length << std::endl;
    }
  }

  std::cerr << nAtLeast5 << " chains have length at least 5." << std::endl;
  std::cerr << nAtLeast10 << " chains have length at least 10." << std::endl;
  std::cerr << nAtLeast25 << " chains have length at least 25." << std::endl;
  std::cerr << nChains << " chains in total, " << nUsedFactors_ << " factors" << std::endl;
}

double FactorChainDualDecomposition::optimize(uint nIter, double start_step_size) {


  std::cerr.precision(10);

  std::cerr << "subgradient optimization" << std::endl;
  std::cerr << nUsedFactors_ << " factors" << std::endl;

  if (!optimize_called_) {
    set_up_chains();

    for (uint v=0; v < nUsedVars_; v++)
      var_[v]->set_up_chains();
  }

  optimize_called_ = true;

  bool projective = true;

  double best_dual = -1e300;
#ifdef PRIMAL_DUAL_STEPSIZE
  double best_primal = 1e300;
#endif

  Math1D::Vector<uint> var_label(nUsedVars_);

  labeling_.resize(nUsedVars_,0);

  Storage1D<Math1D::Vector<uint> > factor_label(nUsedFactors_);

  // double denom = 0.0;
  
  for (uint f=0; f < nUsedFactors_; f++) {
    factor_label[f].resize(factor_[f]->involved_vars().size());
    //   denom += factor_[f]->involved_vars().size();
  }

  std::map<ChainDDVar*,uint> var_num;
  for (uint v=0; v < nUsedVars_; v++)
    var_num[var_[v]] = v;

  std::map<ChainDDFactor*,uint> factor_num;
  for (uint f=0; f < nUsedFactors_; f++)
    factor_num[factor_[f]] = f;

  uint nIncreases = 1;

  double delta = start_step_size;

  size_t effort_per_iter = 0;

  for (uint f=0; f < nUsedFactors_; f++) {

    uint size = factor_[f]->involved_vars().size();

    effort_per_iter += size;
  }

  //store chains for efficient access

  std::vector<std::vector<ChainDDFactor*> > chain;
  std::vector<std::vector<ChainDDVar*> > out_var;
  std::vector<ChainDDVar*> in_var;

  for (uint f=0; f < nUsedFactors_; f++) {

    if (factor_[f]->prev_factor() == 0) {

      chain.push_back(std::vector<ChainDDFactor*>());
      out_var.push_back(std::vector<ChainDDVar*>());
	
      ChainDDFactor* cur_factor = factor_[f];

      ChainDDVar* cur_in_var = 0;
      for (uint k=0; k < cur_factor->involved_vars().size(); k++) {
        if (cur_factor->involved_vars()[k] != cur_factor->next_var()) {
          cur_in_var = cur_factor->involved_vars()[k];
          break;
        }
      }
      in_var.push_back(cur_in_var);
      
      while (cur_factor != 0) {
        chain.back().push_back(cur_factor);
        assert(cur_factor->next_var() != in_var.back());
        //assert(std::find(out_var.begin(),out_var.end(),cur_factor->next_var()) == out_var.end());
        out_var.back().push_back(cur_factor->next_var());
        cur_factor = cur_factor->next_factor();
      }

      //find chain end
      for (uint k=0; k < chain.back().back()->involved_vars().size(); k++) {
      
        if (chain.back().back()->involved_vars()[k] != chain.back().back()->prev_var() 
            && chain.back().back()->involved_vars()[k] != in_var.back()) { //can happen for chains of length 1
          out_var.back().back() = chain.back().back()->involved_vars()[k];
          break;
        }
      }      
    }
  }


  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "iteration #" << iter << std::endl;

    uint nDisagreements = 0;

    //double step_size = start_step_size / iter;
    double step_size = start_step_size / nIncreases;

    double cur_bound = 0.0;

    if (!projective) {
      for (uint v=0; v < nUsedVars_; v++) {
        uint cur_label = 0;
        cur_bound += var_[v]->dual_value(cur_label);
        var_label[v] = cur_label;
      }
    
      std::cerr << "A, intermediate bound: " << cur_bound << std::endl;
    }

    //uint nLongChainsProcessed = 0;


    for (uint k=0; k < chain.size(); k++) {

      if (true) {

        const std::vector<ChainDDFactor*>& cur_chain = chain[k];
        const std::vector<ChainDDVar*>& cur_out_var = out_var[k];
        ChainDDVar* cur_in_var = in_var[k];
        
        uint chain_length = cur_chain.size();

        Math1D::NamedVector<double> forward1(MAKENAME(forward1));
        Math1D::NamedVector<double> forward2(cur_in_var->nLabels(),MAKENAME(forward2));
        for (uint l=0; l < forward2.size(); l++)
          forward2[l] = cur_in_var->cost()[l];

        NamedStorage1D<Math2D::Matrix<uint> > trace(chain_length,MAKENAME(trace));

        //std::cerr << "start forward" << std::endl;

        //std::cerr << "chain: " << chain << std::endl;
        //std::cerr << "chain length: " << chain_length << std::endl;

        //compute forward
        cur_bound += cur_chain[0]->compute_forward(cur_in_var,cur_out_var[0],forward2,forward1,trace[0]);

        for (uint k=1; k < chain_length; k++) {

          Math1D::Vector<double>& last_forward = ((k % 2) == 1) ? forward1 : forward2;
          Math1D::Vector<double>& new_forward = ((k % 2) == 0) ? forward1 : forward2;

          cur_bound += cur_chain[k]->compute_forward(cur_out_var[k-1],cur_out_var[k],last_forward,new_forward,trace[k]);
        }

        //std::cerr << "start traceback" << std::endl;

        //traceback
        Math1D::Vector<double>& total_forward = ((chain_length-1) % 2 == 0) ? forward1 : forward2;

        assert(total_forward.size() == cur_out_var[chain_length-1]->nLabels());

        double best = 1e300;
        uint arg_best = MAX_UINT;
        for (uint l=0; l < total_forward.size(); l++) {
          if (total_forward[l] < best) {

            best = total_forward[l];
            arg_best = l;
          }
        }

        assert(arg_best < MAX_UINT);

        cur_bound += best;

        for (int k=chain_length-1; k >= 0; k--) {

          //std::cerr << "k: " << k << std::endl;

          Math1D::Vector<uint>& labeling = factor_label[factor_num[cur_chain[k]]];	  
#ifdef PRIMAL_DUAL_STEPSIZE
          const Storage1D<ChainDDVar*>& involved_vars = factor_[factor_num[cur_chain[k]]]->involved_vars();
#endif          

          assert(labeling.size() == trace[k].xDim());

          for (uint v=0; v < trace[k].xDim(); v++) {
            labeling[v] = trace[k](v,arg_best);
#ifdef PRIMAL_DUAL_STEPSIZE
            labeling_[var_num[involved_vars[v]]] = labeling[v];
#endif
          }

          //update arg_best
          if (k > 0) {
            for (uint v=0; v < trace[k].xDim(); v++) {
	      
              if (cur_chain[k]->involved_vars()[v] == cur_out_var[k-1]) {
                arg_best = labeling[v];
                break;
              }
            }
          }
        }
      }
    }

    if (cur_bound > best_dual) {
      best_dual = cur_bound;
      delta *= 1.5;
    }
    else {
      nIncreases++;
      delta *= 0.95;
    }

    std::cerr << "cur bound: " << cur_bound << ", best ever: " << best_dual << std::endl;

    //TRIAL: primal-dual bound
#ifdef PRIMAL_DUAL_STEPSIZE
    double cur_primal = 0.0;
    for (uint v=0; v < nUsedVars_; v++)
      cur_primal += var_[v]->cost()[labeling_[v]] * var_[v]->nChains(); //the cost are scaled internally

    for (uint f=0; f < nUsedFactors_; f++) {

      const Storage1D<ChainDDVar*>& involved_vars = factor_[f]->involved_vars();

      Math1D::Vector<uint> labeling(involved_vars.size());
      for (uint v=0; v < involved_vars.size(); v++)
        labeling[v] = labeling_[var_num[involved_vars[v]]];
      
      cur_primal += factor_[f]->cost(labeling);
    }

    best_primal = std::min(best_primal,cur_primal);
    std::cerr << "cur primal: " << cur_primal << ", best primal: " << best_primal << std::endl;

    assert(best_primal >= best_dual - 1e-5);

    Storage1D<Math1D::Vector<double> > gradient(nUsedVars_);
    for (uint v=0; v < nUsedVars_; v++)
      gradient[v].resize(var_[v]->cost().size(),0.0);

    for (uint f=0; f < nUsedFactors_; f++) {
      
      for (uint k=0; k < factor_label[f].size(); k++) {
        
        ChainDDVar* cur_var = factor_[f]->involved_vars()[k];
	
        uint cur_fac_label = factor_label[f][k];
        gradient[var_num[cur_var]][cur_fac_label] += 1.0;
      }
    }
    double grad_sqr_norm = 0.0;
    for (uint v=0; v < nUsedVars_; v++)
      grad_sqr_norm += gradient[v].sqr_norm();

    step_size = start_step_size * (best_primal - cur_bound) / std::max(1.0,grad_sqr_norm);
#endif
    //END_TRIAL


    //TRIAL
    //like in the Torressani et al. paper
    //step_size = (best_dual + delta - cur_bound) / denom;
    //step_size = std::min(step_size,start_step_size);
    //std::cerr << "step size: " << step_size << std::endl;
    //END_TRIAL

    if (!projective) {

      for (uint f=0; f < nUsedFactors_; f++) {
	
        for (uint k=0; k < factor_label[f].size(); k++) {
	  
          ChainDDVar* cur_var = factor_[f]->involved_vars()[k];
	  
          uint cur_fac_label = factor_label[f][k];
          uint cur_var_label = var_label[var_num[cur_var]];
	  
          if (cur_fac_label != cur_var_label) {
	    
            nDisagreements++;

            factor_[f]->get_duals(cur_var)[cur_var_label] += step_size;
            factor_[f]->get_duals(cur_var)[cur_fac_label] -= step_size;
          }
        }
      }
      
      std::cerr << nDisagreements << " disagreements" << std::endl;
      //std::cerr << "total repar cost: " << total_repar << std::endl;
    }
    else {

      for (uint f=0; f < nUsedFactors_; f++) {

        for (uint k=0; k < factor_label[f].size(); k++) {
	  
          ChainDDVar* cur_var = factor_[f]->involved_vars()[k];
	  
          uint cur_fac_label = factor_label[f][k];
          factor_[f]->get_duals(cur_var)[cur_fac_label] -= step_size;
          assert(!isinf(factor_[f]->get_duals(cur_var)[cur_fac_label]));
        }
      }

      //re-project
      for (uint v=0; v < nUsedVars_; v++) {

        if (var_[v]->neighboring_factor().size() > 0) {

          Math1D::Vector<double> sum(var_[v]->nLabels());
	  
          for (uint k=0; k < var_[v]->neighboring_factor().size(); k++) {
            sum += var_[v]->neighboring_factor()[k]->get_duals(var_[v]);
          }
	  
          sum *= 1.0 / var_[v]->neighboring_factor().size();
	  
          for (uint k=0; k < var_[v]->neighboring_factor().size(); k++) {
	    
            var_[v]->neighboring_factor()[k]->get_duals(var_[v]) -= sum;
            for (uint l=0; l < var_[v]->nLabels(); l++)
              assert(!isinf(var_[v]->neighboring_factor()[k]->get_duals(var_[v])[l]));
          }
        }
      }
    }
  }

  if (projective) {

    labeling_.set_constant(MAX_UINT);
    
    for (uint f=0; f < nUsedFactors_; f++) {
      
      for (uint k=0; k < factor_label[f].size(); k++) {
	
        ChainDDVar* cur_var = factor_[f]->involved_vars()[k];
        
        uint vnum = var_num[cur_var];

        if (labeling_[vnum] == MAX_UINT)
          labeling_[vnum] = factor_label[f][k];
      }
    }
  }
  else
    labeling_ = var_label;

  size_t message_effort = 0;

  for (uint f=0; f < nUsedFactors_; f++) {

    uint size = factor_[f]->involved_vars().size();

    message_effort += size;
  }
  
  message_effort *= nIter;
  std::cerr << "message effort: " << message_effort << std::endl;

  return best_dual;
}



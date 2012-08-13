/**** written by Thomas Schoenemann as a private person without employment, August 2011 *****/


#include "factorTRWS.hh"

#include <set>
#include "stl_out.hh"


CumTRWSVar::CumTRWSVar(const Math1D::Vector<float>& cost, uint rank) : cost_(cost), rank_(rank), nChains_(1) {
  cum_cost_.resize(cost.size());
  for (uint k=0; k < cum_cost_.size(); k++) {
    //convert float->double
    cum_cost_[k] = cost[k];
  }
}

void CumTRWSVar::add_cost(const Math1D::Vector<float>& add_cost) {

  if (add_cost.size() != cost_.size()) {
    INTERNAL_ERROR << "cannot add cost due to incompatible vector sizes: " << cost_.size() << " and " << add_cost.size() << std::endl;
    exit(1);
  }

  for (uint i=0; i < cost_.size(); i++)
    cost_[i] += add_cost[i] / float(nChains_);
}

void CumTRWSVar::set_up_chains() {

  nChains_ = 0;

  {
    uint nOutgoing = 0;
    uint nMiddle = 0;
    uint nIncoming = 0;

    for (uint k=0; k < adjacent_factor_.size(); k++) {

      if (adjacent_factor_[k]->min_rank() == rank_)
        nOutgoing++;
      else if (adjacent_factor_[k]->max_rank() == rank_)
        nIncoming++;
      else
        nMiddle++;
    }

    uint nChainlessChains = nMiddle + std::max(nOutgoing,nIncoming);

    nChains_ = nChainlessChains;
  }

  //initialize cost
  for (uint l=0; l < cost_.size(); l++)
    cum_cost_[l] = cost_[l] / nChains_;
}

void CumTRWSVar::add_factor(CumTRWSFactor* factor) {

  uint size = adjacent_factor_.size();
  
  adjacent_factor_.resize(size+1);
  adjacent_factor_[size] = factor;
}

double CumTRWSVar::average(uint& arg_min) {

  const uint nLabels = cost_.size();
  const uint nFactors = adjacent_factor_.size();

  for (uint l=0; l < nLabels; l++)
    cum_cost_[l] = cost_[l]; //convert float -> double

  for (uint k=0; k < nFactors; k++)
    cum_cost_ += adjacent_factor_[k]->reparameterization(this);

  double offs = 1e300;

  for (uint l=0; l < nLabels; l++) {
    if (cum_cost_[l] < offs) {
      offs = cum_cost_[l];
      arg_min = l;
    }
  }

  for (uint l=0; l < nLabels; l++) {
    cum_cost_[l] -= offs;
  }

  cum_cost_ *= 1.0 / nChains_;

  return offs;
}

uint CumTRWSVar::rank() const {
  return rank_;
}

void CumTRWSVar::set_rank(uint rank) {
  rank_ = rank;
}

size_t CumTRWSVar::nLabels() const {
  return cost_.size();
}

const Storage1D<CumTRWSFactor*>& CumTRWSVar::adjacent_factor() const {
  return adjacent_factor_;
}

const Math1D::Vector<double>& CumTRWSVar::cost() const {
  return cum_cost_;
}

/**************************/

CumTRWSFactor::CumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars) :
  involved_var_(involved_vars)  {

  uint nVars = involved_vars.size();
  assert(nVars >= 2);

  reparameterization_.resize(nVars);

  for (uint v=0; v < nVars; v++) {

    uint nLabels = involved_vars[v]->nLabels();
    reparameterization_[v].resize(nLabels,0.0);

    involved_vars[v]->add_factor(this);
  }

  compute_rank_range();
}
  
uint CumTRWSFactor::min_rank() const {
  return min_rank_;
}

uint CumTRWSFactor::max_rank() const {
  return max_rank_;
}

void CumTRWSFactor::compute_rank_range() {

  min_rank_ = MAX_UINT;
  max_rank_ = 0;

  const uint nVars = involved_var_.size();

  for (uint v=0; v < nVars; v++) {

    min_rank_ = std::min(min_rank_,involved_var_[v]->rank());
    max_rank_ = std::max(max_rank_,involved_var_[v]->rank());
  }
}

const Storage1D<CumTRWSVar*>& CumTRWSFactor::involved_vars() const {
  return involved_var_;
}

const Math1D::Vector<double>& CumTRWSFactor::reparameterization(const CumTRWSVar* var) const {

  for (uint k=0; k < involved_var_.size(); k++) {

    if (involved_var_[k] == var)
      return reparameterization_[k];
  }

  assert(false);
  return reparameterization_[0];
}


/********************/

BinaryCumTRWSFactor::BinaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                         const Math2D::Matrix<float>& cost) : 
  CumTRWSFactor(involved_vars), cost_(cost) {
  assert(involved_vars.size() == 2);
}

/*virtual*/ double BinaryCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {
  
  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  //this routine also updates reparameterization_
  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++) {
    param[k] -= involved_var_[k]->cost();
  }

  uint idx = 0;

  if (var == involved_var_[0]) {

    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double hyp = cost_(l1,l2) - param[1][l2];

        if (hyp < best)
          best = hyp;
      }

      reparameterization_[idx][l1] = best;
    }
  }
  else {
    assert(var == involved_var_[1]);

    idx = 1;

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      double best = 1e300;

      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double hyp = cost_(l1,l2) - param[0][l1];

        if (hyp < best)
          best = hyp;
      }

      reparameterization_[idx][l2] = best;
    }
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;
  
  return offs;
}

/********************/

TernaryCumTRWSFactorBase::TernaryCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars) :
  CumTRWSFactor(involved_vars) {
  assert(involved_vars.size() == 3);
}
  
double TernaryCumTRWSFactorBase::compute_reparameterization(CumTRWSVar* var, const Math3D::Tensor<float>& cost) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  //this routine also updates reparameterization_
  uint idx = 0;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->cost();
  
  if (var == involved_var_[0]) {
    
    idx = 0;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - w2 - param[2][l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l1] = best;
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - w1 - param[2][l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l2] = best;
    }    
  }
  else {
    assert(var == involved_var_[2]);
    
    idx = 2;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
          double hyp = cost(l1,l2,l3) - w1 - param[1][l2];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l3] = best;
    }
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;

  return offs;  
}

/********************/

TernaryCumTRWSFactor::TernaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math3D::Tensor<float>& cost) :
  TernaryCumTRWSFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ double TernaryCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) { 

  return TernaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

TernaryCumTRWSRefFactor::TernaryCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math3D::Tensor<float>& cost) :
  TernaryCumTRWSFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ double TernaryCumTRWSRefFactor::compute_reparameterization(CumTRWSVar* var) { 

  return TernaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

SecondDiffCumTRWSFactor::SecondDiffCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, float lambda) :
  CumTRWSFactor(involved_vars), lambda_(lambda) {}

/*virtual*/ double SecondDiffCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  //this routine also updates reparameterization_
  uint idx = 0;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->cost();
  
  if (var == involved_var_[0]) {

    idx = 0;

    Math1D::Vector<double>& message = reparameterization_[idx];

    double base_cost = -param[1].max() - param[2].max() + 3*lambda_;

    for (int l1=0; l1 < int(nLabels1); l1++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {

        const double w2 = param[1][l2];

        const int part = - 2*l2 + l1;

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 + part; //- 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -w2 - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -w2 - param[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l1] = best;
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    Math1D::Vector<double>& message = reparameterization_[idx];

    double base_cost = -param[0].max() - param[2].max() + 3*lambda_;

    for (int l2=0; l2 < int(nLabels2); l2++) {
      
      double best = base_cost;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

        const double w1 = param[0][l1];

        const int part = - 2*l2 + l1;

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 + part; //- 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -w1 - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -w1 - param[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l2] = best;
    }    

  }
  else {
    assert(var == involved_var_[2]);

    idx = 2;

    Math1D::Vector<double>& message = reparameterization_[idx];

    double base_cost = -param[0].max() - param[1].max() + 3*lambda_;

    for (int l3=0; l3 < int(nLabels3); l3++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {

        const double w2 = param[1][l2];

        const int part = l3 - 2*l2;

        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = part + l1; //l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[0][l1] - w2;
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[0][l1] - w2 + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l3] = best;
    }
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;

  return offs;  
}


/********************/

FourthOrderCumTRWSFactor::FourthOrderCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                                   const Storage1D<Math3D::Tensor<float> >& cost) :
  CumTRWSFactor(involved_vars), cost_(cost) {
  assert(involved_vars.size() == 4);
}

/*virtual*/ double FourthOrderCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  double offs = 0.0;
  
  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();
  
  //this routine also updates reparameterization_

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->cost();

  uint idx = 0;

  if (var == involved_var_[0]) {

    idx = 0;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      const Math3D::Tensor<float>& cur_cost = cost_[l1];
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double sum3 = w2 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum3 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l1] = best;
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    Math1D::Vector<double>& message = reparameterization_[idx];
    
    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];
	
        const double w1 = param[0][l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double sum3 = w1 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels3; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum3 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l2] = best;
    }    
  }
  else if (var == involved_var_[2]) {
    
    idx = 2;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];
	
        const double w1 = param[0][l1];

        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double sum2 = w1 + param[1][l2];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum2 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l3] = best;
    }
  }
  else {
    
    assert(var == involved_var_[3]);

    idx = 3;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        const double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double sum2 = w1 + param[1][l2];

          for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum2 - param[2][l3];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l4] = best;
    }    
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;
  
  return offs;
}

/********************/

OneOfNCumTRWSFactor::OneOfNCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars) :
  CumTRWSFactor(involved_vars) {

  for (uint v=0; v < involved_vars.size(); v++)
    assert(involved_vars[v]->nLabels() == 2);
}

/*virtual*/ double OneOfNCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
      param[k].set_constant(0.0);
    }
  }

  Math1D::Vector<double>& message = reparameterization_[idx];

  assert(idx < nVars);

  double best_gain = 1e300;
  uint argmin = MAX_UINT;

  double sum = 0.0;

  for (uint i=0; i < nVars; i++) {

    if (involved_var_[i] == var) {
      message[0] = -param[i][0];
      message[1] = -param[i][1]; 
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

  message[0] += sum + param[argmin][0] - param[argmin][1];
  message[1] += sum;

  double offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  return offs;
}

/********************/

CardinalityCumTRWSFactor::CardinalityCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math1D::Vector<float>& cost) :
  CumTRWSFactor(involved_vars), cost_(cost) {

  for (uint v=0; v < involved_vars.size(); v++)
    assert(involved_vars[v]->nLabels() == 2);
}
  
/*virtual*/ double CardinalityCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  assert(var->nLabels() == 2);

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
      param[k].set_constant(0.0);
    }
  }

  assert(idx < nVars);

  Math1D::Vector<double>& message = reparameterization_[idx];
  message.set_constant(1e300);

  Math1D::NamedVector<double> rel_param(nVars-1,MAKENAME(rel_param));

  double offs = 0.0;

  uint next = 0;
  for (uint k=0; k < nVars; k++) {

    if (k != idx) {
      rel_param[next] = (-param[k][1]) - (-param[k][0]);
      offs += -param[k][0];

      next++;
    }
  }

  std::sort(rel_param.direct_access(), rel_param.direct_access() + nVars-1);
  
  double cum_sum = 0.0;

  for (uint c=0; c < nVars; c++) {

    double hyp0 = cum_sum + cost_[c];
    if (hyp0 < message[0])
      message[0] = hyp0;

    double hyp1 = cum_sum + cost_[c+1];
    if (hyp1 < message[1])
      message[1] = hyp1;

    if (c+1 < nVars) 
      cum_sum += rel_param[c];
  }

  for (uint l=0; l < 2; l++)
    message[l] += offs - param[idx][l];

  double total_offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  return total_offs;
}

/********************/

BILPCumTRWSFactor::BILPCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                                     int rhs_lower, int rhs_upper) :
  CumTRWSFactor(involved_vars), positive_(positive), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

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

/*virtual*/ double BILPCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
      param[k].set_constant(0.0);
    }
  }

  Math1D::Vector<double>& message = reparameterization_[idx];

  assert(idx < nVars);

  /**** forward ****/

  Math3D::NamedTensor<double> forward(range_,2,idx+1,MAKENAME(forward));
  Math2D::Matrix<double> forward_light(range_,idx+1);

  //init
  for (int sum=0; sum < range_; sum++) {

    forward_light(sum,0) = 1e100;
    for (int l=0; l < 2; l++) {
      forward(sum,l,0) = 1e100;
    }
  }

  forward(zero_offset_,0,0) = -param[0][0];
  forward_light(zero_offset_,0) = -param[0][0];
  const int init_mul = (positive_[0]) ? 1 : -1;
  if (int(zero_offset_)+init_mul >= 0
      && int(zero_offset_)+init_mul < range_) {
    forward(zero_offset_+init_mul,1,0) = -param[0][1];
    forward_light(zero_offset_+init_mul,0) = -param[0][1];
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
      forward_light(sum,v) = std::min(forward(sum,0,v), forward(sum,1,v));
    }
  }

  Math2D::Matrix<double> backward_light(range_,nVars);

  /**** backward ****/

  const uint last_var = nVars-1;

  //init
  for (int sum=0; sum < range_; sum++) 
    backward_light(sum,last_var) = 1e100;

  backward_light(zero_offset_,last_var) = -param[last_var][0];
  const int end_mul = (positive_[last_var]) ? 1 : -1;
  if (int(zero_offset_) + end_mul >= 0
      && int(zero_offset_) + end_mul < range_)
    backward_light(zero_offset_ + end_mul,last_var) = -param[last_var][1];

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
      	}

        if (hyp < best_prev)
          best_prev = hyp;
      }

      backward_light(sum,v) = best_prev;
    }
  }

  for (uint l=0; l < 2; l++) {

    double min_msg = 1e300;
    
    for (int s=0; s < (int) range_; s++) {
      
      double hyp = forward(s,l,idx);

      if (idx+1 < nVars) {

        double best_bwd = 1e300;

        const int diff = (s - zero_offset_);

        for (int r=rhs_lower_; r <= rhs_upper_; r++) {
          const int other = r + zero_offset_ - diff; 
	  
          if (other >= 0 && other < (int) range_) {

            best_bwd = std::min(best_bwd,backward_light(other,idx+1));
          }
        }

        hyp += best_bwd;
      }
      else {
        if (s < int(rhs_lower_ + zero_offset_) || s > int(rhs_upper_ + zero_offset_)) 
          hyp = 1e300;
      }

      if (hyp < min_msg)
        min_msg = hyp;
    }

    message[l] = min_msg;
  }

  double offs = message.min();

  if (fabs(offs) > 1e10) {

    std::cerr << "idx: " << idx << "/" << nVars << std::endl;
    std::cerr << "message: " << message << std::endl;
    for (uint v=0; v < nVars; v++) {
      std::cerr << "combined parameters[" << v <<  "]: " << param[v] << std::endl;
      std::cerr << "reparameterization[" << v << "]: " << reparameterization_[v] << std::endl;
    }

    std::cerr << "factor: ";
    for (uint v=0; v < involved_var_.size(); v++)
      std::cerr << involved_var_[v]->rank() << " ";
    std::cerr << std::endl;
    std::cerr << "var params:" << std::endl;
    for (uint v=0; v < involved_var_.size(); v++)
      std::cerr << involved_var_[v]->cost() << std::endl;
  }

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  return offs;
}

/********************/

CumFactorTRWS::CumFactorTRWS(uint nVars, uint nFactors) : var_(nVars), factor_(nFactors), labeling_(nVars,0), rank2var_(nVars,MAX_UINT),
                                                          nUsedVars_(0), nUsedFactors_(0), optimize_called_(false) {}

CumFactorTRWS::~CumFactorTRWS() {

  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];

  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

void CumFactorTRWS::add_var(const Math1D::Vector<float>& cost) {

  assert(!optimize_called_);

  if (nUsedVars_ == var_.size())
    var_.resize(uint(nUsedVars_*1.2)+4);

  assert(nUsedVars_ < var_.size());
  var_[nUsedVars_] = new CumTRWSVar(cost,nUsedVars_);
  rank2var_[nUsedVars_] = nUsedVars_;

  nUsedVars_++;
}

void CumFactorTRWS::add_factor(CumTRWSFactor* fac) {

  assert(!optimize_called_);

  if (factor_.size() == nUsedFactors_)
    factor_.resize(uint(nUsedFactors_*1.2)+4);

  factor_[nUsedFactors_] = fac;
  nUsedFactors_++;
}
  
void CumFactorTRWS::add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost) {

  assert(var1 < nUsedVars_);
  assert(var2 < nUsedVars_);

  Storage1D<CumTRWSVar*> vars(2);
  vars[0] = var_[var1];
  vars[1] = var_[var2];

  add_factor(new BinaryCumTRWSFactor(vars,cost));
}

void CumFactorTRWS::add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {

  Storage1D<CumTRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  CumTRWSFactor* new_fac = 0;

  if (!ref)
    new_fac = new TernaryCumTRWSFactor(vars,cost);
  else
    new_fac = new TernaryCumTRWSRefFactor(vars,cost);

  add_factor(new_fac);
}

void CumFactorTRWS::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  Storage1D<CumTRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  add_factor(new SecondDiffCumTRWSFactor(vars,lambda));
}

void CumFactorTRWS::add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                            const Storage1D<Math3D::Tensor<float> >& cost) {

  Storage1D<CumTRWSVar*> vars(4);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];
  vars[3] = var_[var4];

  add_factor(new FourthOrderCumTRWSFactor(vars,cost));
}

void CumFactorTRWS::add_one_of_n_factor(Math1D::Vector<uint>& var) {

  Storage1D<CumTRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  add_factor(new OneOfNCumTRWSFactor(vars));
}

void CumFactorTRWS::add_cardinality_factor(Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost) {

  if (var.size() == 1) {
    var_[var[0]]->add_cost(cost);
  }
  else {

    Storage1D<CumTRWSVar*> vars(var.size());
  
    for (uint k=0; k < var.size(); k++)
      vars[k] = var_[var[k]];

    add_factor(new CardinalityCumTRWSFactor(vars,cost));
  }
}

void CumFactorTRWS::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                                          int rhs_lower, int rhs_upper) {

  uint nUseful = 0;
  for (uint k=0; k < var.size(); k++) {
    
    const Math1D::Vector<double>& cur_cost = var_[var[k]]->cost();

    if (fabs(cur_cost[0] - cur_cost[1]) < 1e10)
      nUseful++;
    else {
      assert(cur_cost[0] < cur_cost[1]); 
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
    
    Storage1D<CumTRWSVar*> vars(nUseful);
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

    add_factor(new BILPCumTRWSFactor(vars,reduced_positive,rhs_lower,rhs_upper));
  }
  else {
    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;
  }
}

double CumFactorTRWS::optimize(uint nIter) {

  std::cerr << "efficient Cum-scheme" << std::endl;

  if (!optimize_called_) {
    size_t max_order = 0;

    for (uint f=0; f < nUsedFactors_; f++) {
      factor_[f]->compute_rank_range();
      max_order = std::max(max_order,factor_[f]->involved_vars().size());
    }
    
    std::cerr << "maximal order of all factors: " << max_order << std::endl;
    
    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->set_up_chains();
    }
  }
  optimize_called_ = true;

  uint arg_min;

  double bound = 1e300;

  {
    size_t message_effort = 0;
    for (uint f=0; f < nUsedFactors_; f++) {
      message_effort += 2* (factor_[f]->involved_vars().size()-1) * (factor_[f]->involved_vars().size());
    } 
    message_effort *= nIter;
    std::cerr << "predicted message effort: " << message_effort << std::endl;
  }

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "******* iteration " << iter << " *********" << std::endl;

    //forward
    double forward_lower = 0.0;

    for (uint i=0; i < nUsedVars_; i++) {

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors
      
      CumTRWSVar* cur_var = var_[rank2var_[i]];
      
      Storage1D<CumTRWSFactor*> adj_factors = cur_var->adjacent_factor();
      
      for (uint k=0; k < adj_factors.size(); k++) {

        if (i != adj_factors[k]->min_rank()) {

          double cur_offs = adj_factors[k]->compute_reparameterization(cur_var);
	
	  if (i == adj_factors[k]->max_rank())
	    forward_lower += cur_offs;
        }
      }

      forward_lower += cur_var->average(arg_min);
      labeling_[i] = arg_min;
    }

    //std::cerr << "fwd bound: " << forward_lower << std::endl;

    //backward
    double backward_lower = 0.0;
    for (int i = nUsedVars_-1; i >= 0; i--) {

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors
      
      CumTRWSVar* cur_var = var_[rank2var_[i]];

      Storage1D<CumTRWSFactor*> adj_factors = cur_var->adjacent_factor();
      
      for (uint k=0; k < adj_factors.size(); k++) {

        if (i != int(adj_factors[k]->max_rank())) {	
          double cur_offs = adj_factors[k]->compute_reparameterization(cur_var);
	
          if (i == int(adj_factors[k]->min_rank()))
            backward_lower += cur_offs;
        }
      }

      backward_lower += cur_var->average(arg_min);
      labeling_[i] = arg_min;
    }

    std::cerr << "bwd bound: " << backward_lower << std::endl;

    bound = backward_lower;
  }

  size_t message_effort = 0;
  for (uint f=0; f < nUsedFactors_; f++) {
    message_effort += 2* (factor_[f]->involved_vars().size()-1) * (factor_[f]->involved_vars().size());
  } 
  message_effort *= nIter;
  std::cerr << "message effort: " << message_effort << std::endl;

  return bound;
}

const Math1D::Vector<uint>& CumFactorTRWS::labeling() {
  return labeling_;
}



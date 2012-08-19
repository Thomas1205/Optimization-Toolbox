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

  //std::cerr << "avg" << std::endl;

  const uint nLabels = cost_.size();
  const uint nFactors = adjacent_factor_.size();

  for (uint l=0; l < nLabels; l++)
    cum_cost_[l] = cost_[l]; //convert float -> double

  //std::cerr << "cost initialized" << std::endl;

  for (uint k=0; k < nFactors; k++)
    cum_cost_ += adjacent_factor_[k]->reparameterization(this);

  //std::cerr << "cost added" << std::endl;

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


BinaryCumTRWSFactorBase::BinaryCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars)
  : CumTRWSFactor(involved_vars) {
  assert(involved_vars.size() == 2);
}

double BinaryCumTRWSFactorBase::compute_reparameterization(CumTRWSVar* var, const Math2D::Matrix<float>& cost) {

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

        double hyp = cost(l1,l2) - param[1][l2];

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

        double hyp = cost(l1,l2) - param[0][l1];

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

BinaryCumTRWSFactor::BinaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                         const Math2D::Matrix<float>& cost) : 
  BinaryCumTRWSFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ double BinaryCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {
  
  return BinaryCumTRWSFactorBase::compute_reparameterization(var,cost_);

  // double offs = 0.0;

  // const uint nLabels1 = involved_var_[0]->nLabels();
  // const uint nLabels2 = involved_var_[1]->nLabels();

  // //this routine also updates reparameterization_
  // Storage1D<Math1D::Vector<double> > param = reparameterization_;
  // for (uint k=0; k < involved_var_.size(); k++) {
  //   param[k] -= involved_var_[k]->cost();
  // }

  // uint idx = 0;

  // if (var == involved_var_[0]) {

  //   idx = 0;

  //   for (uint l1 = 0; l1 < nLabels1; l1++) {

  //     double best = 1e300;

  //     for (uint l2 = 0; l2 < nLabels2; l2++) {

  //       double hyp = cost_(l1,l2) - param[1][l2];

  //       if (hyp < best)
  //         best = hyp;
  //     }

  //     reparameterization_[idx][l1] = best;
  //   }
  // }
  // else {
  //   assert(var == involved_var_[1]);

  //   idx = 1;

  //   for (uint l2 = 0; l2 < nLabels2; l2++) {

  //     double best = 1e300;

  //     for (uint l1 = 0; l1 < nLabels1; l1++) {

  //       double hyp = cost_(l1,l2) - param[0][l1];

  //       if (hyp < best)
  //         best = hyp;
  //     }

  //     reparameterization_[idx][l2] = best;
  //   }
  // }

  // double msg_offs = reparameterization_[idx].min();
  // offs += msg_offs;

  // for (uint l=0; l < reparameterization_[idx].size(); l++)
  //   reparameterization_[idx][l] -= msg_offs;
  
  // return offs;
}

/********************/

BinaryCumTRWSRefFactor::BinaryCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                               const Math2D::Matrix<float>& cost) : 
  BinaryCumTRWSFactorBase(involved_vars), cost_(cost) {
  assert(involved_vars.size() == 2);
}

/*virtual*/ double BinaryCumTRWSRefFactor::compute_reparameterization(CumTRWSVar* var) {
  
  return BinaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
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
	  
          double hyp = cost(l1,l2,l3) - w2 /*param[1][l2]*/ - param[2][l3];
	  
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
	  
          double hyp = cost(l1,l2,l3) - w1 /*param[0][l1]*/ - param[2][l3];
	  
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
	  
          double hyp = cost(l1,l2,l3) - w1 /*param[0][l1]*/ - param[1][l2];
	  
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
            //hyp = -param[1][l2] - param[2][l3];
            hyp = -w2 - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            //hyp = -param[1][l2] - param[2][l3] + lambda_;
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
            //hyp = -param[0][l1] - param[2][l3];
            hyp = -w1 - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            //hyp = -param[0][l1] - param[2][l3] + lambda_;
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
            //hyp = -param[0][l1] - param[1][l2];
            hyp = -param[0][l1] - w2;
          }
          else if (abs(so_diff) <= 1) {
            //hyp = -param[0][l1] - param[1][l2] + lambda_;
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

FourthOrderCumTRWSFactorBase::FourthOrderCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars)
  : CumTRWSFactor(involved_vars) {
  assert(involved_vars.size() == 4);
}

double FourthOrderCumTRWSFactorBase::compute_reparameterization(CumTRWSVar* var, const Storage1D<Math3D::Tensor<float> >& cost) {

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

      const Math3D::Tensor<float>& cur_cost = cost[l1];
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double sum3 = w2 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              //- param[1][l2] - param[2][l3] - param[3][l4];
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

        const Math3D::Tensor<float>& cur_cost = cost[l1];
	
        const double w1 = param[0][l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double sum3 = w1 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels3; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              //- param[0][l1] - param[2][l3] - param[3][l4];
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

        const Math3D::Tensor<float>& cur_cost = cost[l1];
	
        const double w1 = param[0][l1];

        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double sum2 = w1 + param[1][l2];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              //- param[0][l1] - param[1][l2] - param[3][l4];
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

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double sum2 = w1 + param[1][l2];

          for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              //- param[0][l1] - param[1][l2] - param[2][l3];
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


FourthOrderCumTRWSFactor::FourthOrderCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                                   const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderCumTRWSFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ double FourthOrderCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  return  FourthOrderCumTRWSFactorBase::compute_reparameterization(var,cost_);

  // double offs = 0.0;
  
  // const uint nLabels1 = involved_var_[0]->nLabels();
  // const uint nLabels2 = involved_var_[1]->nLabels();
  // const uint nLabels3 = involved_var_[2]->nLabels();
  // const uint nLabels4 = involved_var_[3]->nLabels();
  
  // //this routine also updates reparameterization_

  // Storage1D<Math1D::Vector<double> > param = reparameterization_;
  // for (uint k=0; k < involved_var_.size(); k++)
  //   param[k] -= involved_var_[k]->cost();

  // uint idx = 0;

  // if (var == involved_var_[0]) {

  //   idx = 0;

  //   Math1D::Vector<double>& message = reparameterization_[idx];

  //   for (uint l1 = 0; l1 < nLabels1; l1++) {
      
  //     double best = 1e300;

  //     const Math3D::Tensor<float>& cur_cost = cost_[l1];
      
  //     for (uint l2 = 0; l2 < nLabels2; l2++) {

  //       const double w2 = param[1][l2];
	
  //       for (uint l3 = 0; l3 < nLabels3; l3++) {

  //         const double sum3 = w2 + param[2][l3];

  //         for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
  //           double hyp = cur_cost(l2,l3,l4) 
  //             //- param[1][l2] - param[2][l3] - param[3][l4];
  //             - sum3 - param[3][l4];
	  
  //           if (hyp < best)
  //             best = hyp;
  //         }
  //       }
  //     }
      
  //     message[l1] = best;
  //   }
  // }
  // else if (var == involved_var_[1]) {

  //   idx = 1;

  //   Math1D::Vector<double>& message = reparameterization_[idx];
    
  //   for (uint l2 = 0; l2 < nLabels1; l2++) {
      
  //     double best = 1e300;
      
  //     for (uint l1 = 0; l1 < nLabels1; l1++) {

  //       const Math3D::Tensor<float>& cur_cost = cost_[l1];
	
  //       const double w1 = param[0][l1];

  //       for (uint l3 = 0; l3 < nLabels3; l3++) {

  //         const double sum3 = w1 + param[2][l3];

  //         for (uint l4 = 0; l4 < nLabels3; l4++) {
	  
  //           double hyp = cur_cost(l2,l3,l4) 
  //             //- param[0][l1] - param[2][l3] - param[3][l4];
  //             - sum3 - param[3][l4];
	  
  //           if (hyp < best)
  //             best = hyp;
  //         }
  //       }
  //     }
      
  //     message[l2] = best;
  //   }    
  // }
  // else if (var == involved_var_[2]) {
    
  //   idx = 2;

  //   Math1D::Vector<double>& message = reparameterization_[idx];

  //   for (uint l3 = 0; l3 < nLabels3; l3++) {
      
  //     double best = 1e300;
      
  //     for (uint l1 = 0; l1 < nLabels1; l1++) {

  //       const Math3D::Tensor<float>& cur_cost = cost_[l1];
	
  //       const double w1 = param[0][l1];

  //       for (uint l2 = 0; l2 < nLabels2; l2++) {

  //         const double sum2 = w1 + param[1][l2];

  //         for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
  //           double hyp = cur_cost(l2,l3,l4) 
  //             //- param[0][l1] - param[1][l2] - param[3][l4];
  //             - sum2 - param[3][l4];
	  
  //           if (hyp < best)
  //             best = hyp;
  //         }
  //       }
  //     }
      
  //     message[l3] = best;
  //   }
  // }
  // else {
    
  //   assert(var == involved_var_[3]);

  //   idx = 3;

  //   Math1D::Vector<double>& message = reparameterization_[idx];

  //   for (uint l4 = 0; l4 < nLabels4; l4++) {
      
  //     double best = 1e300;
      
  //     for (uint l1 = 0; l1 < nLabels1; l1++) {

  //       const Math3D::Tensor<float>& cur_cost = cost_[l1];

  //       const double w1 = param[0][l1];
	
  //       for (uint l2 = 0; l2 < nLabels2; l2++) {

  //         const double sum2 = w1 + param[1][l2];

  //         for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
  //           double hyp = cur_cost(l2,l3,l4) 
  //             //- param[0][l1] - param[1][l2] - param[2][l3];
  //             - sum2 - param[2][l3];
	  
  //           if (hyp < best)
  //             best = hyp;
  //         }
  //       }
  //     }
      
  //     message[l4] = best;
  //   }    
  // }

  // double msg_offs = reparameterization_[idx].min();
  // offs += msg_offs;

  // for (uint l=0; l < reparameterization_[idx].size(); l++)
  //   reparameterization_[idx][l] -= msg_offs;
  
  // return offs;
}

/********************/

FourthOrderCumTRWSRefFactor::FourthOrderCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                                         const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderCumTRWSFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ double FourthOrderCumTRWSRefFactor::compute_reparameterization(CumTRWSVar* var) {

  return  FourthOrderCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

TwoLevelCumTRWSFactor::TwoLevelCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                             const Math2D::Matrix<uint>& exception, double exception_cost, double std_cost) :
  CumTRWSFactor(involved_vars), exception_(exception), exception_cost_(exception_cost), std_cost_(std_cost) {

  const uint nVars = involved_vars.size();
  assert(nVars == exception.xDim());

  if (exception_cost > std_cost) {

    for (uint k1=0; k1 < exception_.yDim(); k1++) {
      for (uint k2=0; k2 < exception_.yDim(); k2++) {

        uint nDiffs = 0;
        for (uint l=0; l < nVars; l++)
          if (exception_(l,k1) != exception_(l,k2))
            nDiffs++;

        if (nDiffs == 1) {
	  
          INTERNAL_ERROR << "when exceptions are costlier, all exceptions need to have a hamming distance of 2 to each other. Exiting.." 
                         << std::endl;
          exit(1);
        }
      }
    }
  }
}

/*virtual*/ double TwoLevelCumTRWSFactor::compute_reparameterization(CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  const uint yDim = exception_.yDim();

  uint idx = MAX_UINT;

  Math1D::Vector<double> mincost(nVars,1e300);
  Math1D::Vector<uint> argmin(nVars,MAX_UINT);

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
      param[k].set_constant(0.0);
    }

    for (uint l=0; l < involved_var_[k]->nLabels(); l++) {

      if (-param[k][l] < mincost[k]) {
        mincost[k] = -param[k][l];
        argmin[k] = l;
      }
    }
  }

  double cost_base = mincost.sum() - mincost[idx];

  Math1D::Vector<double>& message = reparameterization_[idx];

  std::vector<std::vector<uint> > label_exceptions(involved_var_[idx]->nLabels());

  for (uint exc = 0; exc < yDim; exc++) {
    label_exceptions[exception_(idx,exc)].push_back(exc);
  }

  for (uint l=0; l < involved_var_[idx]->nLabels(); l++) {
    
    if (exception_cost_ <= std_cost_) {

      double best_entry = cost_base - param[idx][l] + std_cost_;

      for (uint exc = 0; exc < yDim; exc++) {
        if (exception_(idx,exc) == l) {
          double hyp_entry = exception_cost_;
          for (uint k=0; k < nVars; k++)
            hyp_entry -= param[k][exception_(k,idx)];
				  
          if (hyp_entry < best_entry)
            best_entry = hyp_entry;
        }
      }

      message[l] = best_entry;
    }
    else {

      double std_entry = cost_base - param[idx][l] + std_cost_;
      mincost[idx] = -param[idx][l];
      argmin[idx] = l;

      bool argmin_is_exc = false;

      double exc_entry = 1e300;

      // for (uint exc = 0; exc < yDim; exc++) {
      //   if (exception_(idx,exc) == l) {

      for (uint k=0; k < label_exceptions[l].size(); k++) {

        uint exc = label_exceptions[l][k];
        
        if (true) {

          //check if the standard argmin is an exception
          if (!argmin_is_exc) {
            argmin_is_exc = true;

            for (uint k=0; k < nVars; k++) {

              if (exception_(k,exc) != argmin[k]) {
                argmin_is_exc = false;
                break;
              }
            }
          }

          // double hyp_entry = exception_cost_;
          // for (uint k=0; k < nVars; k++)
          //   hyp_entry -= param[k][exception_(k,exc)];
				  
          // if (hyp_entry < exc_entry)
          //   exc_entry = hyp_entry;
        }
      }

      if (argmin_is_exc) {

        for (uint k=0; k < label_exceptions[l].size(); k++) {

          uint exc = label_exceptions[l][k];

          double hyp_entry = exception_cost_;
          for (uint k=0; k < nVars; k++)
            hyp_entry -= param[k][exception_(k,exc)];
				  
          if (hyp_entry < exc_entry)
            exc_entry = hyp_entry;
        }

        double best_add = 0.0;
        double best_change = 1e300;
        uint arg_b2b = MAX_UINT;

        for (uint k=0; k < nVars; k++) {

          if (k != idx) {

#if 1
            Math1D::Vector<double> temp = param[k];
            temp *= -1.0;
            std::sort(temp.direct_access(),temp.direct_access()+temp.size());

            if (temp[1] - temp[0] < best_change) {
              best_change = temp[1] - temp[0];
              best_add = temp[1];
              arg_b2b = k;
            }
#else
            double cur_best = 1e300;
            double cur_2nd_best = 1e300;

            const uint nLabels = param[k].size();
            for (uint l=0; l < nLabels; l++) {
              double cur = -param[k][l];

              if (cur < cur_best) {
                cur_2nd_best = cur_best;
                cur_best = cur;
              }
              else if (cur < cur_2nd_best)
                cur_2nd_best = cur;
            }
            const double  diff = cur_2nd_best - cur_best;
            if (diff < best_change) {
              best_change = diff;
              best_add = cur_2nd_best;
              arg_b2b = k;
            }
#endif       
          }
        }

        std_entry += best_add - mincost[arg_b2b];
      }

      message[l] = std::min(std_entry, exc_entry);
    }
  }


  const double offs = message.min();

  for (uint k=0; k < message.size(); k++)
    message[k] -= offs;

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

  //DEBUG
  // if (involved_var_.size() == 23 && min_rank_ == 3084 && max_rank_ == 13012) {
  //   std::cerr << "1ON, C, intermediate msg0: " << message[0] << std::endl;
  // }
  //END_DEBUG

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

  //std::cerr << "nVars: " << nVars << std::endl;

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

  //DEBUG
  // if (involved_var_.size() == 2 && involved_var_[0]->rank() == 48905 && involved_var_[1]->rank() == 48906) {

  //   std::cerr << "CARD, idx: " << idx << std::endl;
  //   std::cerr << "cost: " << cost_ << std::endl;

  //   for (uint i=0; i < 2; i++) {
  //     std::cerr << "reparameterization [" << i << "]: " << reparameterization_[i] << std::endl;
  //   }

  //   for (uint i=0; i < 2; i++) {
  //     std::cerr << "params[" << i << "]: " << param[i] << std::endl;
  //   }
    
  //   std::cerr << "computed message: " << message << std::endl;
  // }
  //END_DEBUG

  double total_offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  //std::cerr << "end card. message comp." << std::endl; 

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

  // std::cerr << "positive: " << positive_ << std::endl;
  // std::cerr << "rhs_lower: " << rhs_lower << ", rhs_upper: " << rhs_upper << std::endl;
  // std::cerr << "lower bound: " << lower_bound << std::endl;

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

  //std::cerr << "nVars: " << nVars << std::endl;

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

  //std::cerr << "zero offset: " << zero_offset_ << std::endl;
  //std::cerr << "pos. vector: " << positive_ << std::endl;

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

    //std::cerr << "v: " << v << std::endl;

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

    //std::cerr << "v: " << v << std::endl;

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

    //std::cerr << "l: " << l << std::endl;

    double min_msg = 1e300;
    
    for (int s=0; s < (int) range_; s++) {

      //std::cerr << "s: " << s << std::endl;
      
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

  //std::cerr << "comp offs." << std::endl;

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

    //DEBUG
    exit(1);
    //END_DEBUG
  }

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  //std::cerr << "leaving" << std::endl;

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
  
void CumFactorTRWS::add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref) {

  assert(var1 < nUsedVars_);
  assert(var2 < nUsedVars_);

  Storage1D<CumTRWSVar*> vars(2);
  vars[0] = var_[var1];
  vars[1] = var_[var2];

  CumTRWSFactor* newFac;
  if (ref)
    newFac = new BinaryCumTRWSRefFactor(vars,cost);
  else
    newFac = new BinaryCumTRWSFactor(vars,cost);
      
  add_factor(newFac);
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

  // assert(nUsedFactors_ < factor_.size());

  // factor_[nUsedFactors_] = new SecondDiffCumTRWSFactor(vars,lambda);

  // nUsedFactors_++;
}

void CumFactorTRWS::add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                            const Storage1D<Math3D::Tensor<float> >& cost, bool ref) {

  Storage1D<CumTRWSVar*> vars(4);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];
  vars[3] = var_[var4];

  CumTRWSFactor* newFac;

  if (ref)
    newFac = new FourthOrderCumTRWSRefFactor(vars,cost);
  else
    newFac = new FourthOrderCumTRWSFactor(vars,cost);

  add_factor(newFac);
}

void CumFactorTRWS::add_two_level_factor(Math1D::Vector<uint>& var, Math2D::Matrix<uint>& exceptions,
                                         double exception_cost, double standard_cost) {

  Storage1D<CumTRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  add_factor(new TwoLevelCumTRWSFactor(vars,exceptions,exception_cost,standard_cost));
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

  // if (nUseful != var.size())
  //   std::cerr << "removed " << (var.size() - nUseful) << " / " << var.size() << " vars for BILP node" << std::endl;

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

  //DEBUG (for ICML)
  //std::ofstream of("trws.dat");
  //END_DEBUG

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "******* iteration " << iter << " *********" << std::endl;

    //forward
    double forward_lower = 0.0;

    for (uint i=0; i < nUsedVars_; i++) {

      //std::cerr << "i: " << i << std::endl;

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors
      
      CumTRWSVar* cur_var = var_[rank2var_[i]];
      
      Storage1D<CumTRWSFactor*> adj_factors = cur_var->adjacent_factor();
      
      for (uint k=0; k < adj_factors.size(); k++) {
	
        //std::cerr << "k: " << k << std::endl;

        if (i != adj_factors[k]->min_rank()) {

          double cur_offs = adj_factors[k]->compute_reparameterization(cur_var);
	
	  if (i == adj_factors[k]->max_rank())
	    forward_lower += cur_offs;
        }
      }

      //std::cerr << "calling avrg" << std::endl;

      forward_lower += cur_var->average(arg_min);
      labeling_[i] = arg_min;

      //std::cerr << "done" << std::endl;
    }

    //DEBUG (for ICML)
    //of << ((iter-1)*effort_per_iter + effort_per_iter/2) << " " << forward_lower << std::endl;
    //END_DEBUG

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

    //DEBUG (for ICML)
    //of << (effort_per_iter * iter) << " " << backward_lower << std::endl;
    //END_DEBUG
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


/*******************************************************************************************************/


TRWSVar::TRWSVar(const Math1D::Vector<double>& cost, uint rank) : cost_(cost), rank_(rank) {}

void TRWSVar::add_cost(const Math1D::Vector<double>& add_cost) {
  cost_ += add_cost;
}

void TRWSVar::set_up_chains() {

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

    //DEBUG
    // if (nChains_ != nChainlessChains) {
    //   std::cerr << "ERROR: NOT equivalent: " << nChains_ << " vs " << nChainlessChains << std::endl;
    //   std::cerr << "nOutgoing: " << nOutgoing << ", nIncoming: " << nIncoming << ", nMiddle: " << nMiddle << std::endl;
    // }

    // if (nOutgoing < nIncoming) {

    //   for (uint k=0; k < adjacent_factor_.size(); k++) {

    // 	if (adjacent_factor_[k]->min_rank() == rank_ && adjacent_factor_[k]->prev_var() != this) {
    // 	  std::cerr << "has prev var: " << (adjacent_factor_[k]->prev_var() != 0) << std::endl;
    // 	}
    //   }
    // }
    //END_DEBUG

    nChains_ = nChainlessChains;
  }


  // //NOTE: this is valid only for strict monotonocity (no implementation of relaxed mono. in this class)
  // for (uint k=0; k < adjacent_factor_.size(); k++) {
  //   //count chains, encode roles
  //   if (adjacent_factor_[k]->next_var() == this) {

  //     nChains_++;
  //   }
  //   else if (adjacent_factor_[k]->prev_var() == this) {
  //     //NOTE: this chain is counted elsewhere
  //   }
  //   else {
  //     nChains_++;
  //   }
  // }

  //devide cost by number of chains
  for (uint l=0; l < cost_.size(); l++)
    cost_[l] /= nChains_;
}

void TRWSVar::add_factor(TRWSFactor* factor) {

  uint size = adjacent_factor_.size();
  
  adjacent_factor_.resize(size+1);
  adjacent_factor_[size] = factor;
}

double TRWSVar::average_forward(double& bound, uint& arg_min) {

  Math1D::Vector<double> sum(cost_.size());

  for (uint l=0; l < cost_.size(); l++)
    sum[l] = nChains_ * cost_[l];  

  Math1D::Vector<double> msg(cost_.size(),0.0);

  for (uint k=0; k < adjacent_factor_.size(); k++) {

    uint min_rank = adjacent_factor_[k]->min_rank();

    if (rank_ != min_rank) {

      double cur_offs = adjacent_factor_[k]->compute_message_and_reparameterize(this,msg);

      sum += msg;

      if (rank_ == adjacent_factor_[k]->max_rank())
        bound += cur_offs;
    }
  }

  double offs = 1e300;

  for (uint l=0; l < cost_.size(); l++) {

    if (sum[l] < offs) {
      offs = sum[l];
      arg_min = l;
    }
  }

  for (uint l=0; l < cost_.size(); l++)
    cost_[l] = (sum[l] - offs) / nChains_;
  
  return offs;
}

double TRWSVar::average_backward(double& bound, uint& arg_min) {

  Math1D::Vector<double> sum(cost_.size());

  for (uint l=0; l < cost_.size(); l++)
    sum[l] = nChains_ * cost_[l];  

  Math1D::Vector<double> msg(cost_.size(),0.0);

  for (uint k=0; k < adjacent_factor_.size(); k++) {

    uint max_rank = adjacent_factor_[k]->max_rank();

    if (rank_ != max_rank) {

      double cur_offs = adjacent_factor_[k]->compute_message_and_reparameterize(this,msg);

      sum += msg;

      if (rank_ == adjacent_factor_[k]->min_rank())
        bound += cur_offs;
    }
  }

  double offs = 1e300;

  for (uint l=0; l < cost_.size(); l++) {

    if (sum[l] < offs) {
      offs = sum[l];
      arg_min = l;
    }
  }

  for (uint l=0; l < cost_.size(); l++)
    cost_[l] = (sum[l] - offs) / nChains_;

  return offs;
}

uint TRWSVar::rank() const {
  return rank_;
}

void TRWSVar::set_rank(uint rank) {
  rank_ = rank;
}

size_t TRWSVar::nLabels() const {
  return cost_.size();
}

const Storage1D<TRWSFactor*>& TRWSVar::adjacent_factor() const {
  return adjacent_factor_;
}

const Math1D::Vector<double>& TRWSVar::cost() const {
  return cost_;
}

/*********************/

TRWSFactor::TRWSFactor(const Storage1D<TRWSVar*>& involved_vars) : 
  involved_var_(involved_vars) {

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
  
uint TRWSFactor::min_rank() const {
  return min_rank_;
}

uint TRWSFactor::max_rank() const {
  return max_rank_;
}

void TRWSFactor::compute_rank_range() {

  min_rank_ = MAX_UINT;
  max_rank_ = 0;

  uint nVars = involved_var_.size();

  for (uint v=0; v < nVars; v++) {

    min_rank_ = std::min(min_rank_,involved_var_[v]->rank());
    max_rank_ = std::max(max_rank_,involved_var_[v]->rank());
  }
}

const Storage1D<TRWSVar*>& TRWSFactor::involved_vars() const {
  return involved_var_;
}

const Math1D::Vector<double>& TRWSFactor::reparameterization(const TRWSVar* var) const {

  for (uint k=0; k < involved_var_.size(); k++) {

    if (involved_var_[k] == var)
      return reparameterization_[k];
  }

  assert(false);
  return reparameterization_[0];
}

/***********************************/

BinaryTRWSFactor::BinaryTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                                   const Math2D::Matrix<float>& cost) : TRWSFactor(involved_vars), cost_(cost) {}
  
/*virtual */
double BinaryTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

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

    message.resize(nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double hyp = cost_(l1,l2) - param[1][l2];

        if (hyp < best)
          best = hyp;
      }

      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else {
    assert(var == involved_var_[1]);

    idx = 1;

    message.resize(nLabels2);

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      double best = 1e300;

      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double hyp = cost_(l1,l2) - param[0][l1];

        if (hyp < best)
          best = hyp;
      }

      message[l2] = best - reparameterization_[1][l2];
    }
  }

  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  reparameterization_[idx] += message;
  
  return offs;
}

/**************************************/

TernaryTRWSFactorBase::TernaryTRWSFactorBase(const Storage1D<TRWSVar*>& involved_vars) 
  : TRWSFactor(involved_vars) {}
  
double TernaryTRWSFactorBase::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message,
                                                                 const Math3D::Tensor<float>& cost) {

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

    message.resize(nLabels1);
    
    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - w2 /*param[1][l2]*/ - param[2][l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else if (var == involved_var_[1]) {

    message.resize(nLabels2);
    
    idx = 1;

    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double w1 = param[0][l1];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - w1 - param[2][l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l2] = best - reparameterization_[1][l2];
    }    
  }
  else {
    assert(var == involved_var_[2]);
    
    message.resize(nLabels3);

    idx = 2;

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        double w1 = param[0][l1];	

        for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
          double hyp = cost(l1,l2,l3) - w1 - param[1][l2];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l3] = best - reparameterization_[2][l3];
    }
  }

  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  reparameterization_[idx] += message;
  
  return offs;  
}

/******************/

TernaryTRWSFactor::TernaryTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                                     const Math3D::Tensor<float>& cost) :
  TernaryTRWSFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ double TernaryTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

  return TernaryTRWSFactorBase::compute_message_and_reparameterize(var,message,cost_);
}

/******************/

TernaryTRWSRefFactor::TernaryTRWSRefFactor(const Storage1D<TRWSVar*>& involved_vars,
                                           const Math3D::Tensor<float>& cost) :
  TernaryTRWSFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ double TernaryTRWSRefFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

  return TernaryTRWSFactorBase::compute_message_and_reparameterize(var,message,cost_);
}

/******************/

SecondDiffTRWSFactor::SecondDiffTRWSFactor(const Storage1D<TRWSVar*>& involved_vars, float lambda) :
  TRWSFactor(involved_vars), lambda_(lambda) {}

/*virtual */
double SecondDiffTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

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

    message.resize(nLabels1);

    double base_cost = -param[1].max() - param[2].max() + 3*lambda_;

    for (int l1=0; l1 < int(nLabels1); l1++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {

        double w2 = param[1][l2];

        int part = - 2*l2 + l1;

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 + part;
	  
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

      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    message.resize(nLabels2);

    double base_cost = -param[0].max() - param[2].max() + 3*lambda_;

    for (int l2=0; l2 < int(nLabels2); l2++) {
      
      double best = base_cost;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

        double w1 = param[0][l1];

        int part = - 2*l2 + l1;

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 + part;
	  
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

      message[l2] = best - reparameterization_[1][l2];
    }    

  }
  else {
    assert(var == involved_var_[2]);

    idx = 2;

    message.resize(nLabels3);

    double base_cost = -param[0].max() - param[1].max() + 3*lambda_;

    for (int l3=0; l3 < int(nLabels3); l3++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {

        double w2 = param[1][l2];

        int part = l3 - 2*l2;

        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = part + l1; 
	  
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

      message[l3] = best - reparameterization_[2][l3];
    }
  }


  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  reparameterization_[idx] += message;
  
  return offs;  
}


/**************************************/

FourthOrderTRWSFactor::FourthOrderTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                                             const Storage1D<Math3D::Tensor<float> >& cost) :
  TRWSFactor(involved_vars), cost_(cost) {}

/*virtual */
double FourthOrderTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

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

    message.resize(nLabels1);
    
    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          double w23  = w2 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - w23 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    message.resize(nLabels2);
    
    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        double w1 = param[0][l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          double w13 = w1 + param[2][l3];

          for (uint l4 = 0; l4 < nLabels3; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - w13 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l2] = best - reparameterization_[1][l2];
    }    
  }
  else if (var == involved_var_[2]) {
    
    message.resize(nLabels3);

    idx = 2;

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          double w12 = w1 + param[1][l2];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - w12 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l3] = best - reparameterization_[2][l3];
    }
  }
  else {
    
    assert(var == involved_var_[3]);

    message.resize(nLabels4);

    idx = 3;

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          double w12 = w1 + param[1][l2];

          for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - w12 - param[2][l3];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l4] = best - reparameterization_[3][l4];
    }    
  }

  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  reparameterization_[idx] += message;
  
  return offs;
}

/**************************************/

TwoLevelTRWSFactor::TwoLevelTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                                       const Math2D::Matrix<uint>& exception, double exception_cost, double std_cost) :
  TRWSFactor(involved_vars), exception_(exception), exception_cost_(exception_cost), std_cost_(std_cost) {

  const uint nVars = involved_vars.size();
  assert(nVars == exception.xDim());

  if (exception_cost > std_cost) {

    for (uint k1=0; k1 < exception_.yDim(); k1++) {
      for (uint k2=0; k2 < exception_.yDim(); k2++) {

        uint nDiffs = 0;
        for (uint l=0; l < nVars; l++)
          if (exception_(l,k1) != exception_(l,k2))
            nDiffs++;

        if (nDiffs == 1) {
	  
          INTERNAL_ERROR << "when exceptions are costlier, all exceptions need to have a hamming distance of 2 to each other. Exiting.." 
                         << std::endl;
          exit(1);
        }
      }
    }
  }
}

/*virtual */
double TwoLevelTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Math1D::Vector<double> mincost(nVars,1e300);
  Math1D::Vector<uint> argmin(nVars,MAX_UINT);

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else
      idx = k;

    for (uint l=0; l < involved_var_[k]->nLabels(); l++) {

      if (-param[k][l] < mincost[k]) {
        mincost[k] = -param[k][l];
        argmin[k] = l;
      }
    }
  }

  double cost_base = mincost.sum() - mincost[idx];

  message.resize(involved_var_[idx]->nLabels());

  for (uint l=0; l < involved_var_[idx]->nLabels(); l++) {

    
    if (exception_cost_ <= std_cost_) {

      double best_entry = cost_base - param[idx][l] + std_cost_;

      for (uint exc = 0; exc < exception_.yDim(); exc++) {
        if (exception_(idx,exc) == l) {
          double hyp_entry = exception_cost_;
          for (uint k=0; k < nVars; k++)
            hyp_entry -= param[k][exception_(k,idx)];
				  
          if (hyp_entry < best_entry)
            best_entry = hyp_entry;
        }
      }

      message[l] = best_entry;
    }
    else {

      double std_entry = cost_base - param[idx][l] + std_cost_;
      mincost[idx] = -param[idx][l];
      argmin[idx] = l;

      bool argmin_is_exc = false;

      double exc_entry = 1e300;

      for (uint exc = 0; exc < exception_.yDim(); exc++) {
        if (exception_(idx,exc) == l) {

          //check if the standard argmin is an exception
          if (!argmin_is_exc) {
            argmin_is_exc = true;

            for (uint k=0; k < nVars; k++) {

              if (exception_(k,exc) != argmin[k]) {
                argmin_is_exc = false;
                break;
              }
            }
          }

          double hyp_entry = exception_cost_;
          for (uint k=0; k < nVars; k++)
            hyp_entry -= param[k][exception_(k,exc)];
				  
          if (hyp_entry < exc_entry)
            exc_entry = hyp_entry;
        }
      }      

      if (argmin_is_exc) {

        double best_add = 0.0;
        double best_change = 1e300;
        uint arg_b2b = MAX_UINT;

        for (uint k=0; k < nVars; k++) {

          if (k != idx) {
            Math1D::Vector<double> temp = param[k];
            temp *= -1.0;
            std::sort(temp.direct_access(),temp.direct_access()+temp.size());

            if (temp[1] - temp[0] < best_change) {
              best_change = temp[1] - temp[0];
              best_add = temp[1];
              arg_b2b = k;
            }
          }
        }

        std_entry += best_add - mincost[arg_b2b];
      }

      message[l] = std::min(std_entry, exc_entry);
    }
  }


  double offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  reparameterization_[idx] += message;

  return offs;
}

/**************************************/

OneOfNTRWSFactor::OneOfNTRWSFactor(const Storage1D<TRWSVar*>& involved_vars) :
  TRWSFactor(involved_vars) {}

/*virtual*/ double OneOfNTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

  message.resize(2);

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else
      idx = k;
  }

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

  reparameterization_[idx] += message;

  return offs;
}

/**************************************/

CardinalityTRWSFactor::CardinalityTRWSFactor(const Storage1D<TRWSVar*>& involved_vars, const Math1D::Vector<float>& cost) :
  TRWSFactor(involved_vars), cost_(cost) {}

/*virtual*/ double CardinalityTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {

  assert(var->nLabels() == 2);

  message.resize(2);
  message.set_constant(1e300);

  const uint nVars = involved_var_.size();

  //std::cerr << "nVars: " << nVars << std::endl;

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else
      idx = k;
  }

  assert(idx < nVars);

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

  //DEBUG
  // if (involved_var_.size() == 2 && involved_var_[0]->rank() == 48905 && involved_var_[1]->rank() == 48906) {

  //   std::cerr << "CARD, idx: " << idx << std::endl;
  //   std::cerr << "cost: " << cost_ << std::endl;

  //   for (uint i=0; i < 2; i++) {
  //     std::cerr << "reparameterization [" << i << "]: " << reparameterization_[i] << std::endl;
  //   }

  //   for (uint i=0; i < 2; i++) {
  //     std::cerr << "params[" << i << "]: " << param[i] << std::endl;
  //   }
    
  //   std::cerr << "computed message: " << message << std::endl;
  // }
  //END_DEBUG

  double total_offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  reparameterization_[idx] += message;

  //std::cerr << "end card. message comp." << std::endl; 

  return total_offs;
}

/**************************************/

BILPTRWSFactor::BILPTRWSFactor(const Storage1D<TRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                               int rhs_lower, int rhs_upper)
  : TRWSFactor(involved_vars), positive_(positive), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

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

/*virtual*/ double BILPTRWSFactor::compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) {


  message.resize(2,0.0);

  const uint nVars = involved_var_.size();

  //std::cerr << "nVars: " << nVars << std::endl;

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else
      idx = k;
  }

  assert(idx < nVars);

  /**** forward ****/

  Math3D::Tensor<double> forward(range_,2,idx+1);
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
  if (int(zero_offset_)+end_mul >= 0
      && int(zero_offset_)+end_mul < range_) {
    backward_light(zero_offset_ + end_mul,last_var) = -param[last_var][1];
  }

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
  }

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  reparameterization_[idx] += message;

  return offs;
}


/**************************************/

FactorTRWS::FactorTRWS(uint nVars, uint nFactors) :  var_(nVars), factor_(nFactors), rank2var_(nVars,MAX_UINT), labeling_(nVars,0), 
                                                     nUsedVars_(0), nUsedFactors_(0), constant_energy_(0.0) {}

FactorTRWS::~FactorTRWS() {

  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];

  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

void FactorTRWS::add_var(const Math1D::Vector<double>& cost) {

  //std::cerr << "adding var" << std::endl;

  assert(nUsedVars_ < var_.size());
  var_[nUsedVars_] = new TRWSVar(cost,nUsedVars_);
  rank2var_[nUsedVars_] = nUsedVars_;

  nUsedVars_++;

  //std::cerr << "added" << std::endl;
}
  
void FactorTRWS::add_var(const Math1D::Vector<float>& cost) {

  Math1D::Vector<double> dcost(cost.size());
  for (uint k=0; k < cost.size(); k++)
    dcost[k] = cost[k];

  add_var(dcost);
}

void FactorTRWS::add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost) {

  Storage1D<TRWSVar*> vars(2);
  vars[0] = var_[var1];
  vars[1] = var_[var2];

  assert(var1 < nUsedVars_);
  assert(var2 < nUsedVars_);
  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new BinaryTRWSFactor(vars,cost);

  nUsedFactors_++;
}

void FactorTRWS::add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {

  Storage1D<TRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  assert(nUsedFactors_ < factor_.size());

  if (!ref)
    factor_[nUsedFactors_] = new TernaryTRWSFactor(vars,cost);
  else
    factor_[nUsedFactors_] = new TernaryTRWSRefFactor(vars,cost);

  nUsedFactors_++;
}

void FactorTRWS::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  Storage1D<TRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new SecondDiffTRWSFactor(vars,lambda);

  nUsedFactors_++;
}

void FactorTRWS::add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                         const Storage1D<Math3D::Tensor<float> >& cost) {

  Storage1D<TRWSVar*> vars(4);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];
  vars[3] = var_[var4];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new FourthOrderTRWSFactor(vars,cost);

  nUsedFactors_++;
}

void FactorTRWS::add_two_level_factor(Math1D::Vector<uint>& var, Math2D::Matrix<uint>& exceptions,
                                      double exception_cost, double standard_cost) {

  Storage1D<TRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new TwoLevelTRWSFactor(vars,exceptions,exception_cost,standard_cost);

  nUsedFactors_++;
}

void FactorTRWS::add_one_of_n_factor(Math1D::Vector<uint>& var) {

  Storage1D<TRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new OneOfNTRWSFactor(vars);

  nUsedFactors_++;
}

void FactorTRWS::add_cardinality_factor(Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost)  {

  if (var.size() == 1) {
    Math1D::Vector<double> temp(cost.size());
    for (uint k=0; k < cost.size(); k++)
      temp[k] = cost[k];
    var_[var[0]]->add_cost(temp);
  }
  else {

    Storage1D<TRWSVar*> vars(var.size());
  
    for (uint k=0; k < var.size(); k++)
      vars[k] = var_[var[k]];
    
    assert(nUsedFactors_ < factor_.size());
    
    factor_[nUsedFactors_] = new CardinalityTRWSFactor(vars,cost);
    
    nUsedFactors_++;  
  }
}

void FactorTRWS::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
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
    
    Storage1D<TRWSVar*> vars(nUseful);
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
    
    factor_[nUsedFactors_] = new BILPTRWSFactor(vars,reduced_positive,rhs_lower,rhs_upper);
    
    nUsedFactors_++;  
  }
  else {
    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;
  }
}

double FactorTRWS::optimize(uint nIter) {

  std::cerr << "efficient TS scheme" << std::endl;

  for (uint f=0; f < nUsedFactors_; f++) {
    factor_[f]->compute_rank_range();
  }

  for (uint v=0; v < nUsedVars_; v++) {
    var_[v]->set_up_chains();
  }

  std::cerr.precision(8);

  uint arg_min;

  double bound = -1e300;

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "******* iteration " << iter << " *********" << std::endl;

    //forward
    double forward_lower = 0.0;

    for (uint i=0; i < nUsedVars_; i++) {

      // if ((i%1250) == 0)
      // 	std::cerr << "i: " << i << "/" << nUsedVars_ << std::endl;

      //std::cerr << "A" << std::endl;
      //std::cerr << "associated var: " << rank2var_[i] << std::endl;

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors

      TRWSVar* cur_var = var_[rank2var_[i]];

      //average
      constant_energy_ += cur_var->average_forward(forward_lower,arg_min);

      labeling_[rank2var_[i]] = arg_min;
    }

    forward_lower += constant_energy_;
    //std::cerr << "forward bound: " << forward_lower << std::endl;
    //std::cerr << "constant: " << constant_energy_ << std::endl;


    //backward
    double backward_lower = 0.0;

    for (int i=nUsedVars_-1; i >= 0; i--) {

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors

      TRWSVar* cur_var = var_[rank2var_[i]];

      //average
      constant_energy_ += cur_var->average_backward(backward_lower,arg_min);

      labeling_[rank2var_[i]] = arg_min;
    }
    backward_lower += constant_energy_;
    bound = backward_lower;

    std::cerr << "backward bound: " << backward_lower << std::endl;
  }

  return bound;
}

const Math1D::Vector<uint>& FactorTRWS::labeling() {
  return labeling_;
}

/******************************************************************************/

NaiveTRWSVar::NaiveTRWSVar(const Math1D::Vector<float>& cost, uint rank) :
  cost_(cost), rank_(rank) {}

void NaiveTRWSVar::add_cost(const Math1D::Vector<float>& add_cost) {
  cost_ += add_cost;
}

void NaiveTRWSVar::add_factor(NaiveTRWSFactor* factor) {
  
  uint size = adjacent_factor_.size();
  
  adjacent_factor_.resize(size+1);
  adjacent_factor_[size] = factor;
}

const Math1D::Vector<float>& NaiveTRWSVar::cost() const {

  return cost_;
}

double NaiveTRWSVar::average(uint& arg_min) {
  
  Math1D::Vector<double> sum(cost_.size(),0.0);
  
  double offs = 1e300;

  for (uint c=0; c < chain_.size(); c++) {
    sum += chain_[c].var_parameters_;
  }
  
  //std::cerr << "sum: " << sum << std::endl;

  
  for (uint k=0; k < sum.size(); k++) {
    if (sum[k] < offs) {
      offs = sum[k];
      arg_min = k;
    }
  }

  for (uint k=0; k < sum.size(); k++) {
    sum[k] = (sum[k] - offs) / chain_.size();
  }    

  for (uint c=0; c < chain_.size(); c++)
    chain_[c].var_parameters_ = sum;

  return offs;
}

double NaiveTRWSVar::reparameterize_forward(bool relaxed_monotonicity) {

  Math1D::Vector<double> message;

  double offs = 0.0;

  //std::cerr << chain_.size() << " chains" << std::endl;

  assert(chain_.size() >= 1);

  for (uint c=0; c < chain_.size(); c++) {

    FactorChainElement& cur_element = chain_[c];

    //std::cerr << "before message: " << cur_element.var_parameters_ << std::endl;

    NaiveTRWSFactor* incoming = cur_element.incoming_;
    NaiveTRWSFactor* outgoing = cur_element.outgoing_;

    if (relaxed_monotonicity) {

      assert(incoming->prev_var() != this);

      if (incoming->prev_var() != 0 && 
          ((incoming->prev_var()->rank() < incoming->prev_factor()->max_rank())
           || (incoming->prev_var()->rank() > rank_)) ) {

        bool is_first = true;

        const Storage1D<NaiveTRWSVar*> involved_vars = incoming->involved_vars();

        for (uint i=0; i < involved_vars.size(); i++) {
          uint rank = involved_vars[i]->rank();

          if (rank < rank_ && rank != incoming->prev_var()->rank()) {
            is_first = false;
            break;
          }
        }

        if (is_first) {

          double cur_offs = incoming->prev_factor()->
            compute_message_and_reparameterize(incoming->prev_var(),message,true);
          incoming->prev_var()->add_factor_params(incoming->prev_factor(),message);

          //DEBUG
          //incoming->prev_factor()->set_forward(cur_offs);
          //END_DEBUG	  

          offs += cur_offs;
        }
      }

      if (outgoing != 0 && rank_ > outgoing->min_rank()) {
        outgoing->compute_message_and_reparameterize(this,message,true);
        cur_element.var_parameters_ += message;
      }
    }

    if (!relaxed_monotonicity) {

      if ( incoming->min_rank() < rank_) {

        double cur_offs = incoming->
          compute_message_and_reparameterize(this,message,true);
        cur_element.var_parameters_ += message;
	
        if (incoming->max_rank() == rank_) {
          offs += cur_offs;

          //DEBUG
          //incoming->set_forward(cur_offs);
          //END_DEBUG
        }
      }
      else {

        //this should only happen if the current var is a dead-end node of the factor
	
        if (cur_element.outgoing_ != 0) {
	  
          // std::cerr << "outgoing: " << cur_element.outgoing_ << std::endl;
          // std::cerr << "in min-rank: " << cur_element.incoming_->min_rank() << std::endl;
          // std::cerr << "out min-rank: " << cur_element.outgoing_->min_rank() << std::endl;
          // std::cerr << "rank: " << rank_ << std::endl;
          // if (cur_element.incoming_->prev_var() != 0)
          //   std::cerr << "entering var has rank " << cur_element.incoming_->prev_var()->rank() << std::endl;
          // else
          //   std::cerr << "no entering var" << std::endl;
        }
	
        assert(relaxed_monotonicity || cur_element.outgoing_ == 0);
      }
    }
    else {
      //relaxed monotonicity 

      if (incoming->min_rank() == rank_) {

        // if (incoming->min_rank() == 1717 && incoming->max_rank() == 1720) {
        //   std::cerr << "address " << incoming << ", case B, not storing, rank: " << rank_ << std::endl; 
        // }

        if (outgoing != 0) {
          //variable is the outgoing variable of the first factor in a chain
          // => the variable has been updated in the last backward pass, no need for a message now

          assert(incoming->prev_var() == 0);
          //there will be a later message to the variable, only then will we have the correct forward term
        }
        else {
          //dead end

          if (incoming->prev_var() != 0 && incoming->prev_var()->rank() > rank_) {
            //in this case, a message HAS TO BE PASSED
            // since above a message was passed to incoming->prev_var()

            incoming->compute_message_and_reparameterize(this,message,true);
            cur_element.var_parameters_ += message;
          }
          else {

            //a dead-end of minimal rank that does not fit above has to be the start of a chain
            // => no message needs to be passed

            assert(incoming->prev_var() == 0);
          }
        }
      }
      else {
        //variable is not of minimal rank

        if (outgoing == 0) {
          //the present variable is not a linking variable in the chain

          //here, always a message needs to be passed since there are variables of lower rank in the factor
          double cur_offs = incoming->compute_message_and_reparameterize(this,message,true);
          cur_element.var_parameters_ += message;

          if (incoming->next_var() == 0) {

            //NOTE: sometimes the chain-end is not the highest variable in a factor 
            //  (only if the highest var. is the incoming variable)
            // but then, we have to wait until we reached the highest variable before we can determine
            // the forward share of the factor

            if ( incoming->max_rank() == rank_  ) {
              //end of chain

              assert(incoming->prev_var() != this);

              offs += cur_offs;
	      
              //DEBUG
              //incoming->set_forward(cur_offs);
              //END_DEBUG
            }
          }
        }
        else {
          //linking variable
	  
          if (incoming->max_rank() == rank_ && outgoing->min_rank() == rank_) {

            double cur_offs = incoming->compute_message_and_reparameterize(this,message,true);
            cur_element.var_parameters_ += message;

            offs += cur_offs;

            //DEBUG
            //incoming->set_forward(cur_offs);
            //END_DEBUG
          }
          else {

            if (outgoing->min_rank() != rank_) {
              //in this case the forward message was already passed when visiting the minimal
              // variable of the outgoing factor
            }
            else {	      
              assert(incoming->max_rank() != rank_);

              incoming->compute_message_and_reparameterize(this,message,true);
              cur_element.var_parameters_ += message;
            }

            if (outgoing->next_var() == 0 && outgoing->max_rank() == rank_) {

              //non-standard chain end in the next factor
              //will need to handle this, but it can be done only AFTER averaging
            }
          }
        }
      }
    }

    //std::cerr << "after message: " << cur_element.var_parameters_ << std::endl;

    //outgoing can be ignored (no changes since last visit)
  }

  return offs;
}

double NaiveTRWSVar::special_fwd_cases() {

  Math1D::Vector<double> message;

  double offs = 0.0;

  for (uint c=0; c < chain_.size(); c++) {

    FactorChainElement& cur_element = chain_[c];

    //std::cerr << "before message: " << cur_element.var_parameters_ << std::endl;

    //NaiveTRWSFactor* incoming = cur_element.incoming_;
    NaiveTRWSFactor* outgoing = cur_element.outgoing_;

    if (outgoing != 0 && outgoing->next_var() == 0 && outgoing->max_rank() == rank_) {
      
      NaiveTRWSVar* out_var = 0;
      for (uint i=0; i < outgoing->involved_vars().size(); i++) {
        if (outgoing->involved_vars()[i] != this && 
            (out_var == 0 || out_var->rank() < outgoing->involved_vars()[i]->rank() ) ) {
          out_var = outgoing->involved_vars()[i];
        }
      }
      
      assert(out_var != 0);

      double new_offs = outgoing->compute_message_and_reparameterize(out_var,message,false);		

      //DEBUG
      Math1D::Vector<double> org_msg = message;
      //END_DEBUG

      message += out_var->factor_params(outgoing);
      
      double addon = new_offs + message.min();
      offs += addon;	      
    }
  }

  return offs;
}
  
double NaiveTRWSVar::reparameterize_backward(bool relaxed_monotonicity) {

  Math1D::Vector<double> message;

  double offs = 0.0;

  for (uint c=0; c < chain_.size(); c++) {

    FactorChainElement& cur_element = chain_[c];

    NaiveTRWSFactor* outgoing = cur_element.outgoing_;
    NaiveTRWSFactor* incoming = cur_element.incoming_;

    if (outgoing != 0) {

      if (relaxed_monotonicity) {

        if (outgoing->next_var() != 0 
            && outgoing->next_var()->rank() >  outgoing->next_factor()->min_rank()) {
	  
          assert(rank_ < outgoing->next_var()->rank());
	  
          bool is_first = true;
	  
          const Storage1D<NaiveTRWSVar*> involved_vars = outgoing->involved_vars();
	  
          for (uint i=0; i < involved_vars.size(); i++) {
	    
            uint rank = involved_vars[i]->rank();
	    
            if (rank > rank_ && rank != outgoing->next_var()->rank()) {
              is_first = false;
              break;
            }
          }
	  
          if (is_first) {
            double addon = outgoing->next_factor()->
              compute_message_and_reparameterize(outgoing->next_var(),message,true);
            outgoing->next_var()->add_factor_params(outgoing->next_factor(),message);
            offs += addon;
          }
        }

        if (rank_ < incoming->max_rank()) {
          incoming->compute_message_and_reparameterize(this,message,true);
          cur_element.var_parameters_ += message;
        }
      }

      //this needs to be done in next to all cases
      // possible exception: outgoing is a chain end, where the incoming variable is of highest rank

      if (!relaxed_monotonicity) {

        double cur_offs = outgoing->
          compute_message_and_reparameterize(this,message,true);
        cur_element.var_parameters_ += message;

        offs += cur_offs;
      }
      else {

        if (incoming->max_rank() == rank_) {
          //otherwise the desired message was already sent when processing 
          // the highest ranked var. of the incoming factor

          double cur_offs = outgoing->
            compute_message_and_reparameterize(this,message,true);
          cur_element.var_parameters_ += message;
	 
          if (outgoing->min_rank() == rank_) {
            offs += cur_offs;
          }
        }
      }
    }
    else {
      //outgoing == 0

      if (relaxed_monotonicity) {
	
        if (incoming->next_var() != 0 
            && (incoming->next_var()->rank() > incoming->next_factor()->min_rank()
                || rank_ > incoming->next_var()->rank())) { 

          bool is_first = true;

          const Storage1D<NaiveTRWSVar*> involved_vars = incoming->involved_vars();

          for (uint i=0; i < involved_vars.size(); i++) {
	    
            uint rank = involved_vars[i]->rank();

            if (rank > rank_ && rank != incoming->next_var()->rank()) {
              is_first = false;
              break;
            }
          }

          if (is_first) {

            double addon = incoming->next_factor()->
              compute_message_and_reparameterize(incoming->next_var(),message,true);
            incoming->next_var()->add_factor_params(incoming->next_factor(),message);	

            offs += addon;
          }
        }
      }

      double cur_offs = incoming->
        compute_message_and_reparameterize(this,message,true);
      cur_element.var_parameters_ += message;
      
      if (incoming->max_rank() > rank_) {

      	if (incoming->min_rank() == rank_)
      	  offs += cur_offs;
      }
    }
  }

  return offs;
}

uint NaiveTRWSVar::nChains() const {

  return chain_.size();
}

void NaiveTRWSVar::set_up_chains() {

  for (uint f=0; f < adjacent_factor_.size(); f++) {

    if (adjacent_factor_[f]->prev_factor() != 0)
      assert(adjacent_factor_[f]->prev_factor()->next_factor() == adjacent_factor_[f]);

    if (adjacent_factor_[f]->next_factor() != 0)
      assert(adjacent_factor_[f]->next_factor()->prev_factor() == adjacent_factor_[f]);

    uint size = chain_.size();

    if (adjacent_factor_[f]->prev_var() != this && adjacent_factor_[f]->next_var() != this) {
      
      //the present var. is a dead-end node in this chain
      chain_.resize(size+1);
      chain_[size].incoming_ = adjacent_factor_[f];
      chain_[size].outgoing_ = 0;
      size++;
    }
    else if (adjacent_factor_[f]->prev_var() == this) {

      assert (adjacent_factor_[f]->next_var() != this);

      if (adjacent_factor_[f]->prev_factor() == 0) {
        chain_.resize(size+1);
        chain_[size].incoming_ = adjacent_factor_[f];
        chain_[size].outgoing_ = 0;
        size++;
      }
      //otherwise covered below
    }
    else {

      assert(adjacent_factor_[f]->next_var() == this);
      chain_.resize(size+1);

      chain_[size].incoming_ = adjacent_factor_[f];
      chain_[size].outgoing_ = adjacent_factor_[f]->next_factor();
      size++;

      //DEBUG
      uint max_prev = 0;
      const Storage1D<NaiveTRWSVar*>& first_vars = adjacent_factor_[f]->involved_vars();

      for (uint k=0; k < first_vars.size(); k++) {
        if (first_vars[k] != this)
          max_prev = std::max(max_prev,first_vars[k]->rank());
      }

      uint min_next = MAX_UINT;
      const Storage1D<NaiveTRWSVar*>& second_vars = adjacent_factor_[f]->next_factor()->involved_vars();

      for (uint k=0; k < second_vars.size(); k++) {
        if (second_vars[k] != this)
          min_next = std::min(min_next,second_vars[k]->rank());
      }
      // if (! (max_prev < min_next) ) {
      // 	std::cerr << "max_prev: " << max_prev << std::endl;
      // 	std::cerr << "min_next: " << min_next << std::endl;

      // 	std::cerr << "linking var: " << rank_ << std::endl;
      // 	std::cerr << "prev factor has range " << adjacent_factor_[f]->min_rank() 
      // 		  << " - " << adjacent_factor_[f]->max_rank() << std::endl;
      // 	std::cerr << "next factor has range " << adjacent_factor_[f]->next_factor()->min_rank() 
      // 		  << " - " << adjacent_factor_[f]->next_factor()->max_rank() << std::endl;
      // }
      // assert(max_prev < min_next);
      //END_DEBUG

    }
  }

  for (uint c=0; c < chain_.size(); c++) {

    chain_[c].var_parameters_.resize(cost_.size());
    for (uint k=0; k < cost_.size(); k++)
      chain_[c].var_parameters_[k] = cost_[k] / chain_.size();
  }

}

uint NaiveTRWSVar::rank() const {
  return rank_;
}

void NaiveTRWSVar::set_rank(uint new_rank) {
  rank_ = new_rank;
}

size_t NaiveTRWSVar::nLabels() const {
  return cost_.size();
}

const Storage1D<NaiveTRWSFactor*>& NaiveTRWSVar::adjacent_factor() const {
  return adjacent_factor_;
}

const Math1D::Vector<double>& NaiveTRWSVar::factor_params(const NaiveTRWSFactor* factor) const {

  for (uint k=0; k < chain_.size(); k++) {

    if (chain_[k].incoming_ == factor || chain_[k].outgoing_ == factor)
      return chain_[k].var_parameters_;
  }

  assert(false);
  return chain_[0].var_parameters_;
}

double NaiveTRWSVar::add_factor_params(NaiveTRWSFactor* factor, const Math1D::Vector<double>& message) {

  double offs = 0.0;

  //DEBUG
  bool found = false;
  //END_DEBUG

  for (uint k=0; k < chain_.size(); k++) {

    if (chain_[k].incoming_ == factor || chain_[k].outgoing_ == factor) {
      chain_[k].var_parameters_ += message;

      //DEBUG
      found = true;
      //END_DEBUG

      break;
    }
  }

  //DEBUG
  assert(found);
  //END_DEBUG

  return offs;
}

double NaiveTRWSVar::normalize_factor_params(NaiveTRWSFactor* factor) {

  double offs = 0.0;

  for (uint k=0; k < chain_.size(); k++) {

    if (chain_[k].incoming_ == factor || chain_[k].outgoing_ == factor) {

      offs = chain_[k].var_parameters_.min();
      for (uint l=0; l < chain_[k].var_parameters_.size(); l++)
        chain_[k].var_parameters_[l] -= offs;
    }
  }

  return offs;
}

void NaiveTRWSVar::set_factor_params(NaiveTRWSFactor* factor, const Math1D::Vector<double>& params) {

  for (uint k=0; k < chain_.size(); k++) {

    if (chain_[k].incoming_ == factor || chain_[k].outgoing_ == factor) {
      chain_[k].var_parameters_ = params;
      return;
    }
  }

  assert(false);
}

double NaiveTRWSVar::dual_value(uint& argmin) const {

  Math1D::Vector<double> param(cost_.size(),0.0);
  for (uint k=0; k < cost_.size(); k++)
    param[k] = cost_[k];

  for (uint f = 0; f < adjacent_factor_.size(); f++) {

    param += adjacent_factor_[f]->reparameterization(this);
  }

  double best = 1e300;

  //std::cerr << "param: " << param << std::endl;

  for (uint k=0; k < cost_.size(); k++)  {

    if (param[k] < best) {
      best = param[k];
      argmin = k;
    }
  }

  return best;
}


/*******************/

NaiveTRWSFactor::NaiveTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars) :
  involved_var_(involved_vars), prev_var_(0), prev_factor_(0), next_var_(0), next_factor_(0) {

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

//DEBUG
// void NaiveTRWSFactor::set_forward(double val) {

//   // if (min_rank_ == 3 && max_rank_ == 827) {
//   //   std::cerr << "factor " << this <<  ": request to set forward to " << val << std::endl;
//   // }

//   if (forward_ != 0.0) {
//     std::cerr << "forward: " << forward_ << ", address: " << this << std::endl;
//     std::cerr << "rank-range: " << min_rank() << " - " << max_rank() << std::endl;
//     if (prev_var() != 0) {
//       std::cerr << "in-var: " << prev_var()->rank() << std::endl;
//     }
//     else
//       std::cerr << "no in-var" << std::endl;
//     std::cerr << "out-var: " << next_var() << std::endl;
//   }

//   assert(forward_ == 0.0);
//   forward_ = val;
// }
//END_DEBUG

uint NaiveTRWSFactor::min_rank() const {
  return min_rank_;
}

uint NaiveTRWSFactor::max_rank() const {
  return max_rank_;
}

void NaiveTRWSFactor::compute_rank_range() {

  min_rank_ = MAX_UINT;
  max_rank_ = 0;

  uint nVars = involved_var_.size();

  for (uint v=0; v < nVars; v++) {

    min_rank_ = std::min(min_rank_,involved_var_[v]->rank());
    max_rank_ = std::max(max_rank_,involved_var_[v]->rank());
  }
}

NaiveTRWSVar* NaiveTRWSFactor::prev_var() const {
  return prev_var_;
}

NaiveTRWSVar* NaiveTRWSFactor::next_var() const {
  return next_var_;
}

NaiveTRWSFactor* NaiveTRWSFactor::prev_factor() const {
  return prev_factor_;
}

NaiveTRWSFactor* NaiveTRWSFactor::next_factor() const {
  return next_factor_;
}

void NaiveTRWSFactor::set_prev_var(NaiveTRWSVar* var) {
  prev_var_ = var;
}

void NaiveTRWSFactor::set_next_var(NaiveTRWSVar* var) {
  next_var_ = var;
}

void NaiveTRWSFactor::set_prev_factor(NaiveTRWSFactor* factor) {
  prev_factor_ = factor;
}

void NaiveTRWSFactor::set_next_factor(NaiveTRWSFactor* factor) {
  next_factor_ = factor;
}

const Storage1D<NaiveTRWSVar*>& NaiveTRWSFactor::involved_vars() const {
  return involved_var_;
}

const Math1D::Vector<double>& NaiveTRWSFactor::reparameterization(const NaiveTRWSVar* var) const {

  for (uint v=0; v < involved_var_.size(); v++) {

    if (involved_var_[v] == var)
      return reparameterization_[v];
  }

  assert(false);
  return reparameterization_[0];
}

Math1D::Vector<double>& NaiveTRWSFactor::reparameterization(const NaiveTRWSVar* var)  {

  for (uint v=0; v < involved_var_.size(); v++) {

    if (involved_var_[v] == var)
      return reparameterization_[v];
  }

  assert(false);
  return reparameterization_[0];
}



/*******************/

NaiveBinaryTRWSFactor::NaiveBinaryTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                                             const Math2D::Matrix<float>& cost) :
  NaiveTRWSFactor(involved_vars), cost_(cost) {
}

/*virtual*/
double NaiveBinaryTRWSFactor::cost(const Math1D::Vector<uint>& labeling) const {

  double cost = cost_(labeling[0],labeling[1]);
  cost -= reparameterization_[0][labeling[0]];
  cost -= reparameterization_[1][labeling[1]];

  cost += involved_var_[0]->factor_params(this)[labeling[0]];
  cost += involved_var_[1]->factor_params(this)[labeling[1]];

  return cost;
}

/*virtual*/ 
double NaiveBinaryTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, 
                                                                 Math1D::Vector<double>& message,
                                                                 bool reparameterize) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  //this routine also updates reparameterization_
  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++) {
    param[k] -= involved_var_[k]->factor_params(this);
  }

  uint idx = 0;

  if (var == involved_var_[0]) {

    idx = 0;

    message.resize(nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double hyp = cost_(l1,l2) - param[1][l2];

        if (hyp < best)
          best = hyp;
      }

      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else {
    assert(var == involved_var_[1]);

    idx = 1;

    message.resize(nLabels2);

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      double best = 1e300;

      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double hyp = cost_(l1,l2) - param[0][l1];

        if (hyp < best)
          best = hyp;
      }

      message[l2] = best - reparameterization_[1][l2];
    }
  }

  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  if (reparameterize)
    reparameterization_[idx] += message;
  
  return offs;
}


//used for the subgradient method
/*virtual*/
void NaiveBinaryTRWSFactor::compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) {


  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++) {
    param[k] -= involved_var_[k]->factor_params(this);
  }

  if (out_var == involved_var_[0]) {

    forward.resize(nLabels1);
    trace.resize(2,nLabels1);
    
    for (uint l1 = 0; l1 < nLabels1; l1++) {

      double best = 1e300;
      uint arg_best = MAX_UINT;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double hyp = cost_(l1,l2) - param[1][l2];

        if (hyp < best) {
          best = hyp;
          arg_best = l2;
        }
      }

      forward[l1] = best - reparameterization_[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = arg_best;
    }
  }
  else {
    assert(out_var == involved_var_[1]);

    forward.resize(nLabels2);
    trace.resize(2,nLabels2);

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      double best = 1e300;
      uint arg_best = MAX_UINT;

      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double hyp = cost_(l1,l2) - param[0][l1];

        if (hyp < best) {
          best = hyp;
          arg_best = l1;
        }
      }

      forward[l2] = best - reparameterization_[1][l2];
      trace(0,l2) = arg_best;
      trace(1,l2) = l2;
    }
  }

  //std::cerr << "cost: " << cost_ << std::endl;
  //std::cerr << "new forward: " << forward << std::endl;
}


/*******************/

NaiveTernaryTRWSFactorBase::NaiveTernaryTRWSFactorBase(const Storage1D<NaiveTRWSVar*>& involved_vars) :
  NaiveTRWSFactor(involved_vars)
{
}


  
double NaiveTernaryTRWSFactorBase::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                                      const Math3D::Tensor<float>& cost, bool reparameterize) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  //this routine also updates reparameterization_
  uint idx = 0;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->factor_params(this);
  
  if (var == involved_var_[0]) {

    message.resize(nLabels1);
    
    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - param[1][l2] - param[2][l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else if (var == involved_var_[1]) {

    message.resize(nLabels2);
    
    idx = 1;

    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - param[0][l1] - param[2][l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l2] = best - reparameterization_[1][l2];
    }    
  }
  else {
    assert(var == involved_var_[2]);
    
    message.resize(nLabels3);

    idx = 2;

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
          double hyp = cost(l1,l2,l3) - param[0][l1] - param[1][l2];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l3] = best - reparameterization_[2][l3];
    }
  }

  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  if (reparameterize)
    reparameterization_[idx] += message;
  
  return offs;  
}

void NaiveTernaryTRWSFactorBase::compute_forward_base(const Math3D::Tensor<float>& cost, 
                                                      const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, 
                                                      Math2D::Matrix<uint>& trace) {


  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->factor_params(this);
  
  if (out_var == involved_var_[0]) {

    forward.resize(nLabels1);
    trace.resize(3,nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      uint argbest2 = MAX_UINT;
      uint argbest3 = MAX_UINT;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - param[1][l2] - param[2][l3];
	  
          if (hyp < best) {
            best = hyp;
            argbest2 = l2;
            argbest3 = l3;
          }
        }
      }
      
      forward[l1] = best - reparameterization_[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = argbest2;
      trace(2,l1) = argbest3;
    }
  }
  else if (out_var == involved_var_[1]) {

    forward.resize(nLabels2);
    trace.resize(3,nLabels2);

    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest3 = MAX_UINT;

      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - param[0][l1] - param[2][l3];
	  
          if (hyp < best) {
            best = hyp;
            argbest1 = l1;
            argbest3 = l3;
          }
        }
      }
      
      forward[l2] = best - reparameterization_[1][l2];
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
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
          double hyp = cost(l1,l2,l3) - param[0][l1] - param[1][l2];
	  
          if (hyp < best) {
            best = hyp;
            argbest1 = l1;
            argbest2 = l2;
          }
        }
      }
      
      forward[l3] = best - reparameterization_[2][l3];
      trace(0,l3) = argbest1;
      trace(1,l3) = argbest2;
      trace(2,l3) = l3;
    }
  }

}

/*******************/

NaiveTernaryTRWSFactor::NaiveTernaryTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                                               const Math3D::Tensor<float>& cost) :
  NaiveTernaryTRWSFactorBase(involved_vars), cost_(cost) 
{
}
 
/*virtual*/ 
double NaiveTernaryTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, 
                                                                  Math1D::Vector<double>& message,
                                                                  bool reparameterize) {

  return NaiveTernaryTRWSFactorBase::compute_message_and_reparameterize(var, message, cost_, reparameterize);
}

/*virtual*/ 
void NaiveTernaryTRWSFactor::compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) {

  NaiveTernaryTRWSFactorBase::compute_forward_base(cost_,out_var,forward,trace);
}

/*virtual*/
double NaiveTernaryTRWSFactor::cost(const Math1D::Vector<uint>& labeling) const {

  double cost = cost_(labeling[0],labeling[1],labeling[2]);
  cost -= reparameterization_[0][labeling[0]];
  cost -= reparameterization_[1][labeling[1]];
  cost -= reparameterization_[2][labeling[2]];

  cost += involved_var_[0]->factor_params(this)[labeling[0]];
  cost += involved_var_[1]->factor_params(this)[labeling[1]];
  cost += involved_var_[2]->factor_params(this)[labeling[2]];

  return cost;
}



/*******************/

NaiveTernaryTRWSRefFactor::NaiveTernaryTRWSRefFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                                                     const Math3D::Tensor<float>& cost) :
  NaiveTernaryTRWSFactorBase(involved_vars), cost_(cost) 
{
}

/*virtual */
double NaiveTernaryTRWSRefFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                                     bool reparameterize) {

  return NaiveTernaryTRWSFactorBase::compute_message_and_reparameterize(var, message, cost_, reparameterize);
}

/*virtual*/ 
void NaiveTernaryTRWSRefFactor::compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) {

  NaiveTernaryTRWSFactorBase::compute_forward_base(cost_,out_var,forward,trace);
}

/*virtual*/
double NaiveTernaryTRWSRefFactor::cost(const Math1D::Vector<uint>& labeling) const {

  double cost = cost_(labeling[0],labeling[1],labeling[2]);
  cost -= reparameterization_[0][labeling[0]];
  cost -= reparameterization_[1][labeling[1]];
  cost -= reparameterization_[2][labeling[2]];

  cost += involved_var_[0]->factor_params(this)[labeling[0]];
  cost += involved_var_[1]->factor_params(this)[labeling[1]];
  cost += involved_var_[2]->factor_params(this)[labeling[2]];

  return cost;
}

/*******************/

NaiveSecondDiffTRWSFactor::NaiveSecondDiffTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars, float lambda) :
  NaiveTRWSFactor(involved_vars), lambda_(lambda) {

  assert(involved_vars.size() == 3);
}

/*virtual */
double NaiveSecondDiffTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                                     bool reparameterize) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  //this routine also updates reparameterization_
  uint idx = 0;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->factor_params(this);
  
  if (var == involved_var_[0]) {

    idx = 0;

    message.resize(nLabels1);

    double base_cost = -param[1].max() - param[2].max() + 3*lambda_;

    for (int l1=0; l1 < int(nLabels1); l1++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[1][l2] - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[1][l2] - param[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    message.resize(nLabels2);

    double base_cost = -param[0].max() - param[2].max() + 3*lambda_;

    for (int l2=0; l2 < int(nLabels2); l2++) {
      
      double best = base_cost;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[0][l1] - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[0][l1] - param[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l2] = best - reparameterization_[1][l2];
    }    

  }
  else {
    assert(var == involved_var_[2]);

    idx = 2;

    message.resize(nLabels3);

    double base_cost = -param[0].max() - param[1].max() + 3*lambda_;

    for (int l3=0; l3 < int(nLabels3); l3++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {
        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[0][l1] - param[1][l2];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[0][l1] - param[1][l2] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l3] = best - reparameterization_[2][l3];
    }
  }


  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  if (reparameterize)
    reparameterization_[idx] += message;
  
  return offs;  
}

//used for the subgradient method
/*virtual */
void NaiveSecondDiffTRWSFactor::compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->factor_params(this);

  if (out_var == involved_var_[0]) {

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

          int so_diff = l3 - 2*l2 + l1;
	  
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

      forward[l1] = best - reparameterization_[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = best2;
      trace(2,l1) = best3;
    }
  }
  else if (out_var == involved_var_[1]) {  

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

          int so_diff = l3 - 2*l2 + l1;
	  
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

      forward[l2] = best - reparameterization_[1][l2];
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

          int so_diff = l3 - 2*l2 + l1;
	  
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

      forward[l3] = best - reparameterization_[2][l3];
      trace(0,l3) = best1;
      trace(1,l3) = best2;
      trace(2,l3) = l3;
    }
  }
}

/*virtual */
double NaiveSecondDiffTRWSFactor::cost(const Math1D::Vector<uint>& labeling) const {

  int diff1 = labeling[1] - labeling[0];
  int diff2 = labeling[2] - labeling[1];
  
  int so_diff = diff2 - diff1;
  
  if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
    return 0.0; //no cost
  else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
    return lambda_;

  return 3*lambda_;
}


/*******************/

NaiveFourthOrderTRWSFactor::NaiveFourthOrderTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                                                       const Storage1D<Math3D::Tensor<float> >& cost) :
  NaiveTRWSFactor(involved_vars), cost_(cost) {
}

/*virtual*/ 
double NaiveFourthOrderTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, 
                                                                      Math1D::Vector<double>& message,
                                                                      bool reparameterize) {

  double offs = 0.0;
  
  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();
  
  //this routine also updates reparameterization_

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->factor_params(this);

  uint idx = 0;

  if (var == involved_var_[0]) {

    message.resize(nLabels1);
    
    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[1][l2] - param[2][l3] - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l1] = best - reparameterization_[0][l1];
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    message.resize(nLabels2);
    
    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          for (uint l4 = 0; l4 < nLabels3; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[0][l1] - param[2][l3] - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l2] = best - reparameterization_[1][l2];
    }    
  }
  else if (var == involved_var_[2]) {
    
    message.resize(nLabels3);

    idx = 2;

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[0][l1] - param[1][l2] - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l3] = best - reparameterization_[2][l3];
    }
  }
  else {
    
    assert(var == involved_var_[3]);

    message.resize(nLabels4);

    idx = 3;

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[0][l1] - param[1][l2] - param[2][l3];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l4] = best - reparameterization_[3][l4];
    }    
  }

  double msg_offs = message.min();
  offs += msg_offs;

  for (uint l=0; l < message.size(); l++)
    message[l] -= msg_offs;

  if (reparameterize)
    reparameterization_[idx] += message;
  
  return offs;
}

/*virtual*/
void NaiveFourthOrderTRWSFactor::compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->factor_params(this);

  if (out_var == involved_var_[0]) {

    forward.resize(nLabels1);
    trace.resize(4,nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;
      uint argbest2 = MAX_UINT;
      uint argbest3 = MAX_UINT;
      uint argbest4 = MAX_UINT;

      for (uint l2 = 0; l2 < nLabels2; l2++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[1][l2] - param[2][l3] - param[3][l4];
	  
            if (hyp < best) {
              best = hyp;
              argbest2 = l2;
              argbest3 = l3;
              argbest4 = l4;
            }
          }
        }
      }
      
      forward[l1] = best - reparameterization_[0][l1];
      trace(0,l1) = l1;
      trace(1,l1) = argbest2;
      trace(2,l1) = argbest3;
      trace(3,l1) = argbest4;
    }
  }
  else if (out_var == involved_var_[1]) {

    forward.resize(nLabels2);
    trace.resize(4,nLabels2);
    
    for (uint l2 = 0; l2 < nLabels1; l2++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest3 = MAX_UINT;
      uint argbest4 = MAX_UINT;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          for (uint l4 = 0; l4 < nLabels3; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[0][l1] - param[2][l3] - param[3][l4];
	  
            if (hyp < best) {
              best = hyp;
              argbest1 = l1;
              argbest3 = l3;
              argbest4 = l4;
            }
          }
        }
      }
      
      forward[l2] = best - reparameterization_[1][l2];
      trace(0,l2) = argbest1;
      trace(1,l2) = l2;
      trace(2,l2) = argbest3;
      trace(3,l2) = argbest4;
    }    
  }
  else if (out_var == involved_var_[2]) {
    
    forward.resize(nLabels2);
    trace.resize(4,nLabels2);

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest2 = MAX_UINT;
      uint argbest4 = MAX_UINT;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[0][l1] - param[1][l2] - param[3][l4];
	  
            if (hyp < best) {
              best = hyp;
              argbest1 = l1;
              argbest2 = l2;
              argbest4 = l4;
            }
          }
        }
      }
      
      forward[l3] = best - reparameterization_[2][l3];
      trace(0,l3) = argbest1;
      trace(1,l3) = argbest2;
      trace(2,l3) = l3;
      trace(3,l3) = argbest4;
    }
  }
  else {
    
    assert(out_var == involved_var_[3]);

    forward.resize(nLabels3);
    trace.resize(4,nLabels3);

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      uint argbest1 = MAX_UINT;
      uint argbest2 = MAX_UINT;
      uint argbest3 = MAX_UINT;

      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
            double hyp = cost_[l1](l2,l3,l4) 
              - param[0][l1] - param[1][l2] - param[2][l3];
	  
            if (hyp < best) {
              best = hyp;
              argbest1 = l1;
              argbest2 = l2;
              argbest3 = l3;
            }
          }
        }
      }
      
      forward[l4] = best - reparameterization_[3][l4];
      trace(0,l4) = argbest1;
      trace(1,l4) = argbest2;
      trace(2,l4) = argbest3;
      trace(3,l4) = l4;
    }    
  }

}

/*virtual*/
double NaiveFourthOrderTRWSFactor::cost(const Math1D::Vector<uint>& labeling) const {

  double cost = cost_[labeling[0]](labeling[1],labeling[2],labeling[3]);
  cost -= reparameterization_[0][labeling[0]];
  cost -= reparameterization_[1][labeling[1]];
  cost -= reparameterization_[2][labeling[2]];
  cost -= reparameterization_[3][labeling[3]];

  cost += involved_var_[0]->factor_params(this)[labeling[0]];
  cost += involved_var_[1]->factor_params(this)[labeling[1]];
  cost += involved_var_[2]->factor_params(this)[labeling[2]];
  cost += involved_var_[3]->factor_params(this)[labeling[3]];

  return cost;
}


/*******************/

NaiveTwoLevelTRWSFactor::NaiveTwoLevelTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                                                 const Math2D::Matrix<uint>& exception, double exception_cost, double std_cost) :
  NaiveTRWSFactor(involved_vars), exception_(exception), exception_cost_(exception_cost), std_cost_(std_cost) {

  const uint nVars = involved_vars.size();
  assert(nVars == exception.xDim());

  if (exception_cost > std_cost) {

    for (uint k1=0; k1 < exception_.yDim(); k1++) {
      for (uint k2=0; k2 < exception_.yDim(); k2++) {

        uint nDiffs = 0;
        for (uint l=0; l < nVars; l++)
          if (exception_(l,k1) != exception_(l,k2))
            nDiffs++;

        if (nDiffs == 1) {
	  
          INTERNAL_ERROR << "when exceptions are costlier, all exceptions need to have a hamming distance of 2 to each other. Exiting.." 
                         << std::endl;
          exit(1);
        }
      }
    }
  }

}

/*virtual*/
double NaiveTwoLevelTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                                   bool reparameterize) {


  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Math1D::Vector<double> mincost(nVars,1e300);
  Math1D::Vector<uint> argmin(nVars,MAX_UINT);

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->factor_params(this);
    else
      idx = k;

    for (uint l=0; l < involved_var_[k]->nLabels(); l++) {

      if (-param[k][l] < mincost[k]) {
        mincost[k] = -param[k][l];
        argmin[k] = l;
      }
    }
  }

  double cost_base = mincost.sum() - mincost[idx];

  message.resize(involved_var_[idx]->nLabels());

  for (uint l=0; l < involved_var_[idx]->nLabels(); l++) {

    
    if (exception_cost_ <= std_cost_) {

      double best_entry = cost_base - param[idx][l] + std_cost_;

      for (uint exc = 0; exc < exception_.yDim(); exc++) {
        if (exception_(idx,exc) == l) {
          double hyp_entry = exception_cost_;
          for (uint k=0; k < nVars; k++)
            hyp_entry -= param[k][exception_(k,idx)];
				  
          if (hyp_entry < best_entry)
            best_entry = hyp_entry;
        }
      }

      message[l] = best_entry;
    }
    else {

      double std_entry = cost_base - param[idx][l] + std_cost_;
      mincost[idx] = -param[idx][l];
      argmin[idx] = l;

      bool argmin_is_exc = false;

      double exc_entry = 1e300;

      for (uint exc = 0; exc < exception_.yDim(); exc++) {
        if (exception_(idx,exc) == l) {

          //check if the standard argmin is an exception
          if (!argmin_is_exc) {
            argmin_is_exc = true;

            for (uint k=0; k < nVars; k++) {

              if (exception_(k,exc) != argmin[k]) {
                argmin_is_exc = false;
                break;
              }
            }
          }

          double hyp_entry = exception_cost_;
          for (uint k=0; k < nVars; k++)
            hyp_entry -= param[k][exception_(k,exc)];
				  
          if (hyp_entry < exc_entry)
            exc_entry = hyp_entry;
        }
      }      

      if (argmin_is_exc) {

        double best_add = 0.0;
        double best_change = 1e300;
        uint arg_b2b = MAX_UINT;

        for (uint k=0; k < nVars; k++) {

          if (k != idx) {
            Math1D::Vector<double> temp = param[k];
            temp *= -1.0;
            std::sort(temp.direct_access(),temp.direct_access()+temp.size());

            if (temp[1] - temp[0] < best_change) {
              best_change = temp[1] - temp[0];
              best_add = temp[1];
              arg_b2b = k;
            }
          }
        }

        std_entry += best_add - mincost[arg_b2b];
      }

      message[l] = std::min(std_entry, exc_entry);
    }
  }


  double offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  if (reparameterize)
    reparameterization_[idx] += message;

  return offs;
}

/*virtual*/
void NaiveTwoLevelTRWSFactor::compute_forward(const NaiveTRWSVar* /*out_var*/, Math1D::Vector<double>& /*forward*/, 
					      Math2D::Matrix<uint>& /*trace*/) {

  TODO("compute_forward");
}


/*virtual*/
double NaiveTwoLevelTRWSFactor::cost(const Math1D::Vector<uint>& /*labeling*/) const {
  TODO("cost comp");
}

/*******************/

NaiveOneOfNTRWSFactor::NaiveOneOfNTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars) :
  NaiveTRWSFactor(involved_vars) {}

/*virtual*/
double NaiveOneOfNTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                                 bool reparameterize) {

  message.resize(2);

  // if (involved_var_.size() == 2 && involved_var_[0]->rank() == 48906 && involved_var_[1]->rank() == 48907) {
  //   std::cerr << "1ON, A, intermediate msg0: " << message[0] << std::endl;
  // }

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->factor_params(this);
    else
      idx = k;
  }

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

  if (reparameterize)
    reparameterization_[idx] += message;

  return offs;
}


/*virtual*/
void NaiveOneOfNTRWSFactor::compute_forward(const NaiveTRWSVar* /*out_var*/, Math1D::Vector<double>& /*forward*/, 
					    Math2D::Matrix<uint>& /*trace*/) {

  TODO("compute_forward");
}

/*virtual*/
double NaiveOneOfNTRWSFactor::cost(const Math1D::Vector<uint>& /*labeling*/) const {
  TODO("cost comp");
}


/*******************/

NaiveCardinalityTRWSFactor::NaiveCardinalityTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars, 
                                                       const Math1D::Vector<float>& cost) :
  NaiveTRWSFactor(involved_vars), cost_(cost) {}
  
/*virtual*/ 
double NaiveCardinalityTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                                      bool reparameterize) {


  //std::cerr << "start card. message comp." << std::endl; 

  assert(var->nLabels() == 2);

  message.resize(2);
  message.set_constant(1e300);

  const uint nVars = involved_var_.size();

  //std::cerr << "nVars: " << nVars << std::endl;

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->factor_params(this);
    else
      idx = k;
  }

  assert(idx < nVars);


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

  //DEBUG
  // if (involved_var_.size() == 2 && involved_var_[0]->rank() == 48905 && involved_var_[1]->rank() == 48906) {

  //   std::cerr << "CARD, idx: " << idx << std::endl;
  //   std::cerr << "cost: " << cost_ << std::endl;

  //   for (uint i=0; i < 2; i++) {
  //     std::cerr << "reparameterization [" << i << "]: " << reparameterization_[i] << std::endl;
  //   }

  //   for (uint i=0; i < 2; i++) {
  //     std::cerr << "params[" << i << "]: " << param[i] << std::endl;
  //   }
    
  //   std::cerr << "computed message: " << message << std::endl;
  // }
  //END_DEBUG

  double total_offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  if (reparameterize)
    reparameterization_[idx] += message;

  //std::cerr << "end card. message comp." << std::endl; 

  return total_offs;
}

/*virtual*/
void NaiveCardinalityTRWSFactor::compute_forward(const NaiveTRWSVar* /*out_var*/, Math1D::Vector<double>& /*forward*/, 
						 Math2D::Matrix<uint>& /*trace*/) {

  TODO("compute_forward");
}

/*virtual*/
double NaiveCardinalityTRWSFactor::cost(const Math1D::Vector<uint>& /*labeling*/) const {
  TODO("cost comp");
}


/*******************/

NaiveBILPTRWSFactor::NaiveBILPTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars, 
                                         const Storage1D<bool>& positive,
                                         int rhs_lower, int rhs_upper) :
  NaiveTRWSFactor(involved_vars), positive_(positive), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

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
double NaiveBILPTRWSFactor::compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                               bool reparameterize) {

  //based on [Potetz & Lee CVIU 2007]

  message.resize(2,0.0);

  const uint nVars = involved_var_.size();

  //std::cerr << "nVars: " << nVars << std::endl;

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->factor_params(this);
    else
      idx = k;
  }

  assert(idx < nVars);

  /**** forward ****/

  Math3D::Tensor<double> forward(range_,2,idx+1);
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
  if (int(zero_offset_)+end_mul >= 0
      && int(zero_offset_)+end_mul < range_) {
    backward_light(zero_offset_ + end_mul,last_var) = -param[last_var][1];
  }

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
  }

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  if (reparameterize)
    reparameterization_[idx] += message;

  return offs;
}

/*virtual*/
void NaiveBILPTRWSFactor::compute_forward(const NaiveTRWSVar* /*out_var*/, Math1D::Vector<double>& /*forward*/, 
					  Math2D::Matrix<uint>& /*trace*/) {

  TODO("compute_forward");
}

/*virtual*/
double NaiveBILPTRWSFactor::cost(const Math1D::Vector<uint>& /*labeling*/) const {
  TODO("cost comp");
}


/*******************/

NaiveFactorTRWS::NaiveFactorTRWS(uint nVars, uint nFactors) : 
  var_(nVars), factor_(nFactors), rank2var_(nVars,MAX_UINT), labeling_(nVars,0), 
  nUsedVars_(0), nUsedFactors_(0), constant_energy_(0.0) {}

NaiveFactorTRWS::~NaiveFactorTRWS() {

  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];

  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

const Math1D::Vector<uint>& NaiveFactorTRWS::labeling() {
  return labeling_;
}

void NaiveFactorTRWS::add_var(const Math1D::Vector<float>& cost) {

  //std::cerr << "adding var" << std::endl;

  assert(nUsedVars_ < var_.size());
  var_[nUsedVars_] = new NaiveTRWSVar(cost,nUsedVars_);
  rank2var_[nUsedVars_] = nUsedVars_;

  nUsedVars_++;

  //std::cerr << "added" << std::endl;
}

void NaiveFactorTRWS::add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost){

  //std::cerr << "add2" << std::endl; 

  Storage1D<NaiveTRWSVar*> vars(2);
  vars[0] = var_[var1];
  vars[1] = var_[var2];

  assert(var1 < nUsedVars_);
  assert(var2 < nUsedVars_);
  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new NaiveBinaryTRWSFactor(vars,cost);

  nUsedFactors_++;
  //std::cerr << "added2" << std::endl;
}
 
void NaiveFactorTRWS::add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {

  //std::cerr << "add3" << std::endl; 

  Storage1D<NaiveTRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  assert(nUsedFactors_ < factor_.size());

  if (!ref)
    factor_[nUsedFactors_] = new NaiveTernaryTRWSFactor(vars,cost);
  else
    factor_[nUsedFactors_] = new NaiveTernaryTRWSRefFactor(vars,cost);

  nUsedFactors_++;
}

void NaiveFactorTRWS::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  Storage1D<NaiveTRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new NaiveSecondDiffTRWSFactor(vars,lambda);

  nUsedFactors_++;
}

void NaiveFactorTRWS::add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                              const Storage1D<Math3D::Tensor<float> >& cost) {

  //std::cerr << "add4" << std::endl; 

  Storage1D<NaiveTRWSVar*> vars(4);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];
  vars[3] = var_[var4];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new NaiveFourthOrderTRWSFactor(vars,cost);

  nUsedFactors_++;
}

void NaiveFactorTRWS::add_two_level_factor(Math1D::Vector<uint>& var, Math2D::Matrix<uint>& exceptions,
                                           double exception_cost, double standard_cost) {


  Storage1D<NaiveTRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new NaiveTwoLevelTRWSFactor(vars,exceptions,exception_cost,standard_cost);

  nUsedFactors_++;
}


void NaiveFactorTRWS::add_one_of_n_factor(Math1D::Vector<uint>& var) {

  Storage1D<NaiveTRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++)
    vars[k] = var_[var[k]];

  assert(nUsedFactors_ < factor_.size());

  factor_[nUsedFactors_] = new NaiveOneOfNTRWSFactor(vars);

  nUsedFactors_++;
}

void NaiveFactorTRWS::add_cardinality_factor(Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost) {

  if (var.size() == 1) {
    var_[var[0]]->add_cost(cost);
  }
  else {

    Storage1D<NaiveTRWSVar*> vars(var.size());
  
    for (uint k=0; k < var.size(); k++)
      vars[k] = var_[var[k]];
    
    assert(nUsedFactors_ < factor_.size());
    
    factor_[nUsedFactors_] = new NaiveCardinalityTRWSFactor(vars,cost);
    
    nUsedFactors_++;  
  }
}

void NaiveFactorTRWS::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
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
    
    Storage1D<NaiveTRWSVar*> vars(nUseful);
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
    
    factor_[nUsedFactors_] = new NaiveBILPTRWSFactor(vars,reduced_positive,rhs_lower,rhs_upper);
    
    nUsedFactors_++;  
  }
  else {
    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;
  }
}


double NaiveFactorTRWS::dual_bound(bool check) {

  double dual_bound = constant_energy_;

  for (uint f=0; f < nUsedFactors_; f++) {

    if (factor_[f]->prev_var() == 0) {

      NaiveTRWSVar* in_var = 0;
      for (uint l=0; l < factor_[f]->involved_vars().size(); l++) {
        if (factor_[f]->involved_vars()[l] != factor_[f]->next_var()) {
          if (in_var == 0 || in_var->rank() > factor_[f]->involved_vars()[l]->rank()) 
            in_var = factor_[f]->involved_vars()[l];
        }
      }

      Math1D::Vector<double> vec1 = in_var->factor_params(factor_[f]);
      Math1D::Vector<double> vec2;

      NaiveTRWSFactor* cur_fac = factor_[f];

      for (uint i=0; true; i++) {

        Math1D::Vector<double>& last_forward = ((i%2) == 0) ? vec1 : vec2;
        Math1D::Vector<double>& new_forward = ((i%2) == 0) ? vec2 : vec1;

        double offs = last_forward.min();
        dual_bound += offs;
        for (uint k=0; k < last_forward.size(); k++)
          last_forward[k] -= offs;

        NaiveTRWSVar* out_var = cur_fac->next_var();

        if (cur_fac->next_var() == 0) { //end of chain

          for (uint l=0; l < cur_fac->involved_vars().size(); l++) {
            if (cur_fac->involved_vars()[l] != in_var
                && (out_var == 0 || out_var->rank() < cur_fac->involved_vars()[l]->rank()) ) {
              out_var = cur_fac->involved_vars()[l];
            }
          }
        }

        Math1D::Vector<double> temp = in_var->factor_params(cur_fac);
        in_var->set_factor_params(cur_fac,last_forward);
	
        double addon = cur_fac->compute_message_and_reparameterize(out_var,new_forward,false);
        dual_bound += addon;

        //DEBUG
        // if (check && fabs(addon + offs - cur_fac->forward_) > 0.001) {

        //   std::cerr << "----addon: " << addon << std::endl;	  
        //   std::cerr << "previous node-forward: " << offs << std::endl;
        //   std::cerr << "stored forward: " << cur_fac->forward_ << std::endl;
        //   std::cerr << "factor #" << (i+1) << " in the chain, " 
        // 	    <<  cur_fac->involved_vars().size() << " variables"  << std::endl;
        //   std::cerr << "rank-range " << cur_fac->min_rank() << " - " << cur_fac->max_rank() 
        // 	    << ", address: " << cur_fac << std::endl;
        //   std::cerr << "variables: ";
        //   for (uint i=0; i < cur_fac->involved_vars().size(); i++)
        //     std::cerr << cur_fac->involved_vars()[i]->rank() << ", ";
        //   std::cerr << std::endl;
	
        //   std::cerr << "cast to one-of-n-factor: " << dynamic_cast<NaiveOneOfNTRWSFactor*>(cur_fac) << std::endl;
        //   if (cur_fac->prev_factor() != 0) {
        //     std::cerr << "cast of prev. factor to one-of-n-factor: " 
        // 	      << dynamic_cast<NaiveOneOfNTRWSFactor*>(cur_fac->prev_factor()) << std::endl;
        //     std::cerr << "cast of prev. factor to cardinality factor: " 
        // 	      << dynamic_cast<NaiveCardinalityTRWSFactor*>(cur_fac->prev_factor()) << std::endl;

        //     std::cerr << "prev. stored forward: " << cur_fac->prev_factor()->forward_ << std::endl;
        //   }

        //   std::cerr << "out-var used for dual bound: " << out_var->rank() << std::endl;
        //   std::cerr << "in-var used for dual bound: " << in_var->rank() << std::endl;

        //   if (cur_fac->prev_var() == 0)
        //     std::cerr << "no incoming var" << std::endl;
        //   else {
        //     std::cerr << "incoming var has rank " << cur_fac->prev_var()->rank() << std::endl;
        //     std::cerr << "previous factor has rank-range " << cur_fac->prev_factor()->min_rank()
        // 	      << " - " << cur_fac->prev_factor()->max_rank() << std::endl;
        //   }

        //   if (cur_fac->next_var() != 0) {
        //     std::cerr << "outgoing var has rank " << cur_fac->next_var()->rank() << std::endl;
        //     std::cerr << "following factor " << cur_fac->next_factor()
        // 	      << " has rank-range " << cur_fac->next_factor()->min_rank()
        // 	      << " - " << cur_fac->next_factor()->max_rank() 
        // 	      << " and " << cur_fac->next_factor()->involved_vars().size() << " variables" << std::endl;

        //     std::cerr << "next variables: ";
        //     for (uint i=0; i < cur_fac->next_factor()->involved_vars().size(); i++)
        //       std::cerr << cur_fac->next_factor()->involved_vars()[i]->rank() << ", ";
        //     std::cerr << std::endl;
        //   }
        //   else
        //     std::cerr << "no outgoing var" << std::endl;

        //   std::cerr << "previous forward vector: " << last_forward << std::endl;
        //   std::cerr << "message: " << new_forward << std::endl;

        //   assert(false);
        // }
        //END_DEBUG

        //restore parameters
        in_var->set_factor_params(cur_fac,temp);

        new_forward += out_var->factor_params(cur_fac);

        if (cur_fac->next_var() == 0) { //end of chain
          dual_bound += new_forward.min();
          break;
        }

        in_var = out_var;
        cur_fac = cur_fac->next_factor();
      }
    }
  }

  return dual_bound;
}

void NaiveFactorTRWS::create_chains(bool update_order, bool relaxed_monotonicity) {

  if (update_order) {
    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->set_rank(MAX_UINT);
    }
  }

  uint next_rank = 0;

  uint nChains = 0;
  uint nAtLeast5 = 0;
  uint nAtLeast10 = 0;
  uint nAtLeast25 = 0;

  for (uint f=0; f < nUsedFactors_; f++) {

    //std::cerr << "factor #" << f << "/" << nUsedFactors_ << std::endl;
    
    if (factor_[f]->prev_var() == 0 && factor_[f]->next_var() == 0) {

      uint length = 1;

      std::set<uint> current_ranks;

      //std::cerr << "extend lower" << std::endl;

      bool extension_found = true;

      NaiveTRWSFactor* cur_factor = factor_[f];

#if 1
      //extend lower end
      uint min_rank = factor_[f]->min_rank();

      while (extension_found) {

        extension_found = false;

        const Storage1D<NaiveTRWSVar*>& involved_vars = cur_factor->involved_vars();

        if (update_order) {
	  
          for (double k=0; k < involved_vars.size(); k++) {
            NaiveTRWSVar* var = involved_vars[k];
            if (var->rank() == MAX_UINT) {
              var->set_rank(next_rank);
              next_rank++;
            }
          }

          cur_factor->compute_rank_range();
        }

        for (uint k=0; k < involved_vars.size(); k++) {
          current_ranks.insert(involved_vars[k]->rank());
        }
	

        min_rank = std::min(min_rank,cur_factor->min_rank());

        for (double k=0; k < involved_vars.size(); k++) {
	
          NaiveTRWSVar* var = involved_vars[k];
	  
          uint var_rank = var->rank();

          uint second_rank = (relaxed_monotonicity) ? MAX_UINT : min_rank;
          if (relaxed_monotonicity) {

            for (std::set<uint>::iterator it = current_ranks.begin(); it != current_ranks.end(); it++) {
	      
              uint rank = *it;
              if (rank != var_rank && rank < second_rank)
                second_rank = rank;
            }

            //DEBUG
            //second_rank = min_rank;
            //END_DEBUG
	    
            // for (uint l=0; l < involved_vars.size(); l++) {
            //   if (l != k && involved_vars[l]->rank() < second_rank)
            // 	second_rank = involved_vars[l]->rank();
            // }
          }

          assert(var_rank < MAX_UINT);

          if (cur_factor->next_var() != 0 && var_rank >= cur_factor->next_var()->rank())
            continue;

          if (!relaxed_monotonicity && var_rank != cur_factor->min_rank()) {
            //requires passing more messages => needs to be specifically selected
            continue;
          }

          // //currently the incoming var needs to be the lowest ranked in a factor
          // if (relaxed_monotonicity && second_rank < var_rank)
          //   continue;

          const Storage1D<NaiveTRWSFactor*>& adjacent_factor = var->adjacent_factor();

          for (uint l=0; l < adjacent_factor.size(); l++) {
	    
            NaiveTRWSFactor* hyp_factor = adjacent_factor[l];
            const Storage1D<NaiveTRWSVar*>& hyp_involved_vars = hyp_factor->involved_vars();

            bool is_valid_extension = false;

            if (hyp_factor != cur_factor
                && hyp_factor->prev_var() == 0 && hyp_factor->next_var() == 0) {
	      
              is_valid_extension = true;
	      
              for (uint v=0; v < hyp_involved_vars.size(); v++) {
		
                uint rank = hyp_involved_vars[v]->rank();
		
                //if a variable has not been assigned a rank yet, it cannot be used in the downward direction
                if (rank == MAX_UINT || (rank != var_rank && rank >= second_rank)
                    || cur_factor->next_var() == hyp_involved_vars[v])
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
      uint max_rank = factor_[f]->max_rank();

      //std::cerr << "extend upper" << std::endl;

      cur_factor = factor_[f];
      extension_found = true;

      while (extension_found) {

        extension_found = false;

        const Storage1D<NaiveTRWSVar*>& involved_vars = cur_factor->involved_vars();

        if (update_order) {
	  
          for (double k=0; k < involved_vars.size(); k++) {
            NaiveTRWSVar* var = involved_vars[k];
            if (var->rank() == MAX_UINT) {
              var->set_rank(next_rank);
              next_rank++;
            }
          }

          cur_factor->compute_rank_range();
        }

        for (uint k=0; k < involved_vars.size(); k++) {
          current_ranks.insert(involved_vars[k]->rank());
        }

        max_rank = std::max(max_rank,cur_factor->max_rank());

        for (double k=0; k < involved_vars.size(); k++) {
	
          NaiveTRWSVar* var = involved_vars[k];
	  
          uint var_rank = var->rank();

          uint second_rank = (relaxed_monotonicity) ? 0 : max_rank;

          if (relaxed_monotonicity) {

            for (std::set<uint>::iterator it = current_ranks.begin(); it != current_ranks.end(); it++) {
	      
              uint rank = *it;
              if (rank != var_rank && rank > second_rank)
                second_rank = rank;
            }
	    
            //DEBUG
            //second_rank = max_rank;
            //END_DEBUG


            // if (cur_factor->prev_factor() != 0)
            //   second_rank = cur_factor->prev_factor()->max_rank();

            // for (uint l=0; l < involved_vars.size(); l++) {
            //   if (l != k && involved_vars[l]->rank() > second_rank)
            // 	second_rank = involved_vars[l]->rank();
            // }
          }

          assert(var_rank < MAX_UINT);

          if (cur_factor->prev_var() != 0 && var_rank <= cur_factor->prev_var()->rank())
            continue;

          if (!relaxed_monotonicity && var_rank != cur_factor->max_rank()) {
            //requires passing more messages => needs to be specifically selected
            continue;
          }

          const Storage1D<NaiveTRWSFactor*>& adjacent_factor = var->adjacent_factor();
	  
          for (uint l=0; l < adjacent_factor.size(); l++) {
	    
            NaiveTRWSFactor* hyp_factor = adjacent_factor[l];

            bool is_valid_extension = false;

            if (hyp_factor != cur_factor
                && hyp_factor->prev_var() == 0 && hyp_factor->next_var() == 0) {

              const Storage1D<NaiveTRWSVar*>& hyp_involved_vars = hyp_factor->involved_vars();
	      
              is_valid_extension = true;
	      
              for (uint v=0; v < hyp_involved_vars.size(); v++) {
		
                uint rank = hyp_involved_vars[v]->rank();
		
                if ((rank != var_rank && rank <= second_rank)
                    || cur_factor->prev_var() == hyp_involved_vars[v])
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


      if (length > 15)
        std::cerr << "chain length " << length << std::endl;
    }
  }

  std::cerr << nAtLeast5 << " chains have length at least 5." << std::endl;
  std::cerr << nAtLeast10 << " chains have length at least 10." << std::endl;
  std::cerr << nAtLeast25 << " chains have length at least 25." << std::endl;
  std::cerr << nChains << " chains in total, " << nUsedFactors_ << " factors" << std::endl;

  if (update_order) {
    assert(next_rank <= nUsedVars_); //variables without adjacent factors will not be assigned a rank

    if (next_rank < nUsedVars_)
      rank2var_.set_constant(MAX_UINT);

    for (uint v=0; v < nUsedVars_; v++) {
      if (var_[v]->rank() < MAX_UINT)
        rank2var_[var_[v]->rank()] = v;
    }
  }
}

void NaiveFactorTRWS::create_nonmonotonic_chains() {

  uint nChains = 0;
  uint nAtLeast5 = 0;
  uint nAtLeast10 = 0;
  uint nAtLeast25 = 0;

  for (uint f=0; f < nUsedFactors_; f++) {

    //std::cerr << "factor #" << f << "/" << nUsedFactors_ << std::endl;
    
    if (factor_[f]->prev_var() == 0 && factor_[f]->next_var() == 0) {

      uint length = 1;

      std::set<uint> current_ranks;

      //std::cerr << "extend lower" << std::endl;

      bool extension_found = true;

      NaiveTRWSFactor* cur_factor = factor_[f];

#if 1
      //extend lower end

      while (extension_found) {

        extension_found = false;

        const Storage1D<NaiveTRWSVar*>& involved_vars = cur_factor->involved_vars();

        for (uint k=0; k < involved_vars.size(); k++) {
          current_ranks.insert(involved_vars[k]->rank());
        }
	
        for (double k=0; k < involved_vars.size(); k++) {
	
          NaiveTRWSVar* var = involved_vars[k];

          if (var == cur_factor->next_var())
            continue;

          const Storage1D<NaiveTRWSFactor*>& adjacent_factor = var->adjacent_factor();

          for (uint l=0; l < adjacent_factor.size(); l++) {
	    
            NaiveTRWSFactor* hyp_factor = adjacent_factor[l];
            const Storage1D<NaiveTRWSVar*>& hyp_involved_vars = hyp_factor->involved_vars();

            bool is_valid_extension = false;

            if (hyp_factor != cur_factor
                && hyp_factor->prev_var() == 0 && hyp_factor->next_var() == 0) {
	      
              is_valid_extension = true;
	      
              for (uint v=0; v < hyp_involved_vars.size(); v++) {
		
                uint rank = hyp_involved_vars[v]->rank();
		
                if (hyp_involved_vars[v] != var && current_ranks.find(rank) != current_ranks.end())
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

        const Storage1D<NaiveTRWSVar*>& involved_vars = cur_factor->involved_vars();

        for (uint k=0; k < involved_vars.size(); k++) {
          current_ranks.insert(involved_vars[k]->rank());
        }

        for (double k=0; k < involved_vars.size(); k++) {
	
          NaiveTRWSVar* var = involved_vars[k];

          if (var == cur_factor->prev_var())
            continue;
	  
          const Storage1D<NaiveTRWSFactor*>& adjacent_factor = var->adjacent_factor();
	  
          for (uint l=0; l < adjacent_factor.size(); l++) {
	    
            NaiveTRWSFactor* hyp_factor = adjacent_factor[l];

            bool is_valid_extension = false;

            if (hyp_factor != cur_factor
                && hyp_factor->prev_var() == 0 && hyp_factor->next_var() == 0) {

              const Storage1D<NaiveTRWSVar*>& hyp_involved_vars = hyp_factor->involved_vars();
	      
              is_valid_extension = true;
	      
              for (uint v=0; v < hyp_involved_vars.size(); v++) {
		
                uint rank = hyp_involved_vars[v]->rank();
		
                if (hyp_involved_vars[v] != var && current_ranks.find(rank) != current_ranks.end())
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


      if (length > 15)
        std::cerr << "chain length " << length << std::endl;
    }
  }

  std::cerr << nAtLeast5 << " chains have length at least 5." << std::endl;
  std::cerr << nAtLeast10 << " chains have length at least 10." << std::endl;
  std::cerr << nAtLeast25 << " chains have length at least 25." << std::endl;
  std::cerr << nChains << " chains in total, " << nUsedFactors_ << " factors" << std::endl;
}

double NaiveFactorTRWS::optimize(uint nIter) {

  std::cerr.precision(10);

  bool relaxed_monotonicity = true;

  std::cerr << "optimize" << std::endl;

  create_chains(false,relaxed_monotonicity);

  std::cerr << "chains created" << std::endl;

  //WARNING: if optimize() was called before, the cost will be incorrectly scaled (fix TODO)
  for (uint v=0; v < nUsedVars_; v++)
    var_[v]->set_up_chains();

  double lower_bound = 1e-300;

  uint arg_min;

  for (uint iter=1; iter <= nIter; iter++) {

    //DEBUG
    // for (uint f=0; f < nUsedFactors_; f++) {
    //   factor_[f]->forward_ = 0.0;
    // }
    //END_DEBUG

    std::cerr << "****** iteration " << iter << std::endl;

    //forward
    double forward_lower = 0.0;

    for (uint i=0; i < nUsedVars_; i++) {

      // if ((i%1250) == 0)
      // 	std::cerr << "i: " << i << "/" << nUsedVars_ << std::endl;

      //std::cerr << "A" << std::endl;
      //std::cerr << "associated var: " << rank2var_[i] << std::endl;

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors

      NaiveTRWSVar* cur_var = var_[rank2var_[i]];

      //reparameterize (will automatically derive all relevant messages)
      forward_lower += cur_var->reparameterize_forward(relaxed_monotonicity);

      //std::cerr << "B" << std::endl;

      //average
      constant_energy_ += cur_var->average(arg_min);

      labeling_[rank2var_[i]] = arg_min;

      if (relaxed_monotonicity) {

        forward_lower += cur_var->special_fwd_cases();
      }
    }

    forward_lower += constant_energy_;
    std::cerr << "forward bound: " << forward_lower << std::endl;
    std::cerr << "constant: " << constant_energy_ << std::endl;

    lower_bound = forward_lower;

    //double independent_bound = dual_bound(true);
    //std::cerr << "independent bound: " << independent_bound << std::endl;

    //backward
    double backward_lower = 0.0;

    for (int i=nUsedVars_-1; i >= 0; i--) {

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors

      NaiveTRWSVar* cur_var = var_[rank2var_[i]];

      //reparameterize (will automatically derive all relevant messages)
      backward_lower += cur_var->reparameterize_backward(relaxed_monotonicity);

      //average
      constant_energy_ += cur_var->average(arg_min);

      labeling_[rank2var_[i]] = arg_min;
    }

    backward_lower += constant_energy_;
  }

  return lower_bound; 
}

double NaiveFactorTRWS::subgradient_opt(uint nIter, double start_step_size) {


  std::cerr.precision(10);

  bool relaxed_monotonicity = false;

  std::cerr << "subgradient optimization" << std::endl;

  std::cerr << nUsedFactors_ << " factors" << std::endl;

  //create_chains(false,relaxed_monotonicity);
  create_nonmonotonic_chains();

  std::cerr << "chains created" << std::endl;

  for (uint v=0; v < nUsedVars_; v++)
    var_[v]->set_up_chains();

  bool projective = true;

  double best_dual = -1e300;

  Math1D::Vector<uint> var_label(nUsedVars_);

  Storage1D<Math1D::Vector<uint> > factor_label(nUsedFactors_);
  
  for (uint f=0; f < nUsedFactors_; f++) {
    factor_label[f].resize(factor_[f]->involved_vars().size());
  }

  std::map<NaiveTRWSVar*,uint> var_num;
  for (uint v=0; v < nUsedVars_; v++)
    var_num[var_[v]] = v;

  std::map<NaiveTRWSFactor*,uint> factor_num;
  for (uint f=0; f < nUsedFactors_; f++)
    factor_num[factor_[f]] = f;

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "iteration #" << iter << std::endl;

    uint nDisagreements = 0;

    double step_size = start_step_size / iter;

    double cur_bound = 0.0;

    if (!projective) {
      for (uint v=0; v < nUsedVars_; v++) {
        uint cur_label = 0;
        cur_bound += var_[v]->dual_value(cur_label);
        var_label[v] = cur_label;
      }
    
      std::cerr << "A, intermediate bound: " << cur_bound << std::endl;
    }

    double factor_bound = 0.0;

    //uint nLongChainsProcessed = 0; 

    for (uint f=0; f < nUsedFactors_; f++) {

      if (factor_[f]->prev_var() == 0) {

        std::vector<NaiveTRWSFactor*> chain;
        std::vector<NaiveTRWSVar*> out_var;
	
        NaiveTRWSFactor* cur_factor = factor_[f];

        while (cur_factor != 0) {
          chain.push_back(cur_factor);
          out_var.push_back(cur_factor->next_var());
          cur_factor = cur_factor->next_factor();
        }

        uint chain_length = chain.size();

        //find chain end
        for (uint k=0; k < chain.back()->involved_vars().size(); k++) {
	  
          if (chain.back()->involved_vars()[k] != chain.back()->prev_var()) {
            out_var.back() = chain.back()->involved_vars()[k];
            break;
          }
        }

        //if (chain.size() > 1)
        //  nLongChainsProcessed++;
	
        Math1D::NamedVector<double> forward1(MAKENAME(forward1));
        Math1D::NamedVector<double> forward2(MAKENAME(forward2));

        NamedStorage1D<Math2D::Matrix<uint> > trace(chain_length,MAKENAME(trace));

        //std::cerr << "start forward" << std::endl;

        //compute forward
        chain[0]->compute_forward(out_var[0],forward1,trace[0]);
        forward1 += out_var[0]->factor_params(chain[0]);

        for (uint k=1; k < chain_length; k++) {

          Math1D::Vector<double>& last_forward = ((k % 2) == 1) ? forward1 : forward2;
          Math1D::Vector<double>& new_forward = ((k % 2) == 0) ? forward1 : forward2;

          //save 
          Math1D::Vector<double> temp = out_var[k-1]->factor_params(chain[k-1]);
          out_var[k-1]->set_factor_params(chain[k-1],last_forward);

          chain[k]->compute_forward(out_var[k],new_forward,trace[k]);
          new_forward += out_var[k]->factor_params(chain[k]);

          // if (nLongChainsProcessed == 1 && k <= 10) {
          //   std::cerr << "k: " << k << std::endl;
          //   std::cerr << "forward: " << new_forward << std::endl;
          //   std::cerr << "in var involved in " << out_var[k-1]->nChains() << " chains" << std::endl;
          //   exit(0);
          // }

          //reset
          out_var[k-1]->set_factor_params(chain[k-1],temp);
        }

        //std::cerr << "start traceback" << std::endl;

        //traceback
        Math1D::Vector<double>& total_forward = ((chain_length-1) % 2 == 0) ? forward1 : forward2;

        // if (chain.size() > 1 && nLongChainsProcessed <= 10) {

        //   std::cerr << "long chain #" << nLongChainsProcessed << ", length " << chain_length << std::endl;
        //   std::cerr << "total forward: " << total_forward << std::endl;
        // }

        if (chain_length == 1) {
          //std::cerr << "total forward: " << total_forward << std::endl;
          //std::cerr << "forward2: " << forward2 << std::endl;
          // for (uint v=0; v < chain[0]->involved_vars().size(); v++) {
          //   //std::cerr << "in param" << chain[0]->involved_vars()[v]->factor_params(chain[0]) << std::endl;
          //   std::cerr << "repar " << chain[0]->reparameterization( chain[0]->involved_vars()[v] ) << std::endl;
          // }
        }

        assert(total_forward.size() == out_var[chain_length-1]->nLabels());

        double best = 1e300;
        uint arg_best = MAX_UINT;
        for (uint l=0; l < total_forward.size(); l++) {
          if (total_forward[l] < best) {

            best = total_forward[l];
            arg_best = l;
          }
        }

        assert(arg_best < MAX_UINT);

        //std::cerr << "adding " << best << std::endl;
        factor_bound += best;

        cur_bound += best;

        for (int k=chain_length-1; k >= 0; k--) {

          //std::cerr << "k: " << k << std::endl;

          Math1D::Vector<uint>& labeling = factor_label[factor_num[chain[k]]];	  

          assert(labeling.size() == trace[k].xDim());

          for (uint v=0; v < trace[k].xDim(); v++)
            labeling[v] = trace[k](v,arg_best);

          //update arg_best
          if (k > 0) {
            for (uint v=0; v < trace[k].xDim(); v++) {
	      
              if (chain[k]->involved_vars()[v] == out_var[k-1]) {
                arg_best = labeling[v];
              }
            }
          }
        }

        if (chain_length == 1) {
          double check_cost = chain[0]->cost(factor_label[factor_num[chain[0]]]);
          //std::cerr << "check cost: " << check_cost << std::endl;
          assert(fabs( check_cost - best ) < 0.01);
        }

        //std::cerr << "traceback finished" << std::endl;
      }
    }

    std::cerr << "factor bound: " << factor_bound << std::endl;

    if (cur_bound > best_dual)
      best_dual = cur_bound;

    std::cerr << "cur bound: " << cur_bound << ", best ever: " << best_dual << std::endl;

    if (!projective) {
      //double total_repar = 0.0;

      for (uint f=0; f < nUsedFactors_; f++) {
	
        for (uint k=0; k < factor_label[f].size(); k++) {
	  
          NaiveTRWSVar* cur_var = factor_[f]->involved_vars()[k];
	  
          uint cur_fac_label = factor_label[f][k];
          uint cur_var_label = var_label[var_num[cur_var]];
	  
          //std::cerr << "factor label: " << cur_fac_label << std::endl;
          //std::cerr << "var label: " << cur_var_label << std::endl;
	  
          //total_repar += factor_[f]->reparameterization(cur_var)[cur_var_label];
          //total_repar -= factor_[f]->reparameterization(cur_var)[cur_fac_label];
	  
          if (cur_fac_label != cur_var_label) {
	    
            nDisagreements++;

            factor_[f]->reparameterization(cur_var)[cur_var_label] += step_size;
            factor_[f]->reparameterization(cur_var)[cur_fac_label] -= step_size;
          }
        }
      }
      
      std::cerr << nDisagreements << " disagreements" << std::endl;
      //std::cerr << "total repar cost: " << total_repar << std::endl;
    }
    else {
      for (uint f=0; f < nUsedFactors_; f++) {

        for (uint k=0; k < factor_label[f].size(); k++) {
	  
          NaiveTRWSVar* cur_var = factor_[f]->involved_vars()[k];
	  
          uint cur_fac_label = factor_label[f][k];
          factor_[f]->reparameterization(cur_var)[cur_fac_label] -= step_size;
        }
      }
      
      for (uint v=0; v < nUsedVars_; v++) {

        Math1D::Vector<double> sum(var_[v]->nLabels());
	
        for (uint k=0; k < var_[v]->adjacent_factor().size(); k++) {
          sum += var_[v]->adjacent_factor()[k]->reparameterization(var_[v]);
        }
	
        sum *= 1.0 / var_[v]->adjacent_factor().size();

        for (uint k=0; k < var_[v]->adjacent_factor().size(); k++) {
	  
          var_[v]->adjacent_factor()[k]->reparameterization(var_[v]) -= sum;
        }
      }
    }
  }

  if (projective) {

    for (uint f=0; f < nUsedFactors_; f++) {
      
      for (uint k=0; k < factor_label[f].size(); k++) {
	
        NaiveTRWSVar* cur_var = factor_[f]->involved_vars()[k];

        labeling_[var_num[cur_var]] = factor_label[f][k];
      }
    }
  }
  else
    labeling_ = var_label;

  return best_dual;

}

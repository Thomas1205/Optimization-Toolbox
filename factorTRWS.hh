/**** written by Thomas Schoenemann. Started as a private person without employment, August 2011 *****/
/**** continued at the University of Pisa, Italy, October - December 2011 ****/
/**** and at the Univeristy of Düsseldorf, Germany, 2012 ****/

/*** this class implements TRWS with factors of arbitrarily high order and singleton separators ***/

#ifndef FACTORTRWS_HH
#define FACTORTRWS_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"
#include "vardim_storage.hh"

class CumTRWSFactor;

/*** variable class ***/
class CumTRWSVar {
public:

  CumTRWSVar(const Math1D::Vector<float>& cost, uint rank);

  void add_cost(const Math1D::Vector<float>& add_cost);

  void set_up_chains();

  void add_factor(CumTRWSFactor* factor);

  double average(uint& arg_min);

  uint rank() const;

  void set_rank(uint rank);

  size_t nLabels() const;

  const Storage1D<CumTRWSFactor*>& adjacent_factor() const;

  const Math1D::Vector<double>& cost() const;

protected:

  Math1D::Vector<float> cost_;

  Math1D::Vector<double> cum_cost_;

  Storage1D<CumTRWSFactor*> adjacent_factor_;

  uint rank_;
  uint nChains_;
};

/******************************/

/*** abstract base class for factors ***/
/* abstract */ class CumTRWSFactor {
public:

  CumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars);

  virtual double compute_reparameterization(CumTRWSVar* var) = 0;
  
  uint min_rank() const;

  uint max_rank() const;

  void compute_rank_range();

  const Storage1D<CumTRWSVar*>& involved_vars() const;

  const Math1D::Vector<double>& reparameterization(const CumTRWSVar* var) const;

protected:

  Storage1D<CumTRWSVar*> involved_var_;

  Storage1D<Math1D::Vector<double> > reparameterization_;

  uint min_rank_;
  uint max_rank_;  
};

/******************************/

class GenericCumTRWSFactor : public CumTRWSFactor {
public:

  GenericCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const VarDimStorage<float>& cost);

  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const VarDimStorage<float> cost_;
};

/******************************/

/*** abstract base class for ternary factors ***/
/*abstract*/ class BinaryCumTRWSFactorBase : public CumTRWSFactor {
public:

  BinaryCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars);

  double compute_reparameterization(CumTRWSVar* var, const Math2D::Matrix<float>& cost);
};

/*** binary factor with costs stored explicitly ***/
class BinaryCumTRWSFactor : public BinaryCumTRWSFactorBase {
public:

  BinaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                      const Math2D::Matrix<float>& cost);
  
  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Math2D::Matrix<float> cost_;
};

/*** binary factor where a reference to the cost is stored (saves memory if you have many similar factors) ***/
class BinaryCumTRWSRefFactor : public BinaryCumTRWSFactorBase {
public:

  BinaryCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                         const Math2D::Matrix<float>& cost);
  
  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Math2D::Matrix<float>& cost_;
};


/****************/

/*** abstract base class for ternary factors ***/
/*abstract*/ class TernaryCumTRWSFactorBase : public CumTRWSFactor {
public:

  TernaryCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars);
  
  double compute_reparameterization(CumTRWSVar* var, const Math3D::Tensor<float>& cost);
};

/****************/

/*** ternary factor where costs are stored explicitly ***/
class TernaryCumTRWSFactor : public TernaryCumTRWSFactorBase {
public:

  TernaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math3D::Tensor<float>& cost);
  
  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Math3D::Tensor<float> cost_;
};

/****************/

/*** ternary factor where a reference to the cost is stored (saves memory if you have many similar factors) ***/
class TernaryCumTRWSRefFactor : public TernaryCumTRWSFactorBase {
public:

  TernaryCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math3D::Tensor<float>& cost);
  
  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Math3D::Tensor<float>& cost_;
};

/****************/

/*** ternary factor where the cost are based on second differences ***/
class SecondDiffCumTRWSFactor : public CumTRWSFactor {
public:

  SecondDiffCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, float lambda);

  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const float lambda_;
};

/****************/

/*** abstract base class for a 4th order factor ***/
/*abstract*/ class FourthOrderCumTRWSFactorBase : public CumTRWSFactor {
public:

  FourthOrderCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars);

  double compute_reparameterization(CumTRWSVar* var, const Storage1D<Math3D::Tensor<float> >& cost);
};

/*** 4th order factor with costs stored explicitly ***/
class FourthOrderCumTRWSFactor : public FourthOrderCumTRWSFactorBase {
public:

  FourthOrderCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                           const Storage1D<Math3D::Tensor<float> >& cost);

  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Storage1D<Math3D::Tensor<float> > cost_;
};


/*** 4th order factor where a reference to the cost is stored (saves memory if you have many similar factors) ***/
class FourthOrderCumTRWSRefFactor : public FourthOrderCumTRWSFactorBase {
public:

  FourthOrderCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                              const Storage1D<Math3D::Tensor<float> >& cost);

  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Storage1D<Math3D::Tensor<float> >& cost_;
};


/****************/

/*** 1-of-N constraint, all variables must be binary ***/
class OneOfNCumTRWSFactor : public CumTRWSFactor {
public:

  OneOfNCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars);

  virtual double compute_reparameterization(CumTRWSVar* var);
};

/****************/

/*** cardinality factor, all variables must be binary ***/
class CardinalityCumTRWSFactor : public CumTRWSFactor {
public:

  CardinalityCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math1D::Vector<float>& cost);
  
  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  const Math1D::Vector<float> cost_;
};

/****************/

/*** integer linear constraint factor for binary variables **/
class BILPCumTRWSFactor : public CumTRWSFactor {
public:

  BILPCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                    int rhs_lower = 0, int rhs_upper = 0);
  
  virtual double compute_reparameterization(CumTRWSVar* var);

protected:

  const Storage1D<bool> positive_;
  short rhs_lower_;
  short rhs_upper_;

  short range_;
  short zero_offset_;
};

/******************************/

/*** this is the main class ***/
class CumFactorTRWS {
public:

  /*** you presently need to provide upper bounds on the number of variables and factors you will add ***/
  CumFactorTRWS(uint nVars, uint nFactors);

  ~CumFactorTRWS();

  uint add_var(const Math1D::Vector<float>& cost);
  
  void add_generic_factor(const Math1D::Vector<uint>& vars, const VarDimStorage<float>& cost);

  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref=false);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  void add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                               const Storage1D<Math3D::Tensor<float> >& cost, bool ref=false);

  void add_one_of_n_factor(Math1D::Vector<uint>& var);

  void add_cardinality_factor(Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost);

  void add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                             int rhs_lower = 0, int rhs_upper = 0);

  double optimize(uint nIter);

  const Math1D::Vector<uint>& labeling();

protected:
  
  void add_factor(CumTRWSFactor* fac);

  Storage1D<CumTRWSVar*> var_;
  Storage1D<CumTRWSFactor*> factor_;

  Math1D::Vector<uint> labeling_;

  Math1D::Vector<uint> rank2var_;

  uint nUsedVars_;
  uint nUsedFactors_;

  bool optimize_called_;
};

#endif

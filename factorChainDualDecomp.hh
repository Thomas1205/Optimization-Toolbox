/******* written by Thomas Schoenemann as an employee of the University of Pisa, Italy, 2011 *****/

/**** this class implements the subgradient method for dual decomposition with factor chains *****/

#ifndef FACTOR_CHAIN_SUBGRADIENT_HH
#define FACTOR_CHAIN_SUBGRADIENT_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"

class ChainDDFactor;

//variable class
class ChainDDVar {
public:

  ChainDDVar(const Math1D::Vector<float>& cost);

  void add_cost(const Math1D::Vector<float>& cost);
  
  void add_factor(ChainDDFactor* factor);
  
  const Math1D::Vector<float>& cost() const;
  
  uint nLabels() const;

  const Storage1D<ChainDDFactor*>& neighboring_factor() const;

  uint nChains() const;

  void set_up_chains();

  double dual_value(uint& arg_min);

protected:

  Math1D::Vector<float> cost_;
  
  Storage1D<ChainDDFactor*> neighboring_factor_; 
};

//abstract base class for a factor
/* abstract */ class ChainDDFactor {
public:

  ChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace) = 0;

  virtual double cost(const Math1D::Vector<uint>& labeling) = 0;

  Math1D::Vector<double>& get_duals(ChainDDVar* var);

  ChainDDVar* prev_var() const;

  ChainDDVar* next_var() const;

  ChainDDFactor* prev_factor() const;

  ChainDDFactor* next_factor() const;

  void set_prev_var(ChainDDVar* var);

  void set_next_var(ChainDDVar* var);

  void set_prev_factor(ChainDDFactor* factor);

  void set_next_factor(ChainDDFactor* factor);

  const Storage1D<ChainDDVar*>& involved_vars();

protected:

  ChainDDVar* prev_var_;
  ChainDDVar* next_var_;

  ChainDDFactor* prev_factor_;
  ChainDDFactor* next_factor_;

  Storage1D<ChainDDVar*> involved_var_;
  Storage1D< Math1D::Vector<double> > dual_var_;
};

/***************************************/

// abstract base class for a binary factor
/*abstract */ class BinaryChainDDFactorBase: public ChainDDFactor {
public:
  
  BinaryChainDDFactorBase(const Storage1D<ChainDDVar*>& involved_vars);

  double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                         const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                         Math2D::Matrix<uint>& trace, const Math2D::Matrix<float>& cost);
};

// binary factor with costs stored explicitly
class BinaryChainDDFactor: public BinaryChainDDFactorBase {
public:
  
  BinaryChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math2D::Matrix<float>& cost);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Math2D::Matrix<float> cost_;
};

// binary factor where a reference to the cost is stored (saves memory if you have many similar factors)
class BinaryChainDDRefFactor: public BinaryChainDDFactorBase {
public:
  
  BinaryChainDDRefFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math2D::Matrix<float>& cost);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Math2D::Matrix<float>& cost_;
};

/***********/

// abstract base class for a ternary factor
class TernaryChainDDFactorBase : public ChainDDFactor {
public:

  TernaryChainDDFactorBase(const Storage1D<ChainDDVar*>& involved_vars);

  double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                         const Math1D::Vector<double>& prev_forward, const Math3D::Tensor<float>& cost,
                         Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);
};

// ternary factor with costs stored explicitly
class TernaryChainDDFactor: public TernaryChainDDFactorBase {
public:
  
  TernaryChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math3D::Tensor<float>& cost);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Math3D::Tensor<float> cost_;
};


// ternary factor storing only a reference to the cost (saves memory if you have many similar factors)
class TernaryChainDDRefFactor: public TernaryChainDDFactorBase {
public:
  
  TernaryChainDDRefFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math3D::Tensor<float>& tensor);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Math3D::Tensor<float>& cost_;
};

// ternary factor with cost based on second order differences
class SecondDiffChainDDFactor : public ChainDDFactor {
public:

  SecondDiffChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, float lambda);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);
  
  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const float lambda_;
};

/*****/

//abstract base class for a 4th order factor
/*abstract*/ class FourthOrderChainDDFactorBase : public ChainDDFactor {
public:

  FourthOrderChainDDFactorBase(const Storage1D<ChainDDVar*>& involved_vars);

  double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                         const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                         Math2D::Matrix<uint>& trace, const Storage1D<Math3D::Tensor<float> >& cost);
};

// 4th order factor with costs stored explicitly
class FourthOrderChainDDFactor : public FourthOrderChainDDFactorBase {
public:
  
  FourthOrderChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, 
                           const Storage1D<Math3D::Tensor<float> >& cost);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Storage1D<Math3D::Tensor<float> > cost_;
};

// 4th order factor storing only a reference to the cost (saves memory if you have many similar factors)
class FourthOrderChainDDRefFactor : public FourthOrderChainDDFactorBase {
public:
  
  FourthOrderChainDDRefFactor(const Storage1D<ChainDDVar*>& involved_vars, 
                           const Storage1D<Math3D::Tensor<float> >& cost);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Storage1D<Math3D::Tensor<float> >& cost_;
};

/******/

//1-of-N constraints (a special case of cardinality potentials), all variables must be binary
class OneOfNChainDDFactor : public ChainDDFactor {
public:

  OneOfNChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);
};

//cardinality factor, all variables must be binary 
class CardinalityChainDDFactor : public ChainDDFactor {
public:

  CardinalityChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Math1D::Vector<float>& cost);

  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:
  const Math1D::Vector<float> cost_;
};

//integer linear constraint factor for binary variables
class BILPChainDDFactor : public ChainDDFactor {
public:

  BILPChainDDFactor(const Storage1D<ChainDDVar*>& involved_vars, const Storage1D<bool>& positive,
                    int rhs_lower = 0, int rhs_upper = 0);
  
  virtual double compute_forward(const ChainDDVar* incoming, const ChainDDVar* outgoing,
                                 const Math1D::Vector<double>& prev_forward, Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling);

protected:

  const Storage1D<bool> positive_;
  short rhs_lower_;
  short rhs_upper_;

  short range_;
  short zero_offset_;
};



/***************************************/

// this is the main class
class FactorChainDualDecomposition {
public:

  /*** you can provide upper bounds on the number of variables and factors you want to add ***/
  FactorChainDualDecomposition(uint nVars, uint nFactors);

  ~FactorChainDualDecomposition();

  void add_var(const Math1D::Vector<float>& cost);

  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref=false);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  void add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                               const Storage1D<Math3D::Tensor<float> >& cost, bool ref=false);

  // all variables must be binary
  void add_one_of_n_factor(const Math1D::Vector<uint>& var);

  // all variables must be binary
  void add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost);

  // all variables must be binary
  void add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                             int rhs_lower = 0, int rhs_upper = 0);


  double optimize(uint nIter, double start_step_size = 1.0);

  const Math1D::Vector<uint>& labeling();

protected:

  void add_factor(ChainDDFactor* fac);

  void set_up_chains();

  Storage1D<ChainDDVar*> var_;
  Storage1D<ChainDDFactor*> factor_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedFactors_;

  bool optimize_called_;
};

#endif

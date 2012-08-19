/**** written by Thomas Schoenemann as a private person without employment, August 2011 *****/

/*** this class implements TRWS with factors of arbitrarily high order and singleton separators ***/

#ifndef FACTORTRWS_HH
#define FACTORTRWS_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"

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

//don't publish this!!
class TwoLevelCumTRWSFactor : public CumTRWSFactor {
public:

  TwoLevelCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                        const Math2D::Matrix<uint>& exception, double exception_cost, double std_cost);

  virtual double compute_reparameterization(CumTRWSVar* var);

protected:
  
  const Math2D::Matrix<uint> exception_;

  const double exception_cost_;
  const double std_cost_;
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

  void add_var(const Math1D::Vector<float>& cost);
  
  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref=false);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  void add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                               const Storage1D<Math3D::Tensor<float> >& cost, bool ref=false);

  void add_two_level_factor(Math1D::Vector<uint>& var, Math2D::Matrix<uint>& exceptions,
                            double exception_cost, double standard_cost);

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


/******************************/

//encodes what role a factor plays for a particular variable
//enum FactorRole {RoleChainStart,RoleChainEnd, RoleChainDeadEnd, RoleChainInLink, RoleChainOutLink};

class TRWSFactor;

class TRWSVar {
public:

  TRWSVar(const Math1D::Vector<double>& cost, uint rank);

  void add_cost(const Math1D::Vector<double>& add_cost);

  void set_up_chains();

  void add_factor(TRWSFactor* factor);

  double average_forward(double& bound, uint& arg_min);

  double average_backward(double& bound, uint& arg_min);

  uint rank() const;

  void set_rank(uint rank);

  size_t nLabels() const;

  const Storage1D<TRWSFactor*>& adjacent_factor() const;

  const Math1D::Vector<double>& cost() const;

protected:
  Math1D::Vector<double> cost_;

  Storage1D<TRWSFactor*> adjacent_factor_;

  uint rank_;
  uint nChains_;
};

class TRWSFactor {
public:

  TRWSFactor(const Storage1D<TRWSVar*>& involved_vars);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message) = 0;
  
  uint min_rank() const;

  uint max_rank() const;

  void compute_rank_range();

  // TRWSVar* prev_var() const;

  // TRWSVar* next_var() const;

  // TRWSFactor* prev_factor() const;

  // TRWSFactor* next_factor() const;

  // void set_prev_var(TRWSVar* var);

  // void set_next_var(TRWSVar* var);

  // void set_prev_factor(TRWSFactor* factor);

  // void set_next_factor(TRWSFactor* factor);

  const Storage1D<TRWSVar*>& involved_vars() const;

  const Math1D::Vector<double>& reparameterization(const TRWSVar* var) const;

protected:

  Storage1D<TRWSVar*> involved_var_;

  Storage1D<Math1D::Vector<double> > reparameterization_;

  // TRWSVar* prev_var_;
  // TRWSFactor* prev_factor_;

  // TRWSVar* next_var_;
  // TRWSFactor* next_factor_;

  uint min_rank_;
  uint max_rank_;  
};

/****************/

class BinaryTRWSFactor : public TRWSFactor {
public:

  BinaryTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                   const Math2D::Matrix<float>& cost);
  
  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:
  Math2D::Matrix<float> cost_;
};

/******************/

class TernaryTRWSFactorBase : public TRWSFactor {
public:
  
  TernaryTRWSFactorBase(const Storage1D<TRWSVar*>& involved_vars);
  
  double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message,
                                            const Math3D::Tensor<float>& cost);
};

class TernaryTRWSFactor : public TernaryTRWSFactorBase {
public:

  TernaryTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                    const Math3D::Tensor<float>& cost);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:
  Math3D::Tensor<float> cost_;
};

class TernaryTRWSRefFactor : public TernaryTRWSFactorBase {
public:
  
  TernaryTRWSRefFactor(const Storage1D<TRWSVar*>& involved_vars,
                       const Math3D::Tensor<float>& cost);
  
  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);
  
protected:
  const Math3D::Tensor<float>& cost_;
};

/***********************/

class SecondDiffTRWSFactor : public TRWSFactor {
public:

  SecondDiffTRWSFactor(const Storage1D<TRWSVar*>& involved_vars, float lambda);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:
  float lambda_;
};

/***********************/

class FourthOrderTRWSFactor : public TRWSFactor {
public:

  FourthOrderTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                        const Storage1D<Math3D::Tensor<float> >& cost);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:
  Storage1D<Math3D::Tensor<float> > cost_;
};

/***********************/

class TwoLevelTRWSFactor : public TRWSFactor {
public:

  TwoLevelTRWSFactor(const Storage1D<TRWSVar*>& involved_vars,
                     const Math2D::Matrix<uint>& exception, double exception_cost, double std_cost);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:
  
  Math2D::Matrix<uint> exception_;

  double exception_cost_;
  double std_cost_;
};

/***********************/

class OneOfNTRWSFactor : public TRWSFactor {
public:

  OneOfNTRWSFactor(const Storage1D<TRWSVar*>& involved_vars);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);
};

/***********************/

class CardinalityTRWSFactor : public TRWSFactor {
public:

  CardinalityTRWSFactor(const Storage1D<TRWSVar*>& involved_vars, const Math1D::Vector<float>& cost);
  
  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:
  Math1D::Vector<float> cost_;
};

/***********************/

class BILPTRWSFactor : public TRWSFactor {
public:

  BILPTRWSFactor(const Storage1D<TRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                 int rhs_lower = 0, int rhs_upper = 0);

  virtual double compute_message_and_reparameterize(TRWSVar* var, Math1D::Vector<double>& message);

protected:

  Storage1D<bool> positive_;
  short rhs_lower_;
  short rhs_upper_;

  short range_;
  short zero_offset_;
};


/***********************/

class FactorTRWS {
public:

  FactorTRWS(uint nVars, uint nFactors);

  ~FactorTRWS();

  void add_var(const Math1D::Vector<double>& cost);

  void add_var(const Math1D::Vector<float>& cost);
  
  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  void add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                               const Storage1D<Math3D::Tensor<float> >& cost);

  void add_two_level_factor(Math1D::Vector<uint>& var, Math2D::Matrix<uint>& exceptions,
                            double exception_cost, double standard_cost);

  void add_one_of_n_factor(Math1D::Vector<uint>& var);

  void add_cardinality_factor(Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost);

  void add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                             int rhs_lower = 0, int rhs_upper = 0);

  //void create_chains(bool update_order = false);

  double optimize(uint nIter);

  const Math1D::Vector<uint>& labeling();

protected:
  Storage1D<TRWSVar*> var_;
  Storage1D<TRWSFactor*> factor_;

  Math1D::Vector<uint> rank2var_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedFactors_;

  double constant_energy_;
};


/************************************************************************************/


class NaiveTRWSFactor;

struct FactorChainElement; 

class NaiveTRWSVar {
public:

  NaiveTRWSVar(const Math1D::Vector<float>& cost, uint rank);

  void add_cost(const Math1D::Vector<float>& add_cost);
  
  void set_up_chains();

  void add_factor(NaiveTRWSFactor* factor);

  double average(uint& arg_min);

  double reparameterize_forward(bool relaxed_monotonicity);

  double special_fwd_cases();

  double reparameterize_backward(bool relaxed_monotonicity);

  const Math1D::Vector<double>& factor_params(const NaiveTRWSFactor* factor) const;

  double add_factor_params(NaiveTRWSFactor* factor, const Math1D::Vector<double>& message);

  double normalize_factor_params(NaiveTRWSFactor* factor);

  void set_factor_params(NaiveTRWSFactor* factor, const Math1D::Vector<double>& params);

  uint rank() const;

  void set_rank(uint rank);

  size_t nLabels() const;

  uint nChains() const;

  const Storage1D<NaiveTRWSFactor*>& adjacent_factor() const;

  const Math1D::Vector<float>& cost() const;

  //used for the subgradient method
  double dual_value(uint& argmin) const;

protected:
  Math1D::Vector<float> cost_;

  Storage1D<NaiveTRWSFactor*> adjacent_factor_;

  Storage1D<FactorChainElement> chain_;

  uint rank_;
};


/* abstract*/ class NaiveTRWSFactor {
public:
  
  NaiveTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize) = 0;

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace) = 0;

  virtual double cost(const Math1D::Vector<uint>& labeling) const = 0;

  uint min_rank() const;

  uint max_rank() const;

  void compute_rank_range();

  NaiveTRWSVar* prev_var() const;

  NaiveTRWSVar* next_var() const;

  NaiveTRWSFactor* prev_factor() const;

  NaiveTRWSFactor* next_factor() const;

  void set_prev_var(NaiveTRWSVar* var);

  void set_next_var(NaiveTRWSVar* var);

  void set_prev_factor(NaiveTRWSFactor* factor);

  void set_next_factor(NaiveTRWSFactor* factor);

  const Storage1D<NaiveTRWSVar*>& involved_vars() const;

  const Math1D::Vector<double>& reparameterization(const NaiveTRWSVar* var) const;

  Math1D::Vector<double>& reparameterization(const NaiveTRWSVar* var);

  //DEBUG
  // void set_forward(double val);

  // double forward_;
  //END_DEBUG

protected:
  Storage1D<NaiveTRWSVar*> involved_var_;

  Storage1D<Math1D::Vector<double> > reparameterization_;

  NaiveTRWSVar* prev_var_;
  NaiveTRWSFactor* prev_factor_;

  NaiveTRWSVar* next_var_;
  NaiveTRWSFactor* next_factor_;

  uint min_rank_;
  uint max_rank_;
};

class NaiveBinaryTRWSFactor : public NaiveTRWSFactor {
public:
  NaiveBinaryTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                        const Math2D::Matrix<float>& cost);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:
  Math2D::Matrix<float> cost_;
};

class NaiveTernaryTRWSFactorBase : public NaiveTRWSFactor {
public:
  
  NaiveTernaryTRWSFactorBase(const Storage1D<NaiveTRWSVar*>& involved_vars);
  
  double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                            const Math3D::Tensor<float>& cost, bool reparameterize);

  //used for the subgradient method
  void compute_forward_base(const Math3D::Tensor<float>& cost, 
                            const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, 
                            Math2D::Matrix<uint>& trace);
};

class NaiveTernaryTRWSFactor : public NaiveTernaryTRWSFactorBase {
public:

  NaiveTernaryTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                         const Math3D::Tensor<float>& cost);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:
  Math3D::Tensor<float> cost_;
};

class NaiveTernaryTRWSRefFactor : public NaiveTernaryTRWSFactorBase {
public:
  
  NaiveTernaryTRWSRefFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                            const Math3D::Tensor<float>& cost);
  
  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);
  
  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);
  
  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:
  const Math3D::Tensor<float>& cost_;
};


class NaiveSecondDiffTRWSFactor : public NaiveTRWSFactor {
public:

  NaiveSecondDiffTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars, float lambda);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:
  float lambda_;
};


class NaiveFourthOrderTRWSFactor : public NaiveTRWSFactor {
public:
  NaiveFourthOrderTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                             const Storage1D<Math3D::Tensor<float> >& cost);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:
  Storage1D<Math3D::Tensor<float> > cost_;
};


class NaiveTwoLevelTRWSFactor : public NaiveTRWSFactor {
public:

  NaiveTwoLevelTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars,
                          const Math2D::Matrix<uint>& exception, double exception_cost, double std_cost);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:
  
  Math2D::Matrix<uint> exception_;

  double exception_cost_;
  double std_cost_;
};

class NaiveOneOfNTRWSFactor : public NaiveTRWSFactor {
public:

  NaiveOneOfNTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;
};

class NaiveCardinalityTRWSFactor : public NaiveTRWSFactor {
public:

  NaiveCardinalityTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars, const Math1D::Vector<float>& cost);
  
  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:

  Math1D::Vector<float> cost_;
};

class NaiveBILPTRWSFactor : public NaiveTRWSFactor {
public:

  NaiveBILPTRWSFactor(const Storage1D<NaiveTRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                      int rhs_lower = 0, int rhs_upper = 0);

  virtual double compute_message_and_reparameterize(NaiveTRWSVar* var, Math1D::Vector<double>& message,
                                                    bool reparameterize);

  //used for the subgradient method
  virtual void compute_forward(const NaiveTRWSVar* out_var, Math1D::Vector<double>& forward, Math2D::Matrix<uint>& trace);

  virtual double cost(const Math1D::Vector<uint>& labeling) const;

protected:

  Storage1D<bool> positive_;
  short rhs_lower_;
  short rhs_upper_;

  short range_;
  short zero_offset_;
};


struct FactorChainElement {
public:
  NaiveTRWSFactor* incoming_;
  NaiveTRWSFactor* outgoing_;

  Math1D::Vector<double> var_parameters_;
};

class NaiveFactorTRWS {
public:

  NaiveFactorTRWS(uint nVars, uint nFactors);

  ~NaiveFactorTRWS();

  void add_var(const Math1D::Vector<float>& cost);
  
  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  void add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                               const Storage1D<Math3D::Tensor<float> >& cost);

  void add_two_level_factor(Math1D::Vector<uint>& var, Math2D::Matrix<uint>& exceptions,
                            double exception_cost, double standard_cost);

  void add_one_of_n_factor(Math1D::Vector<uint>& var);

  void add_cardinality_factor(Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost);

  void add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                             int rhs_lower = 0, int rhs_upper = 0);

  void create_chains(bool update_order = false, bool relaxed_monotonicity = false);

  void create_nonmonotonic_chains();

  //returns the best found lower bound
  double optimize(uint nIter);

  //returns the best found lower bound
  double subgradient_opt(uint nIter, double start_step_size = 0.1);

  double dual_bound(bool check);

  const Math1D::Vector<uint>& labeling();

protected:
  Storage1D<NaiveTRWSVar*> var_;
  Storage1D<NaiveTRWSFactor*> factor_;

  Math1D::Vector<uint> rank2var_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedFactors_;

  double constant_energy_;
};

#endif

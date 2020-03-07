/******* written by Thomas Schoenemann. Started as a private person without employment, July 2011, *****/
/******* continued at the University of Pisa, Italy, October - December 2011 ***/
/******* continued at the University of Düsseldorf, Germany, October - December 2011 ***/

/***** implements solving linear programming relaxations with singleton separators ****/
/***** via MPLP, MSD and subgradient optimization for single factors in negative log-space ****/
/***** NOTE: subgradient optimization for single factors is generally very slow, so not recommended. ****/

#ifndef FACTOR_DUAL_OPT_HH
#define FACTOR_DUAL_OPT_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"
#include "vardim_storage.hh"

class DualFactorNode;

enum DualBCAMode {DUAL_BCA_MODE_MPLP, DUAL_BCA_MODE_MSD};

//variable class
class DualVariableNode {
public:

  DualVariableNode(const Math1D::Vector<float>& cost);

  void add_cost(const Math1D::Vector<float>& add_cost) noexcept;

  void add_factor(DualFactorNode* node) noexcept;

  double* get_dual_vars(const DualFactorNode* node) noexcept;

  const double* get_dual_vars(const DualFactorNode* node) const noexcept;

  double* get_dual_var_start() noexcept;

  const double* get_dual_var_start() const noexcept;

  void compute_message(const DualFactorNode* node, Math1D::Vector<double>& msg) const noexcept;

  double cost(uint label) const noexcept;

  const Math1D::Vector<float>& cost() const noexcept;

  double dual_value(uint& arg_min) const noexcept;

  uint nLabels() const noexcept;

  void init_dual_vars() noexcept;

  const Storage1D<DualFactorNode*>& neighboring_factor() const noexcept;

protected:

  Storage1D<DualFactorNode*> neighboring_factor_;

  Math2D::Matrix<double> dual_var_;

  Math1D::Vector<float> cost_;
};

/***************************/

// abstract base class for factors
/*abstract*/ class DualFactorNode {
public:

  DualFactorNode(const Storage1D<DualVariableNode*>& participating_vars);

  virtual ~DualFactorNode();

  virtual void update_duals(DualBCAMode mode) noexcept = 0;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept = 0;

  virtual double dual_value() const noexcept;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept = 0;

  //factors can modify the current labeling, e.g. to satisfy some constraints
  virtual bool process_labeling(Math1D::Vector<uint>& labels) const noexcept;

  const Storage1D<DualVariableNode*>& participating_nodes() const noexcept;

protected:

  Storage1D<DualVariableNode*> participating_var_;
};

/**************************/

//class for a factor of arbitrary order (run time will scale exponentially with the factor size)
class GenericDualFactorNode : public DualFactorNode {
public:

  GenericDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const VarDimStorage<float>& cost);

  //CAUTION: currently only MPLP supported, other modes will be ignored
  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:

  const VarDimStorage<float> cost_;
};

/**************************/

// abstract base class for a binary factor
/* abstract */ class BinaryDualFactorNodeBase : public DualFactorNode {
public:

  BinaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(const Math2D::Matrix<float>& cost, DualBCAMode mode) noexcept;

  double dual_value(const Math2D::Matrix<float>& cost) const noexcept;

  double compute_minimizer(const Math2D::Matrix<float>& cost, Math1D::Vector<uint>& min_labels) const noexcept;
};

/**************************/

// binary factor with costs stored explicitly
class BinaryDualFactorNode : public BinaryDualFactorNodeBase {
public:

  BinaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math2D::Matrix<float>& cost);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Math2D::Matrix<float> cost_;
};

/**************************/

// binary factor storing only a reference to the cost (saves memory if you have many similar factors)
class BinaryDualRefFactorNode : public BinaryDualFactorNodeBase {
public:

  BinaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math2D::Matrix<float>& cost);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Math2D::Matrix<float>& cost_;
};

/**************************/

//Potts factor (for two variables, i.e. not generalized Potts)
class PottsDualFactorNode: public DualFactorNode {
public:

  PottsDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const float lambda_;
};

/**************************/

//abstract base class for a ternary factor
/*abstract*/ class TernaryDualFactorNodeBase : public DualFactorNode {
public:

  TernaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(const Math3D::Tensor<float>& cost, DualBCAMode mode) noexcept;

  double dual_value(const Math3D::Tensor<float>& cost) const noexcept;

  double compute_minimizer(const Math3D::Tensor<float>& cost, Math1D::Vector<uint>& min_labels) const noexcept;
};

/**************************/

// ternary factor with costs stored explicitly
class TernaryDualFactorNode : public TernaryDualFactorNodeBase {
public:

  TernaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math3D::Tensor<float>& cost);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Math3D::Tensor<float> cost_;
};

/**************************/

// ternary factor storing only a reference to the cost (saves memory if you have many similar factors)
class TernaryDualRefFactorNode : public TernaryDualFactorNodeBase {
public:

  TernaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math3D::Tensor<float>& cost);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Math3D::Tensor<float>& cost_;
};

// ternary factor with costs based on second differences
class SecondDiffDualFactorNode : public DualFactorNode {
public:

  SecondDiffDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:

  const float lambda_;
};

/**************************/

//abstract base class for a fourth-order factor
/*abstract*/ class FourthOrderDualFactorNodeBase : public DualFactorNode {
public:

  FourthOrderDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(const Storage1D<Math3D::Tensor<float> >& cost, DualBCAMode mode) noexcept;

  double dual_value(const Storage1D<Math3D::Tensor<float> >& cost) const noexcept;

  double compute_minimizer(const Storage1D<Math3D::Tensor<float> >& cost, Math1D::Vector<uint>& min_labels) const noexcept;
};

/**************************/

//fourth order factor with costs stored explicitly
class FourthOrderDualFactorNode : public FourthOrderDualFactorNodeBase {
public:

  FourthOrderDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Storage1D<Math3D::Tensor<float> >& cost);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Storage1D<Math3D::Tensor<float> > cost_;
};

/**************************/

// fourth-order factor storing only a reference to the cost (saves memory if you have many similar factors)
class FourthOrderDualRefFactorNode : public FourthOrderDualFactorNodeBase {
public:

  FourthOrderDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Storage1D<Math3D::Tensor<float> >& cost);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Storage1D<Math3D::Tensor<float> >& cost_;
};


/**************************/

class GeneralizedPottsDualFactorNode : public DualFactorNode {
public:

  GeneralizedPottsDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda);

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:

  float lambda_;
};

/**************************/

//1-of-N constraints (a special case of cardinality potentials), all variables must be binary
class OneOfNDualFactorNode : public DualFactorNode {
public:

  //all variables must be binary
  OneOfNDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars);

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;
};

/**************************/

class CardinalityDualFactorNodeBase : public DualFactorNode {
public:

  CardinalityDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(DualBCAMode mode, const Math1D::Vector<float>& card_cost) noexcept;

  double dual_value(const Math1D::Vector<float>& card_cost) const noexcept;

  double compute_minimizer(Math1D::Vector<uint>& min_labels, const Math1D::Vector<float>& card_cost) const noexcept;

};

//cardinality factor, all variables must be binary
class CardinalityDualFactorNode : public CardinalityDualFactorNodeBase {
public:

  CardinalityDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math1D::Vector<float>& card_cost);

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Math1D::Vector<float> card_cost_;
};

/**************************/

//as above, but storing only a reference to the cost (saves memory if you have many similar factors)
class CardinalityDualFactorRefNode : public CardinalityDualFactorNodeBase {
public:

  CardinalityDualFactorRefNode(const Storage1D<DualVariableNode*>& participating_vars, const Math1D::Vector<float>& card_cost);

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  const Math1D::Vector<float>& card_cost_;
};

/**************************/

class NonbinaryCardinalityDualFactorNodeBase : public DualFactorNode {
public:

  NonbinaryCardinalityDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(DualBCAMode mode, const Math1D::Vector<float>& card_cost, const Math1D::Vector<uint>& level) noexcept;

  double dual_value(const Math1D::Vector<float>& card_cost, const Math1D::Vector<uint>& level) const noexcept;

  double compute_minimizer(Math1D::Vector<uint>& min_labels, const Math1D::Vector<float>& card_cost,
                           const Math1D::Vector<uint>& level) const noexcept;
};

/*****/

class NonbinaryCardinalityDualFactorNode : public NonbinaryCardinalityDualFactorNodeBase {
public:

  NonbinaryCardinalityDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                                     const Math1D::Vector<float>& card_cost, const Math1D::Vector<uint>& level);

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:

  const Math1D::Vector<float> card_cost_;
  const Math1D::Vector<uint> level_;
};

/**************************/


//special case of a cardinality potential where cost are 0-infty
class AllPosBILPConstraintDualFactorNode: public DualFactorNode {
public:

  //all variables must be binary
  AllPosBILPConstraintDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, int rhs_lower = 0, int rhs_upper = 0);

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:
  int rhs_lower_;
  int rhs_upper_;
};

/**************************/

//integer linear constraint factor for binary variables
class BILPConstraintDualFactorNode : public DualFactorNode {
public:

  //all variables must be binary
  BILPConstraintDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                               const Storage1D<bool>& positive, int rhs_lower = 0, int rhs_upper = 0);

  virtual double cost(const Math1D::Vector<uint>& labels) const noexcept override;

  virtual void update_duals(DualBCAMode mode) noexcept override;

  virtual double dual_value() const noexcept override;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const noexcept override;

protected:

  uint nPos_;

  int rhs_lower_;
  int rhs_upper_;
};

/**************************/

//this is the main class
class FactorDualOpt {
public:

  FactorDualOpt(uint nVars = 0, uint nFactors=0);

  //delete all allocated owned var. and factor nodes
  ~FactorDualOpt();

  /***** Setting up a factor graph ******/

  //return var. number
  uint add_var(const Math1D::Vector<float>& unary_cost);

  DualVariableNode* get_var(uint var_num);

  //returns factor number
  uint add_generic_factor(const Math1D::Vector<uint> var, const VarDimStorage<float>& cost);

  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  //returns factor number
  uint add_generic_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref=false);

  //returns factor number
  uint add_potts_factor(uint var1, uint var2, double lambda);

  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  //returns factor number
  uint add_generic_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  //returns factor number
  uint add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  //returns factor number
  uint add_generic_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                       const Storage1D<Math3D::Tensor<float> >& cost, bool ref=false);

  //returns factor number
  uint add_generalized_potts_factor(const Math1D::Vector<uint>& var, double lambda);

  //all participating variables must be binary
  //returns factor number
  uint add_one_of_N_factor(const Math1D::Vector<uint>& var);

  //all participating variables must be binary
  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  //returns factor number
  uint add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost, bool ref=false);

  uint add_nonbinary_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost,
                                        const Math1D::Vector<uint>& level);

  //all participating variables must be binary
  //returns factor number
  uint add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                             int rhs_lower = 0, int rhs_upper = 0);

  //CAUTION: the passed factor is deleted together with all truly owned ones.
  //returns factor number
  uint pass_in_factor_node(DualFactorNode* factor);

  /**** run inference ***/

  //returns the computed lower bound
  double dual_bca(uint nIter, DualBCAMode mode = DUAL_BCA_MODE_MPLP, bool init=true, bool quiet = false, bool terminate_when_no_progress = false);

  //returns the best found lower bound
  double subgradient_opt(uint nIter, double start_step_size = 0.1);

  /**** get the state of the solver ****/

  double labeling_energy();

  const Math1D::Vector<uint>& labeling();

  void icm(uint iter);

  uint best_of_n(uint fac_num);

protected:

  //void process_labeling();

  uint add_factor(DualFactorNode* node);

  Storage1D<DualVariableNode*> var_;
  Storage1D<DualFactorNode*> factor_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedFactors_;
};


#endif

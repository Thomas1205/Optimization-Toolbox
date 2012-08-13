/******* written by Thomas Schoenemann as a private person without employment, July 2011 *****/

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

  void add_cost(const Math1D::Vector<float>& add_cost);
  
  void add_factor(DualFactorNode* node);
  
  double* get_dual_vars(const DualFactorNode* node);
  
  const double* get_dual_vars(const DualFactorNode* node) const;

  double* get_dual_var_start();

  const double* get_dual_var_start() const;

  void compute_message(const DualFactorNode* node, Math1D::Vector<double>& msg) const;

  double cost(uint label) const;

  const Math1D::Vector<float>& cost() const;

  double dual_value(uint& arg_min) const;
  
  uint nLabels() const;

  void init_dual_vars();

  const Storage1D<DualFactorNode*>& neighboring_factor() const;

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
  
  virtual void update_duals(DualBCAMode mode) = 0;
  
  virtual double cost(const Math1D::Vector<uint>& labels) const = 0;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const = 0;

  //factors can modify the current labeling, e.g. to satisfy some constraints
  virtual bool process_labeling(Math1D::Vector<uint>& labels) const;

  const Storage1D<DualVariableNode*>& participating_nodes() const;
  
protected:
  
  Storage1D<DualVariableNode*> participating_var_;
};

/**************************/

//class for a factor of arbitrary order (run time will scale exponentially with the factor size)
class GenericDualFactorNode : public DualFactorNode {
public:

  GenericDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const VarDimStorage<float>& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  
  VarDimStorage<float> cost_;
};

/**************************/

// abstract base class for a binary factor
/* abstract */ class BinaryDualFactorNodeBase : public DualFactorNode {
public:

  BinaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(const Math2D::Matrix<float>& cost, DualBCAMode mode);
  
  double dual_value(const Math2D::Matrix<float>& cost) const;

  double compute_minimizer(const Math2D::Matrix<float>& cost, Math1D::Vector<uint>& min_labels) const;
};

/**************************/

// binary factor with costs stored explicitly
class BinaryDualFactorNode : public BinaryDualFactorNodeBase {
public:

  BinaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math2D::Matrix<float>& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  Math2D::Matrix<float> cost_;
};

/**************************/

// binary factor storing only a reference to the cost (saves memory if you have many similar factors)
class BinaryDualRefFactorNode : public BinaryDualFactorNodeBase {
public:

  BinaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math2D::Matrix<float>& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  const Math2D::Matrix<float>& cost_;
};

/**************************/

//Potts factor (for two variables, i.e. not generalized Potts)
class PottsDualFactorNode: public DualFactorNode {
public:

  PottsDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda);

  virtual void update_duals(DualBCAMode mode);

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  float lambda_;
};

/**************************/

//abstract base class for a ternary factor
/*abstract*/ class TernaryDualFactorNodeBase : public DualFactorNode {
public: 

  TernaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(const Math3D::Tensor<float>& cost, DualBCAMode mode);
  
  double dual_value(const Math3D::Tensor<float>& cost) const;

  double compute_minimizer(const Math3D::Tensor<float>& cost, Math1D::Vector<uint>& min_labels) const;
};

/**************************/

// ternary factor with costs stored explicitly
class TernaryDualFactorNode : public TernaryDualFactorNodeBase {
public:

  TernaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math3D::Tensor<float>& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;
  
  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  Math3D::Tensor<float> cost_;
};

/**************************/

// ternary factor storing only a reference to the cost (saves memory if you have many similar factors)
class TernaryDualRefFactorNode : public TernaryDualFactorNodeBase {
public:

  TernaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, const Math3D::Tensor<float>& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;
  
  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  const Math3D::Tensor<float>& cost_;
};

// ternary factor with costs based on second differences
class SecondDiffDualFactorNode : public DualFactorNode {
public:

  SecondDiffDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;
  
  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;
  
protected:

  float lambda_;
};

/**************************/

//abstract base class for a fourth-order factor
/*abstract*/ class FourthOrderDualFactorNodeBase : public DualFactorNode {
public: 

  FourthOrderDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars);

  void update_duals(const Storage1D<Math3D::Tensor<float> >& cost, DualBCAMode mode);
  
  double dual_value(const Storage1D<Math3D::Tensor<float> >& cost) const;

  double compute_minimizer(const Storage1D<Math3D::Tensor<float> >& cost, 
                           Math1D::Vector<uint>& min_labels) const;

};

/**************************/

//fourth order factor with costs stored explicitly
class FourthOrderDualFactorNode : public FourthOrderDualFactorNodeBase {
public:

  FourthOrderDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                            const Storage1D<Math3D::Tensor<float> >& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;
  
protected:
  Storage1D<Math3D::Tensor<float> > cost_;
};

/**************************/

// fourth-order factor storing only a reference to the cost (saves memory if you have many similar factors)
class FourthOrderDualRefFactorNode : public FourthOrderDualFactorNodeBase {
public:

  FourthOrderDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                               const Storage1D<Math3D::Tensor<float> >& cost);

  virtual void update_duals(DualBCAMode mode);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;
  
protected:
  const Storage1D<Math3D::Tensor<float> >& cost_;
};

/**************************/

//1-of-N constraints (a special case of cardinality potentials), all variables must be binary
class OneOfNDualFactorNode : public DualFactorNode {
public:

  //all variables must be binary
  OneOfNDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars);

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual void update_duals(DualBCAMode mode);

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;
};

/**************************/

//cardinality factor, all variables must be binary 
class CardinalityDualFactorNode : public DualFactorNode {
public:

  CardinalityDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                            const Math1D::Vector<float>& card_cost);

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual void update_duals(DualBCAMode mode);

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  const Math1D::Vector<float> card_cost_;
};

/**************************/

//integer linear constraint factor for binary variables
class BILPConstraintDualFactorNode: public DualFactorNode {
public:

  //all variables must be binary
  BILPConstraintDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                               const Storage1D<bool>& positive, int rhs_lower = 0, int rhs_upper = 0);

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual void update_duals(DualBCAMode mode);

  virtual double dual_value() const;

  virtual double compute_minimizer(Math1D::Vector<uint>& min_labels) const;

protected:
  const Storage1D<bool> positive_;
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

  //NOTE: after calling this routine, owned vars can no longer be created
  uint pass_in_var_node(DualVariableNode* var);

  uint add_generic_factor(const Math1D::Vector<uint> var, VarDimStorage<float>& cost);

  uint add_generic_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref=false);
  
  uint add_potts_factor(uint var1, uint var2, double lambda);

  //return factor number
  uint add_generic_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  uint add_second_diff_factor(uint var1, uint var2, uint var3, float lambda);

  //return factor number
  uint add_generic_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                       const Storage1D<Math3D::Tensor<float> >& cost, bool ref=false);


  //all participating variables must be binary
  uint add_one_of_N_factor(const Math1D::Vector<uint>& var);

  //all participating variables must be binary
  uint add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost);

  //all participating variables must be binary
  uint add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                             int rhs_lower = 0, int rhs_upper = 0);

  //NOTE: after calling this routine, owned factors can no longer be created
  uint pass_in_factor_node(DualFactorNode* factor);

  /**** run inference ***/

  //returns the computed lower bound
  double dual_bca(uint nIter, DualBCAMode mode = DUAL_BCA_MODE_MPLP, bool init=true, bool quiet = false);

  //returns the best found lower bound
  double subgradient_opt(uint nIter, double start_step_size = 0.1);

  /**** get the state of the solver ****/

  double labeling_energy();

  const Math1D::Vector<uint>& labeling();

  void icm(uint iter);

protected:

  uint add_factor(DualFactorNode* node, bool owned = true);

  uint first_shared_var_;
  uint first_shared_factor_;

  Storage1D<DualVariableNode*> var_;
  Storage1D<DualFactorNode*> factor_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedFactors_;
};


#endif

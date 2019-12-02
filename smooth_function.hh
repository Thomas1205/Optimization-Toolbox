/**** written by Thomas Schoenemann as a private person without employment, October 2013 ****/


#ifndef SMOOTH_FUNCTION_HH
#define SMOOTH_FUNCTION_HH

#include "vector.hh"
#include "storage_util.hh"

/*** this class represents a function of <code> n_ </code> variables that is continuous and differentiable everywhere ***/
/*abstract*/ class SmoothFunction {
public:

  SmoothFunction(uint n);

  void set_max_alpha(double max_alpha);

  void set_backtracking_factor(double backtrack_factor);

  void set_cutoff_factor(double cutoff_factor);

  void set_cutoff_sqr_grad_norm(double cutoff_sqr_grad_norm);

  //returns the computed function value
  double minimize_nlcg(Math1D::Vector<double>& x, uint nIter);

  //returns the computed function value
  double minimize_l_bfgs(Math1D::Vector<double>& x, int L, uint nIter);

protected:

  virtual double function_value(const Math1D::Vector<double>& x) const = 0;

  virtual double hyp_function_value(const Math1D::Vector<double>& x, const Math1D::Vector<double>& addon, double alpha) const;

  virtual void compute_gradient(const Math1D::Vector<double>& x, Math1D::Vector<double>& grad) const = 0;


  const uint n_; //the number of variables of the function

  bool quiet_;

  double max_alpha_; //the upper range for the line search
  double backtrack_factor_; //the backtracking factor in the line search
  double cutoff_factor_; //how much progress we need to make to accept a step
  double grad_norm_cutoff_; //stop the minimization when the squared gradient norm is less than this value
};


#endif

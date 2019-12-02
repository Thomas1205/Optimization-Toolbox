/**** written by Thomas Schoenemann as a private person without employment, October 2013 ****/

#include "smooth_function.hh"

SmoothFunction::SmoothFunction(uint n) :
  n_(n), quiet_(false), max_alpha_(1.0), backtrack_factor_(0.75), cutoff_factor_(0.05), grad_norm_cutoff_(1-4)
{}

void SmoothFunction::set_max_alpha(double max_alpha)
{
  max_alpha_ = max_alpha;
}

void SmoothFunction::set_backtracking_factor(double backtrack_factor)
{
  backtrack_factor_ = backtrack_factor;
}

void SmoothFunction::set_cutoff_factor(double cutoff_factor)
{
  cutoff_factor_ = cutoff_factor;
}

void SmoothFunction::set_cutoff_sqr_grad_norm(double cutoff_sqr_grad_norm)
{
  grad_norm_cutoff_ = cutoff_sqr_grad_norm;
}


/*virtual*/
double SmoothFunction::hyp_function_value(const Math1D::Vector<double>& x, const Math1D::Vector<double>& addon, double alpha) const
{


  Math1D::Vector<double> hyp_x(n_);
  for (uint k=0; k < n_; k++)
    hyp_x[k] = x[k] + alpha*addon[k];

  return function_value(hyp_x);
}


//returns the computed function value
double SmoothFunction::minimize_nlcg(Math1D::Vector<double>& x, uint nIter)
{


  Math1D::Vector<double> search_direction(n_);

  Math1D::Vector<double> gradient[2];
  gradient[0].resize(n_);
  gradient[1].resize(n_);

  // init search direction with the current gradient
  compute_gradient(x,gradient[0]);

  double cur_energy = function_value(x);

  search_direction = gradient[0];
  negate(search_direction);

  double prev_grad_norm = gradient[0].sqr_norm();

  if (!quiet_)
    std::cerr.precision(10);

  uint cur_idx = 0;

  double max_alpha = max_alpha_;

  for (uint iter=1; iter <= nIter; iter++) {

    if (!quiet_)
      std::cerr << "### iteration " << iter << ", energy: " << cur_energy << std::endl;

    const Math1D::Vector<double>& prev_gradient = gradient[cur_idx];
    cur_idx = 1-cur_idx;
    Math1D::Vector<double>& cur_gradient = gradient[cur_idx];

    //a) line search along the current search direction
    double best_alpha = max_alpha;

    double best_energy = hyp_function_value(x,search_direction,best_alpha);

    double alpha = best_alpha;

    //std::cerr << "alpha: " << alpha << ", hyp energy: " << best_energy << std::endl;

    double cutoff_offset = 0.0;
    for (uint k=0; k < n_; k++)
      cutoff_offset += prev_gradient[k] * search_direction[k];
    cutoff_offset *= cutoff_factor_;

    assert(cutoff_offset < 0.0);

    bool wolfe_satisfied = false;

    while (alpha > 1e-12) {

      alpha *= backtrack_factor_;

      double hyp_energy = hyp_function_value(x,search_direction,alpha);

      // std::cerr << "alpha: " << alpha << ", hyp energy: " << hyp_energy << ", threshold: "
      // 		<< (cur_energy + alpha*cutoff_offset) << ", satisfied: "
      // 		<< (hyp_energy <= cur_energy + alpha*cutoff_offset) << std::endl;

      //Wolfe conditions, part 1:
      if (hyp_energy <= cur_energy + alpha*cutoff_offset) {
        wolfe_satisfied = true; //TODO: second part of strong Wolfe conditions

        if (hyp_energy < best_energy) {
          best_energy = hyp_energy;
          best_alpha = alpha;
        }
        else
          break;
      }
      else if (wolfe_satisfied)
        break;
    }

    //b) update <code> x </code> and the gradient at <code> x </code>
    if (best_energy <= cur_energy) {

      for (uint k=0; k < n_; k++)
        x[k] += best_alpha * search_direction[k];

      cur_energy = best_energy;

      if (best_alpha == max_alpha)
        max_alpha *= 1.5;
      else if (best_alpha < backtrack_factor_*backtrack_factor_*max_alpha)
        max_alpha *= backtrack_factor_;
    }
    else if (!quiet_)
      std::cerr << "WARNING: not updating" << std::endl;

    compute_gradient(x,cur_gradient);

    double grad_norm = cur_gradient.sqr_norm();

    if (!quiet_)
      std::cerr << "sqr grad norm: " << grad_norm << std::endl;

    if (grad_norm < grad_norm_cutoff_)
      break; //problem solved to sufficient accuracy

    double beta_fr = grad_norm / prev_grad_norm; //Fletcher-Reeves variant

    double beta_pr = 0.0;
    double beta_hs = 0.0;
    double beta_hs_denom = 0.0;
    for (uint k=0; k < n_; k++) {

      const double diff = cur_gradient[k] - prev_gradient[k];
      beta_pr += diff * cur_gradient[k]; //this is also the numerator of HS
      beta_hs_denom += diff * search_direction[k];
    }
    beta_hs = beta_pr / beta_hs_denom;

    beta_pr /= prev_grad_norm; //Polack-Ribiere variant

    beta_pr = std::max(0.0,beta_pr);
    //beta_pr = std::max(-beta_fr,beta_pr);
    beta_hs = std::max(0.0,beta_hs); //Hestenes-Stiefel variant
    //beta_hs = std::max(-beta_fr,beta_hs);

    //if (beta_pr > beta_fr)
    //  beta_pr = beta_fr; //PR+

    //if (beta_hs > beta_fr)
    //  beta_hs = beta_fr; //HS+

    //double beta = beta_fr;
    //double beta = beta_pr;
    double beta = beta_hs;

    //experimental safe guard

    if (alpha*cutoff_offset > -0.01) {
      //if (true) {

      if (!quiet_) {
        std::cerr << "RESET" << std::endl;
        std::cerr << "best alpha: " << best_alpha << std::endl;
      }
      beta = 0.0;
    }

    for (uint k=0; k < n_; k++)
      search_direction[k] = (-cur_gradient[k]) + beta * search_direction[k];

    prev_grad_norm = grad_norm;
  }

  return cur_energy;
}

//returns the computed function value
double SmoothFunction::minimize_l_bfgs(Math1D::Vector<double>& x, int L, uint nIter)
{

  assert(L >= 1);

  Storage1D<Math1D::Vector<double> > grad_diff(L);
  Storage1D<Math1D::Vector<double> > step(L);

  Math1D::Vector<double> rho(L);

  for (int l=0; l < L; l++) {

    grad_diff[l].resize(n_);
    step[l].resize(n_);
  }

  if (!quiet_)
    std::cerr.precision(10);

  Math1D::Vector<double> search_direction(n_);

  Math1D::Vector<double> cur_gradient(n_);

  double cur_energy = function_value(x);

  int start_iter = 1;

  for (int iter=1; iter <= int(nIter); iter++) {

    if (!quiet_)
      std::cerr << "### iteration " << iter << ", energy: " << cur_energy << std::endl;

    compute_gradient(x,cur_gradient);

    double sqr_grad_norm = cur_gradient.sqr_norm();
    if (sqr_grad_norm < grad_norm_cutoff_)
      break; //problem solved to sufficient accuracy

    double cur_curv = 0.0;

    if (iter > 1) {
      //update grad_diff and rho
      uint cur_l = (iter-1) % L;
      Math1D::Vector<double>& cur_grad_diff = grad_diff[cur_l];
      const Math1D::Vector<double>& cur_step = step[cur_l];

      double cur_rho = 0.0;

      for (uint k=0; k < n_; k++) {

        //cur_grad_diff was set to minus the previous gradient at the end of the previous iteration
        cur_grad_diff[k] += cur_gradient[k];
        cur_rho += cur_grad_diff[k] * cur_step[k];
      }

      cur_curv = cur_rho / cur_grad_diff.sqr_norm();

      if (cur_curv <= 0.0) {
        //this can happen if the function is not (strictly) convex, since we do not enforce Wolfe part 2
        //  (would need Algorithm 3.5 from [Nocedal & Wright] to enforce that, backtracking line search is NOT enough)
        // Our solution is to simply restart L-BFGS

        start_iter = iter;
      }

      rho[cur_l] = 1.0 / cur_rho;
    }

    //a) determine the search direction via L-BFGS
    search_direction = cur_gradient;

    if (iter > start_iter) {

      const int cur_first_iter = std::max<int>(start_iter,iter-L);

      Math1D::Vector<double> alpha(L);

      //first loop in Algorithm 7.4 from [Nocedal & Wright]
      for (int prev_iter = iter-1; prev_iter >= cur_first_iter; prev_iter--) {

        uint prev_l = prev_iter % L;

        const Math1D::Vector<double>& cur_step = step[prev_l];
        const Math1D::Vector<double>& cur_grad_diff = grad_diff[prev_l];

        double cur_alpha = 0.0;
        for (uint k=0; k < n_; k++) {
          cur_alpha += search_direction[k] * cur_step[k];
        }
        cur_alpha *= rho[prev_l];
        alpha[prev_l] = cur_alpha;

        for (uint k=0; k < n_; k++) {
          search_direction[k] -= cur_alpha * cur_grad_diff[k];
        }
      }

      //we use a scaled identity as base matrix (q=r=search_direction)
      search_direction *= cur_curv;

      //second loop in Algorithm 7.4 from [Nocedal & Wright]
      for (int prev_iter = cur_first_iter; prev_iter < iter; prev_iter++) {

        uint prev_l = prev_iter % L;

        const Math1D::Vector<double>& cur_step = step[prev_l];
        const Math1D::Vector<double>& cur_grad_diff = grad_diff[prev_l];

        double beta = 0.0;
        for (uint k=0; k < n_; k++) {
          beta += search_direction[k] * cur_grad_diff[k];
        }
        beta *= rho[prev_l];

        const double gamma = alpha[prev_l] - beta;

        for (uint k=0; k < n_; k++) {
          search_direction[k] += cur_step[k] * gamma;
        }
      }

    }
    negate(search_direction);

    //b) line search along the current search direction
    double best_alpha = 1.0;
    double best_energy = hyp_function_value(x,search_direction,best_alpha);

    double alpha = best_alpha;

    //std::cerr << "alpha: " << alpha << ", hyp energy: " << best_energy << std::endl;

    double cutoff_offset = 0.0;
    for (uint k=0; k < n_; k++)
      cutoff_offset += cur_gradient[k] * search_direction[k];
    cutoff_offset *= 0.05; //0.001;

    if (cutoff_offset >= 0.0) {
      //this can happen if the function is not strongly convex and we do not enforce part 2 of the Wolfe conditions.
      // in such a case we switch to a steepest descent iteration

      search_direction = cur_gradient;
      negate(search_direction);

      cutoff_offset = - cur_gradient.sqr_norm();
    }

    cutoff_offset *= 0.05; //0.001;


    bool wolfe_satisfied = false;

    while (alpha > 1e-12) {

      alpha *= backtrack_factor_;

      double hyp_energy = hyp_function_value(x,search_direction,alpha);

      // std::cerr << "alpha: " << alpha << ", hyp energy: " << hyp_energy << ", threshold: "
      // 		<< (cur_energy + alpha*cutoff_offset) << ", satisfied: "
      // 		<< (hyp_energy <= cur_energy + alpha*cutoff_offset) << std::endl;

      //Wolfe conditions, part 1:
      if (hyp_energy <= cur_energy + alpha*cutoff_offset) {
        wolfe_satisfied = true;

        //NOTE: part 2 of the Wolfe conditions is difficult to enforce, see Algorithm 3.5 in [Nocedal & Wright]

        if (hyp_energy < best_energy) {
          best_energy = hyp_energy;
          best_alpha = alpha;
        }
        else
          break;
      }
      else if (wolfe_satisfied)
        break;
    }

    //c) update the variables and the step vectors
    if (best_energy <= cur_energy) {

      uint cur_l = (iter % L);

      Math1D::Vector<double>& cur_step = step[cur_l];
      Math1D::Vector<double>& cur_grad_diff = grad_diff[cur_l];

      for (uint k=0; k < n_; k++) {
        double step = best_alpha * search_direction[k];
        cur_step[k] = step;
        x[k] += step;

        //prepare for the next iteration
        cur_grad_diff[k] = -cur_gradient[k];
      }

      cur_energy = best_energy;
    }
    else {
      INTERNAL_ERROR << " failed to get descent, sqr grad norm: " << sqr_grad_norm << std::endl;
      exit(1);
    }

  }

  return cur_energy;
}

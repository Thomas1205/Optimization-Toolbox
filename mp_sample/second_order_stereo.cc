/***** written by Thomas Schoenemann as a private person without employment, September 2011 ****/

#include "second_order_stereo.hh"

#include "factorMPBP.hh"
#include "factorDualOpt.hh"
#include "factorTRWS.hh"
#include "factorChainDualDecomp.hh"
#include "separatorDualOpt.hh"
#include "sepTRWS.hh"
#include "separatorChainDualDecomp.hh"


double disparity_energy(const Math3D::Tensor<double>& data_term, double so_weight, 
                        const Math2D::Matrix<uint>& disparity) {

  double energy = 0.0;
  
  const uint xDim = data_term.xDim();
  const uint yDim = data_term.yDim();

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      int cur_label = disparity(x,y);
      energy += data_term(x,y,cur_label);

      if (x > 0 && x+1 < xDim) {

        int diff1 = cur_label - int(disparity(x-1,y));
        int diff2 = int(disparity(x+1,y)) - cur_label;

        int so_diff = diff2 - diff1;

        if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
          ; //no cost
        else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
          energy += so_weight;
        else
          energy += 3*so_weight;
      }

      if (y > 0 && y+1 < yDim) {

        int diff1 = cur_label - int(disparity(x,y-1));
        int diff2 = int(disparity(x,y+1)) - cur_label;

        int so_diff = diff2 - diff1;
	  
        if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
          ; //no cost
        else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
          energy += so_weight;
        else
          energy += 3*so_weight;
      }
    }
  }


  return energy;
}

double second_order_stereo_msg_passing(const Math3D::Tensor<double>& data_term, 
                                       double so_weight, Math2D::Matrix<uint>& disparity,
                                       std::string mode) {

  std::cerr << "mode: " << mode << std::endl;

  assert(mode == "trws" || mode == "mplp" || mode == "msd" || mode == "chain-dd" || mode == "factor-dd"
         || mode == "sep-msd" || mode == "sep-trws" || mode == "sep-chain-dd" || mode == "mpbp");

  const int nLabels = data_term.zDim();

  Math3D::Tensor<float> trans_cost(nLabels,nLabels,nLabels,0.0);

  for (int l1=0; l1 < nLabels; l1++) {
    for (int l2=0; l2 < nLabels; l2++) {
      for (int l3=0; l3 < nLabels; l3++) {

        int diff1 = l2 - l1;
        int diff2 = l3 - l2;

        int so_diff = diff2 - diff1;

        if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
          trans_cost(l1,l2,l3) = 0.0; //no cost
        else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
          trans_cost(l1,l2,l3) = so_weight;
        else
          trans_cost(l1,l2,l3) = 3*so_weight;
      }
    }
  }

  const uint xDim = data_term.xDim();
  const uint yDim = data_term.yDim();

  disparity.resize(xDim,yDim);

  int do_fac = (mode == "mplp" || mode == "msd" || mode == "factor-dd") ? 1 : 0;
  int trws_fac = (mode == "trws") ? 1 : 0;
  int dd_fac = (mode == "chain-dd") ? 1 : 0;
  int sep_fac = (mode == "sep-msd") ? 1 : 0;
  int sep_trws_fac = (mode == "sep-trws") ? 1 : 0;
  int sep_dd_fac = (mode == "sep-chain-dd") ? 1 : 0;
  int bp_fac = (mode == "mpbp") ? 1 : 0;

  FactorDualOpt facDO(do_fac * xDim*yDim, do_fac * 2*xDim*yDim);
  SeparatorDualOptimization sepDO(sep_fac * xDim*yDim, sep_fac* (yDim*(xDim-1)+xDim*(yDim-1)), sep_fac * 2*xDim*yDim, false);
  CumFactorTRWS facTRWS(trws_fac * xDim*yDim, trws_fac * 2*xDim*yDim);
  FactorMPBP facMPBP(bp_fac *xDim*yDim, bp_fac*2*xDim*yDim);
  AllInclusiveSepCumTRWS sepTRWS(sep_trws_fac * xDim*yDim, sep_trws_fac* (yDim*(xDim-1)+xDim*(yDim-1)), sep_trws_fac * 2*xDim*yDim);
  FactorChainDualDecomposition dual_decomp(dd_fac * xDim*yDim, dd_fac * 2*xDim*yDim);
  SeparatorChainDualDecomposition sepDD(sep_dd_fac * xDim*yDim, sep_dd_fac* (yDim*(xDim-1)+xDim*(yDim-1)), sep_dd_fac * 2*xDim*yDim);

  Math1D::Vector<float> unaries(data_term.zDim());

  // construct energy
  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      for (uint l=0; l < data_term.zDim(); l++) {
        unaries[l] = data_term(x,y,l);
      }

      if (trws_fac == 1)
        facTRWS.add_var(unaries);
      else if (sep_trws_fac == 1)
        sepTRWS.add_var(unaries);
      else if (dd_fac == 1)
        dual_decomp.add_var(unaries);
      else if (do_fac == 1)
        facDO.add_var(unaries);
      else if (sep_dd_fac == 1)
        sepDD.add_var(unaries);
      else if (bp_fac == 1)
        facMPBP.add_var(unaries);
      else
        sepDO.add_var(unaries);
    }
  }

  uint vsep_offs = 10000;

  if (sep_fac == 1 || sep_trws_fac == 1 || sep_dd_fac == 1) {

    //add separators
#if 1
    std::cerr << "adding separators" << std::endl;
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim-1; x++) {
        if (sep_fac == 1)
          sepDO.add_separator(y*xDim+x, y*xDim+x+1);
        else if (sep_dd_fac == 1)
          sepDD.add_pair_separator(y*xDim+x, y*xDim+x+1);
        else
          sepTRWS.add_pair_separator(y*xDim+x, y*xDim+x+1);
      }
    }

    vsep_offs = yDim*(xDim-1);
    for (uint y=0; y < yDim-1; y++) {
      for (uint x=0; x < xDim; x++) {
        if (sep_fac == 1)
          sepDO.add_separator(y*xDim+x, (y+1)*xDim+x);
        else if (sep_dd_fac == 1)
          sepDD.add_pair_separator(y*xDim+x, (y+1)*xDim+x);
        else
          sepTRWS.add_pair_separator(y*xDim+x, (y+1)*xDim+x);        
      }
    }
    
    std::cerr << "added" << std::endl;
#endif
  }

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      uint cur_id = y*xDim+x;

      if (x > 0 && x+1 < xDim) {

        if (trws_fac == 1) {
          facTRWS.add_second_diff_factor(cur_id-1,cur_id,cur_id+1,so_weight);
        }
        else if (dd_fac == 1) {
          dual_decomp.add_second_diff_factor(cur_id-1,cur_id,cur_id+1,so_weight);
        }
        else if (do_fac == 1) {
          facDO.add_second_diff_factor(cur_id-1,cur_id,cur_id+1,so_weight);
        }
        else if (bp_fac == 1) {
          facMPBP.add_generic_ternary_factor(cur_id-1,cur_id,cur_id+1,trans_cost);
        }
        else {
           Storage1D<uint> sep(2);
           sep[0] = y*(xDim-1) + x-1;
           sep[1] = y*(xDim-1) + x;
          //Storage1D<uint> sep(0);

          if (sep_fac == 1)
            sepDO.add_generic_ternary_factor(cur_id-1,cur_id,cur_id+1,sep,trans_cost);
          else if (sep_dd_fac == 1)
            sepDD.add_ternary_factor(cur_id-1,cur_id,cur_id+1,sep,trans_cost);
          else
            sepTRWS.add_ternary_factor(cur_id-1,cur_id,cur_id+1,sep,trans_cost);
        }
      }

      if (y > 0 && y+1 < yDim) {

        if (trws_fac == 1) {
          facTRWS.add_second_diff_factor(cur_id-xDim,cur_id,cur_id+xDim,so_weight);
        }
        else if (dd_fac == 1) {
          dual_decomp.add_second_diff_factor(cur_id-xDim,cur_id,cur_id+xDim,so_weight);
        }
        else if (do_fac == 1) {
          facDO.add_second_diff_factor(cur_id-xDim,cur_id,cur_id+xDim,so_weight);
        }
        else if (bp_fac == 1) {
          facMPBP.add_generic_ternary_factor(cur_id-xDim,cur_id,cur_id+xDim,trans_cost);
        }
        else {
           Storage1D<uint> sep(2);
           sep[0] = vsep_offs + (y-1)*xDim + x;
           sep[1] = vsep_offs + y*xDim + x;
	   //Storage1D<uint> sep(0);

          if (sep_fac == 1)
            sepDO.add_generic_ternary_factor(cur_id-xDim,cur_id,cur_id+xDim,sep,trans_cost);
          else if (sep_dd_fac == 1)
            sepDD.add_ternary_factor(cur_id-xDim,cur_id,cur_id+xDim,sep,trans_cost);
          else
            sepTRWS.add_ternary_factor(cur_id-xDim,cur_id,cur_id+xDim,sep,trans_cost);
        }
      }
    }
  }  

  double lower_bound = 0.0;

  if (trws_fac == 1) {
    lower_bound = facTRWS.optimize(250);
    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = facTRWS.labeling()[k];
  }
  else if (sep_trws_fac == 1) {
    lower_bound = sepTRWS.optimize(250);
    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = sepTRWS.labeling()[k];
  }
  else if (dd_fac == 1) {
    lower_bound = dual_decomp.optimize(500,20.0); //with 1/t
    //lower_bound = dual_decomp.optimize(500,0.5); //with primal-dual
    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = dual_decomp.labeling()[k];
  }
  else  if (do_fac == 1) {
    if (mode == "mplp")
      lower_bound = facDO.dual_bca(500,DUAL_BCA_MODE_MPLP,true,false);
    else if (mode == "factor-dd")
      lower_bound = facDO.subgradient_opt(1000,0.001);
    else
      lower_bound = facDO.dual_bca(500,DUAL_BCA_MODE_MSD,true,false);

    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = facDO.labeling()[k];
  }
  else if (bp_fac == 1) {

    facMPBP.mpbp(500);
    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = facMPBP.labeling()[k];    
  }
  else if (sep_dd_fac == 1) {

    lower_bound = sepDD.optimize(1500,2.0);

    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = sepDD.labeling()[k];    
  }
  else {
    lower_bound = sepDO.optimize(500,DUAL_BCA_MODE_MSD,false);

    for (uint k=0; k < disparity.size(); k++)
      disparity.direct_access(k) = sepDO.labeling()[k];    
  }

  std::cerr << "lower bound: " << lower_bound << std::endl;
  std::cerr << "discrete energy: " << disparity_energy(data_term,so_weight,disparity) << std::endl;

  return lower_bound;
}

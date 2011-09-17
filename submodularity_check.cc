/************** written by Thomas Schoenemann at the University of Bonn, 2007. Slightly adapted later in Lund ********/
/************** usage granted by Daniel Cremers *******/


#include "submodularity_check.hh"


bool is_submodular(BinFunction& bin_func) {

  int order = bin_func.order();

  assert(order < 8*((int) sizeof(int)));

  int bound = 1;
  bound = bound << order;

  for (int lower = 0; lower < (bound -1); lower++) {
    double l_val = bin_func.value(lower);
    for (int upper = lower+1; upper < bound; upper++) {

      double u_val = bin_func.value(upper);

      int intersect = lower & upper;
      int uni = lower | upper;

      double int_val = bin_func.value(intersect);
      double union_val = bin_func.value(uni);

      if ((l_val + u_val) < (int_val + union_val)) {
	std::cerr << "function is not submodular: f(" << lower << ") = " << l_val <<  ", f(" << upper << ") = "
		  << u_val << ", f(" << intersect << ") = " << int_val << ", f(" << uni << ") = " << union_val << std::endl;
	return false;
      }


    }
  }

  return true;
}

bool satisfies_freedman_conditions(BinFunction& bin_func) {

  int order = bin_func.order();

  /**** checking binary cliques ****/
  for (int c1 = 0; c1 < order; c1++) {
    for (int c2 = c1+1; c2 < order; c2++) {

      int arg = (1 << c1) + (1 << c2);
      double val = bin_func.value(arg) + bin_func.value(0)
	- bin_func.value(1 << c1) - bin_func.value(1 << c2);

      if (val > 0.0) {
	std::cerr << "freedman condition violated for pairwise clique " << c1 << "," << c2 
		  << ": value " << val << " is positive" << std::endl;
	return false;
      }
    }
  }

  int bound = 1;
  bound = bound << order;

  // 7 is the lower clique with more than two elements
  for (int clique_val = 7; clique_val < bound; clique_val++) {

    double sum = 0.0;
    int nBitsSet = 0;
    for (int i=0; i < order; i++) {
      if (clique_val & (1 << i))
	nBitsSet++;
    }

    if (nBitsSet >= 3) {
      //std::cerr << "clique_val: " << clique_val << std::endl;

      int init_sign = (nBitsSet % 2 == 0) ? 1 : -1; 

      for (int val =0; val < bound; val++) {
		
	if ((val & clique_val) == val) {
	  int sign = init_sign;
		    
	  for (int i=0; i < order; i++) {
	    if (val & (1 << i))
	      sign *= -1;
	  }
		    
	  sum += sign * bin_func.value(val);
	  //std::cerr << "intermediate result: " << result << std::endl;
	}
      }

      //std::cerr << "result: " << result << std::endl;
      if (sum > 0.0) {
	std::cerr << "function violates Freedman-condition for order " << nBitsSet 
		  << ": clique-val " << clique_val << " has positive sum: " << sum << std::endl;
	return false;
      }
    }
	
  }
    
  return true;
}

bool is_supermodular(BinFunction& bin_func) {

  int order = bin_func.order();

  assert(order < 8*((int) sizeof(int)));

  int bound = 1;
  bound = bound << order;

  for (int lower = 0; lower < (bound -1); lower++) {
    double l_val = bin_func.value(lower);
    for (int upper = lower+1; upper < bound; upper++) {

      double u_val = bin_func.value(upper);

      int intersect = lower & upper;
      int uni = lower | upper;

      double int_val = bin_func.value(intersect);
      double union_val = bin_func.value(uni);

      if ((l_val + u_val) > (int_val + union_val)) {
	std::cerr << "function is not supermodular: f(" << lower << ") = " << l_val <<  ", f(" << upper << ") = "
		  << u_val << ", f(" << intersect << ") = " << int_val << ", f(" << uni << ") = " << union_val << std::endl;
	return false;
      }


    }
  }

  return true;
}


bool is_posimodular(BinFunction& bin_func) {

  int order = bin_func.order();

  assert(order < 8*((int) sizeof(int)));

  int bound = 1;
  bound = bound << order;

  for (int lower = 0; lower < (bound -1); lower++) {
    double l_val = bin_func.value(lower);
    for (int upper = lower+1; upper < bound; upper++) {

      double u_val = bin_func.value(upper);

      int fwd_diff = lower;

      for (int k=0; k < order; k++) {
	int bit = 1 << k;

	if ( (fwd_diff & bit) && (upper & bit))
	  fwd_diff -= bit;
      }

      int back_diff = upper;

      for (int k=0; k < order; k++) {
	int bit = 1 << k;

	if ( (back_diff & bit) && (lower & bit))
	  back_diff -= bit;
      }

      double fwd_val = bin_func.value(fwd_diff);
      double back_val = bin_func.value(back_diff);

      if ((l_val + u_val) < (fwd_val + back_val)) {
	std::cerr << "function is not posimodular: f(" << lower << ") = " << l_val <<  ", f(" << upper << ") = "
		  << u_val << ", f(" << fwd_diff << ") = " << fwd_val << ", f(" << back_diff << ") = " << back_val << std::endl;
	return false;
      }


    }
  }

  return true;
}


bool is_negimodular(BinFunction& bin_func) {

  int order = bin_func.order();

  assert(order < 8*((int) sizeof(int)));

  int bound = 1;
  bound = bound << order;

  for (int lower = 0; lower < (bound -1); lower++) {
    double l_val = bin_func.value(lower);
    for (int upper = lower+1; upper < bound; upper++) {

      double u_val = bin_func.value(upper);

      int fwd_diff = lower;

      for (int k=0; k < order; k++) {
	int bit = 1 << k;

	if ( (fwd_diff & bit) && (upper & bit))
	  fwd_diff -= bit;
      }

      int back_diff = upper;

      for (int k=0; k < order; k++) {
	int bit = 1 << k;

	if ( (back_diff & bit) && (lower & bit))
	  back_diff -= bit;
      }

      double fwd_val = bin_func.value(fwd_diff);
      double back_val = bin_func.value(back_diff);

      if ((l_val + u_val) > (fwd_val + back_val)) {
	std::cerr << "function is not negimodular: f(" << lower << ") = " << l_val <<  ", f(" << upper << ") = "
		  << u_val << ", f(" << fwd_diff << ") = " << fwd_val << ", f(" << back_diff << ") = " << back_val << std::endl;
	return false;
      }
    }
  }

  return true;
}


/************** written by Thomas Schoenemann at the University of Bonn, 2007. Slightly adapted later in Lund ********/
/************** usage granted by Daniel Cremers *******/

#ifndef SUBMODULARITY_CHECK_HH
#define SUBMODULARITY_CHECK_HH

#include "function.hh"

bool is_submodular(BinFunction& bin_func);

bool satisfies_freedman_conditions(BinFunction& bin_func);

bool is_supermodular(BinFunction& bin_func);

bool is_posimodular(BinFunction& bin_func);

bool is_negimodular(BinFunction& bin_func);

#endif

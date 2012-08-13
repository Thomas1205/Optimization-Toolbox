/************** written by Thomas Schoenemann at the University of Bonn, 2007. Slightly adapted later in Lund ********/
/************** usage granted by Daniel Cremers *******/


#include "function.hh"
#include "application.hh"
#include "stringprocessing.hh"
#include "submodularity_check.hh"

#include <set>

int main(int argc, char* argv[]) {

  if (argc == 1 || strings_equal(argv[1],"-h")) {

    std::cerr << "usage: " << argv[0] << std::endl
	      << "-order <uint>" << std::endl
	      << "-func <filename>" << std::endl
	      << "-invert <comma separated list of bit positions> : also check submodularity when the given bits are inverted" << std::endl;
  }

  const int nParam = 3;

  ParamDescr p_descr[nParam] = {{"-order",mandWithValue,0,""},{"-func",mandInFilename,0,""},{"-invert",optWithValue,0,""}};
    
  Application app(argc,argv,p_descr,nParam);

  int order = convert<int>(app.getParam("-order"));

  int nValues = 1 << order;

  ExplicitBinFunction bin_func(order);
  bin_func.read(app.getParam("-func"));

  std::cerr << std::endl << std::endl << "*************** checking for sub-modularity" << std::endl;

  bool is_submod = is_submodular(bin_func);
  if (is_submod)
    std::cerr << "the function is submodular" << std::endl;
  else
    std::cerr << "the function is not submodular" << std::endl;
    
  if (order > 2) {

    bool sf = satisfies_freedman_conditions(bin_func);
    if (sf)
      std::cerr << "satiesfies Freedman conditions => graph-representable" << std::endl;
    else
      std::cerr << "does not satisfy the Freedman conditions" << std::endl;

    std::cerr << std::endl << "********* checking for supermodularity" << std::endl;
    bool is_supermod = is_supermodular(bin_func);
    if (is_supermod)
      std::cerr << "the function is supermodular" << std::endl;
    else
      std::cerr << "the function is not supermodular" << std::endl;

    std::cerr << std::endl << "********** checking for posi-modularity" << std::endl;
    bool is_posmod = is_posimodular(bin_func);
    if (is_posmod) 
      std::cerr << "the function is posi-modular" << std::endl;
    else
      std::cerr << "the function is not posi-modular" << std::endl;
	
    std::cerr << std::endl << "********** checking for negi-modularity" << std::endl;
    bool is_negmod = is_negimodular(bin_func);
    if (is_negmod) 
      std::cerr << "the function is negi-modular" << std::endl;
    else
      std::cerr << "the function is not negi-modular" << std::endl;

    std::cerr << std::endl;
  }

  if (app.is_set("-invert")) {

    std::string s = app.getParam("-invert");
    std::vector<std::string> num_strings;
    tokenize(s,num_strings,',');

    int bit_string  = 0;

    for (size_t k=0; k < num_strings.size(); k++) {
      unsigned int pos = convert<unsigned int>(num_strings[k]);

      bit_string |= 1 << pos;
    }

    std::cerr << std::endl << "bit_string: " << bit_string << std::endl;

    ExplicitBinFunction func2(order);
    for (int i=0; i < nValues; i++) {

      int j = i ^ bit_string; //xor = invert positions

      // 	    std::cerr << "i: " << i << std::endl;
      // 	    std::cerr << "j: " << j << std::endl;

      func2.set_val(i,bin_func.value(j));
    }

    std::cerr << std::endl << "**********checking function with inverted bits" << std::endl;
    bool is_submod = is_submodular(func2);
    if (is_submod)
      std::cerr << "the function is submodular" << std::endl;
    else
      std::cerr << "the function is not submodular" << std::endl;

    if (order > 2) {
      std::cerr << std::endl << "********checking function with inverted bits for supermodularity" << std::endl;
	    
      bool is_supermod = is_supermodular(func2);
      if (is_supermod)
	std::cerr << "the function is supermodular" << std::endl;
      else
	std::cerr << "the function is not supermodular" << std::endl;

      std::cerr << std::endl << "********** checking function with inverted bits for posi-modularity" << std::endl;
      bool is_posmod = is_posimodular(func2);
      if (is_posmod) 
	std::cerr << "the function is posi-modular" << std::endl;
      else
	std::cerr << "the function is not posi-modular" << std::endl;
	    
      std::cerr << std::endl << "********** checking function with inverted bits for negi-modularity" << std::endl;
      bool is_negmod = is_negimodular(func2);
      if (is_negmod) 
	std::cerr << "the function is negi-modular" << std::endl;
      else
	std::cerr << "the function is not negi-modular" << std::endl;


      //DEBUG
      // 	    std::ofstream out("supermodular_cornerpoints.txt");
      // 	    out << "negative of the supermodular cornerpoint function: " << std::endl << std::endl << std::endl;
      // 	    out << "ijkl    value    value+2" << std::endl;
      // 	    for (int i=0; i < nValues; i++) {
      // 		out << ((i & 8) >> 3) << ((i & 4) >> 2) << ((i & 2) >> 1) << (i & 1) << "     " 
      // 		    << neg_func.value(i) << "         "  << (neg_func.value(i) + 2.0) << std::endl;
      // 	    }
      // 	    out.close();
      //END_DEBUG
     
    }
	
  }

}

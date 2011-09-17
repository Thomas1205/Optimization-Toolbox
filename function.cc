/************** written by Thomas Schoenemann at the University of Bonn, 2007. Slightly adapted later in Lund ********/
/************** usage granted by Daniel Cremers *******/

#include "function.hh"
#include <fstream>

BinFunction::BinFunction(int order) : order_(order) {}

BinFunction::~BinFunction() {}

/*virtual*/ int BinFunction::order() {
  return order_;
}

/*virtual*/ double BinFunction::value(const Math1D::Vector<int>& arg) { // first entry in arg means lowest bit, etc...

  /*** compute integer argument ***/
  int int_arg = 0;

  int base = 1;
  for (int i=0; i < (int) arg.size(); i++) {

    if (i != 0)
      base *= 2;
    assert(arg[i] == 0 || arg[i] == 1);
    int_arg += arg[i] * base;
  }

  return value(int_arg);

}


ExplicitBinFunction::ExplicitBinFunction(int order) : BinFunction(order) {

  int size = 1;
  size = size << order;
  value_.resize(size);
  value_.set_constant(0.0);
}

//read function values from file, file contains lines of form bit1...bitn <value>
void ExplicitBinFunction::read(std::string filename) {

  std::ifstream file(filename.c_str());

  std::string bit_pattern = "";
  double val;

  while (file >> bit_pattern >> val) {

    assert(((int) bit_pattern.size()) == order_);

    /*** convert bit-pattern to integer ***/
    int base = 1;
    int arg = 0;
    for (int i=0; i < order_; i++) {
      if (i != 0)
	base *= 2;

      if (bit_pattern[i] == '1')
	arg += base;
      else
	assert(bit_pattern[i] == '0');
    }
    value_[arg] = val;	
  }
}

/*virtual*/ double ExplicitBinFunction::value(int arg) {

  assert(arg < (int) value_.size());
  return value_[arg];
}

void ExplicitBinFunction::set_val(int arg, double value) {
  value_[arg] = value;
}

double ExplicitBinFunction::multilinear_coefficient(std::vector<size_t> one_bits) {
    
  int clique_order = one_bits.size();
#if 0
  double result = 0.0;

  //std::cerr << "max: " << (1 <<clique_order) << std::endl;

  for (int clique_arg = 0; clique_arg < (1 << clique_order); clique_arg++) {

    int sign = (clique_order % 2 == 0) ? 1 : -1;
    int arg = 0;

    for (int i=0; i < clique_order; i++) {
      if (clique_arg & (1 << i)) {
	sign *= -1;
	arg += 1 << one_bits[i];
      }
    }

    result += sign * value_[arg];
  }

  return result;
#else

  int one_arg = 0;
  for (int k=0; k < (int) one_bits.size(); k++)
    one_arg += 1 << one_bits[k];

  double result = value_[one_arg];
  for (int arg=0; arg < (1 << order_); arg++) {
    if (arg != one_arg && (arg & one_arg) == arg) {

      std::vector<size_t> v;
      for (int i = 0; i < order_; i++) {
	if (arg & (1 << i))
	  v.push_back(i);
      }
      result -= multilinear_coefficient(v);
    }
  }

  return result;
#endif
}

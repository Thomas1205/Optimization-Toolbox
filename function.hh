/************** written by Thomas Schoenemann at the University of Bonn, 2007. Slightly adapted later in Lund ********/
/************** usage granted by Daniel Cremers *******/


#ifndef FUNCTION_HH
#define FUNCTION_HH

#include "vector.hh"

#include <vector>

// function of binary variables, can assume any real value
/*abstract */ class BinFunction {
public:
  BinFunction(int order);

  ~BinFunction();

  virtual int order();

  virtual double value(const Math1D::Vector<int>& arg); // first entry in arg means lowest bit, etc...

  virtual double value(int arg) = 0; // bit pattern of arg defines values of variables

protected:
  int order_;
};

class ExplicitBinFunction : public BinFunction {
public:
  ExplicitBinFunction(int order);

  //read function values from file, file contains lines of form bit1...bitn <value>
  void read(std::string filename);

  void set_val(int arg, double value);

  virtual double value(int arg);

  double multilinear_coefficient(std::vector<size_t> one_bits);

private:
  Math1D::Vector<double> value_;
};


#endif

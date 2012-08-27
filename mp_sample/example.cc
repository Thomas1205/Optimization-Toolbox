
//this example was developed in collaboration with Vladimir Kolmogorov

#include "factorMPBP.hh"
#include "factorDualOpt.hh"
#include "factorTRWS.hh"
#include "factorChainDualDecomp.hh"
#include "separatorDualOpt.hh"
#include "sepTRWS.hh"
#include "separatorChainDualDecomp.hh"

// Example: minimize function   f(x,y,z) = x + 2y + 3z + |x-y| - |x+y-2z|   where
//   x \in {0,1}
//   y \in {0,1,2}
//   z \in {0,1}

const int node_num = 3;

const int K[node_num] = { 2, 3, 2 }; // number of labels for nodes 0,1,2
float f_unary(int i, int label) { return (i+1)*label; }

const int factor_num = 2;
float f0(float x, float y) { return fabs(x-y); } // first factor       |x-y|
float f1(float x, float y, float z) { return -fabs(x+y-2*z); } // second factor       -|x+y-2z|







/////////////////////////////////////////////////////////////////////////////////

void Run_facMSD(){

  std::cerr << "***** MSD ******" << std::endl;

  int index, i, a;

  FactorDualOpt facMSD(node_num, factor_num);

  for (i=0; i<node_num; i++) {

    Math1D::Vector<float> unaries(K[i]); 

    for (a=0; a<K[i]; a++)
      {
        unaries[a] = f_unary(i, a); 
      }
    facMSD.add_var(unaries);
  }

  // add term f0 
  Math2D::Matrix<float> cost0(K[0],K[1]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      cost0(x,y) = f0(x,y);

  facMSD.add_generic_binary_factor(0,1,cost0);

  // add term f1
  Math3D::Tensor<float> cost1(K[0],K[1],K[2]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      for (int z=0; z < K[2]; z++)
        cost1(x,y,z) = f1(x,y,z);

  facMSD.add_generic_ternary_factor(0,1,2,cost1);

  // call optimizer
  double bound = facMSD.dual_bca(10,DUAL_BCA_MODE_MSD);

  //print solution and bound - TODO
  const Math1D::Vector<uint>& solution = facMSD.labeling();
}

void Run_facMPLP(){

  std::cerr << "***** MPLP ******" << std::endl;

  int index, i, a;

  FactorDualOpt facMPLP(node_num, factor_num);

  for (i=0; i<node_num; i++) {

    Math1D::Vector<float> unaries(K[i]); 

    for (a=0; a<K[i]; a++)
      {
        unaries[a] = f_unary(i, a); 
      }
    facMPLP.add_var(unaries);
  }

  // add term f0 
  Math2D::Matrix<float> cost0(K[0],K[1]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      cost0(x,y) = f0(x,y);

  facMPLP.add_generic_binary_factor(0,1,cost0);

  // add term f1
  Math3D::Tensor<float> cost1(K[0],K[1],K[2]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      for (int z=0; z < K[2]; z++)
        cost1(x,y,z) = f1(x,y,z);

  facMPLP.add_generic_ternary_factor(0,1,2,cost1);

  // call optimizer
  double bound = facMPLP.dual_bca(10,DUAL_BCA_MODE_MPLP);

  //print solution and bound - TODO
  const Math1D::Vector<uint>& solution = facMPLP.labeling();
}


void Run_sepMSD(){

  //TODO
}

void Run_facTRWS()
{

  std::cerr << "***** TRWS ******" << std::endl;

  int index, i, a;

  CumFactorTRWS facTRWS(node_num, factor_num);

  for (i=0; i<node_num; i++) {

    Math1D::Vector<float> unaries(K[i]); 

    for (a=0; a<K[i]; a++)
      {
        unaries[a] = f_unary(i, a); 
      }
    facTRWS.add_var(unaries);
  }

  // add term f0 
  Math2D::Matrix<float> cost0(K[0],K[1]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      cost0(x,y) = f0(x,y);

  facTRWS.add_binary_factor(0,1,cost0);

  // add term f1
  Math3D::Tensor<float> cost1(K[0],K[1],K[2]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      for (int z=0; z < K[2]; z++)
        cost1(x,y,z) = f1(x,y,z);

  facTRWS.add_ternary_factor(0,1,2,cost1);

  // call optimizer
  double bound = facTRWS.optimize(10);

  //print solution and bound - TODO
  const Math1D::Vector<uint>& solution = facTRWS.labeling();
}

void Run_facMPBP(){

  std::cerr << "***** Belief Propagation ******" << std::endl;

  int index, i, a;

  FactorMPBP facMPBP(node_num, factor_num);

  for (i=0; i<node_num; i++) {

    Math1D::Vector<float> unaries(K[i]); 

    for (a=0; a<K[i]; a++)
      {
        unaries[a] = f_unary(i, a); 
      }
    facMPBP.add_var(unaries);
  }

  // add term f0 
  Math2D::Matrix<float> cost0(K[0],K[1]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      cost0(x,y) = f0(x,y);

  facMPBP.add_generic_binary_factor(0,1,cost0);

  // add term f1
  Math3D::Tensor<float> cost1(K[0],K[1],K[2]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      for (int z=0; z < K[2]; z++)
        cost1(x,y,z) = f1(x,y,z);

  facMPBP.add_generic_ternary_factor(0,1,2,cost1);

  // call optimizer
  facMPBP.mpbp(10);

  //print solution and bound - TODO
  const Math1D::Vector<uint>& solution = facMPBP.labeling();
}

void Run_sepTRWS(){}

void Run_dual_decomp(){


  std::cerr << "***** DD/SG ******" << std::endl;

  int index, i, a;

  FactorChainDualDecomposition facDD(node_num, factor_num);

  for (i=0; i<node_num; i++) {

    Math1D::Vector<float> unaries(K[i]); 

    for (a=0; a<K[i]; a++)
      {
        unaries[a] = f_unary(i, a); 
      }
    facDD.add_var(unaries);
  }

  // add term f0 
  Math2D::Matrix<float> cost0(K[0],K[1]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      cost0(x,y) = f0(x,y);

  facDD.add_binary_factor(0,1,cost0);

  // add term f1
  Math3D::Tensor<float> cost1(K[0],K[1],K[2]);
  for (int x=0; x < K[0]; x++)
    for (int y=0; y < K[1]; y++)
      for (int z=0; z < K[2]; z++)
        cost1(x,y,z) = f1(x,y,z);

  facDD.add_ternary_factor(0,1,2,cost1);

  // call optimizer
  double bound = facDD.optimize(10,1.0);

  //print solution and bound - TODO
  const Math1D::Vector<uint>& solution = facDD.labeling();
}

void Run_sepDD(){}


int main()
{
  Run_facMSD(); // one-line description of what it does - TODO
  Run_facMPLP();
  Run_sepMSD();
  Run_facTRWS();
  Run_facMPBP();
  Run_sepTRWS();
  Run_dual_decomp();
  Run_sepDD();

  return 0;
}

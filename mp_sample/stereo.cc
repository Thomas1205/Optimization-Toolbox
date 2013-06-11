/**** written by Thomas Schoenemann as an employee of Lund University, January 2010 ****/

#include "application.hh"
#include "grayimage.hh"
#include "timing.hh"

#include "second_order_stereo.hh"

int main(int argc, char** argv) {


  if (argc == 1 || (argc == 2 && strings_equal(argv[1],"-h"))) {

    std::cerr << "USAGE: " << argv[0] << std::endl
              << " -first <pgm> : left stereo image" << std::endl
              << " -second <pgm> : right stereo image" << std::endl
              << " -lambda <double> : smoothness weight" << std::endl
              << " -method ( chain-dd | factor-dd | trws | mplp | msd | sep-msd | sep-trws | sep-chain-dd | mpbp):" 
              << "     the computational method to use, default: trws " << std::endl
              << " -max-disp <uint> : maximal allowed disparity" << std::endl 
              << " -o <filename(pgm)> : image where the disparity map is written" << std::endl;
  }

  const int nParams = 6;
  ParamDescr  params[nParams] = {{"-first",mandInFilename,0,""},{"-second",mandInFilename,0,""},
                                 {"-lambda",optWithValue,1,"1.0"},{"-o",mandOutFilename,0,""},
                                 {"-max-disp",mandWithValue,0,""},{"-method",optWithValue,1,"trws"}};

  Application app(argc,argv,params,nParams);

  Math2D::NamedGrayImage<float> first(app.getParam("-first"),MAKENAME(first));
  Math2D::NamedGrayImage<float> second(app.getParam("-second"),MAKENAME(second));

  uint max_disp = convert<uint>(app.getParam("-max-disp"));

  double lambda = convert<double>(app.getParam("-lambda"));

  uint xDim = first.xDim();
  uint yDim = first.yDim();

  assert(second.xDim() == xDim);
  assert(second.yDim() == yDim);

  uint nLabels = max_disp+1;
  
  Math3D::NamedTensor<double> data_term(xDim,yDim,nLabels,MAKENAME(data_term));

  double max_cost = -MAX_DOUBLE;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      double min_cost = MAX_DOUBLE;

      for (uint l=0; l < nLabels; l++) {

        int x_target = ((int) x) - ((int) l);

        double cost = 0.0;

        if (x_target >= 0) {

          float diff = first(x,y) - second(x_target,y);

          cost = std::abs(diff);
        }

        data_term(x,y,l) = cost;

        if (cost > max_cost)
          max_cost = cost;

        if (cost < min_cost)
          min_cost = cost;
      }
    }
  }

  std::cerr << "maximal data term: " << max_cost << std::endl;

  Math2D::NamedMatrix<uint> disparity(xDim,yDim,0,MAKENAME(disparity));

  std::string method_string = app.getParam("-method"); 

  std::clock_t tStartComp, tEndComp;
  
  tStartComp = std::clock();

  std::string method = "trws";
  if (method_string == "mplp" || method_string == "msd" || method_string == "chain-dd" 
      || method_string == "factor-dd" || method_string == "sep-msd" || method_string == "sep-trws" 
      || method_string == "sep-chain-dd" || method_string == "mpbp")
    method = method_string;
  
  second_order_stereo_msg_passing(data_term, lambda, disparity, method);
  
  tEndComp = std::clock();


  std::cerr << "computation took " << diff_seconds(tEndComp,tStartComp) << " seconds" << std::endl;

  uint factor = 1;
  while (max_disp*factor < 255) 
    factor++;
  factor--;

  disparity *= factor;
  
  disparity.savePGM(app.getParam("-o"),255);
}

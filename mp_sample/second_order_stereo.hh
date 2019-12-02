/***** written by Thomas Schoenemann as a private person without employment, September 2011 ****/

#ifndef SECOND_ORDER_STEREO_HH
#define SECOND_ORDER_STEREO_HH

#include "tensor.hh"
#include "matrix.hh"

double second_order_stereo_msg_passing(const Math3D::Tensor<double>& data_term, 
                                       double so_weight, Math2D::Matrix<uint>& disparity,
                                       std::string mode = "trws");

#endif

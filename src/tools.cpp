#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
   rmse << 0,0,0,0;

   for(int i=0; i < estimations.size(); ++i)
   {
     VectorXd temp = (estimations[i] - ground_truth[i]);
     temp = temp.array()*temp.array();
     rmse += temp;

   }

   rmse = rmse/estimations.size();
   rmse = rmse.array().sqrt();

   return rmse;
}

#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.   
   */
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
      std::cout << "invalid inputs for RMSE" << std::endl;
      return rmse;
  }

  // accumulate squared residuals
  vector<VectorXd> residuals2;
  for (int i=0; i < estimations.size(); ++i) {
    VectorXd r = estimations[i] - ground_truth[i];
    VectorXd r2 = (r.array())*(r.array());
    residuals2.push_back(r2);
  }

  // calculate the mean
  VectorXd sum(4);
  sum << 0,0,0,0;
  for (int i=0; i < residuals2.size(); ++i) {
    sum = sum + residuals2[i];
  }
  VectorXd mean = sum.array()/(residuals2.size());
  
  // calculate the squared root
  rmse = mean.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  if (px == 0 && py == 0){
      std::cout << "division by zero in Jacobian" << std::endl;
      return Hj;
  }
  // compute the Jacobian matrix
  else{
      float p2 = px*px + py*py;
      float ps = sqrt(p2);
      float vp = vx*py - vy*px;
      Hj << px/ps, py/ps, 0, 0,
            -py/p2, px/p2, 0, 0,
            (py*vp)/(p2*ps), (px*-vp)/(p2*ps), px/ps, py/ps;
  }

  return Hj;
}

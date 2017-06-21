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

    // TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	if (estimations.size() == 0) {
	    cerr << "Nothing to estimate. Vector empty" << endl;
	    return rmse;
	}
	
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
	if (estimations.size() != ground_truth.size()) {
	    cerr << "The vector sizes differ!" << endl;
	}

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        for(int j=0; j<4; j++) {
            float d = estimations[i][j] - ground_truth[i][j];
		    rmse[j] += d*d;
        }
        
	}

	//calculate the mean
	// ... your code here
	rmse = rmse.array()/estimations.size();
	

	//calculate the squared root
	// ... your code here
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 

	//check division by zero
	float sq_sum = px*px + py*py;
	if (sq_sum == 0) {
	    std::cerr << "division by zero, can't compute Jacobian Matrix" << endl;
	    return Hj;
	}
	
	//compute the Jacobian matrix
	float sq_root = sqrt(sq_sum);
	float th_by_tw = pow(sq_sum, 3/2);
	Hj << px/sq_root, py/sq_root, 0, 0, 
	        -py/sq_sum, px/sq_sum, 0, 0,
	        py*(vx*py - vy*px)/th_by_tw, px*(vy*px - vx*py)/th_by_tw, px/sq_root, py/sq_root;
	        

	return Hj;
}

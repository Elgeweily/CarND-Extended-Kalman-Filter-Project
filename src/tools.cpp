#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd accum(4);
	accum << 0, 0, 0, 0;

	VectorXd res(4);
	VectorXd ressq(4);
	VectorXd mean(4);
	VectorXd rmse(4);

	if (estimations.size() == 0) {
		cout << "error...no estimations" << endl;
	}

	if (estimations.size() != ground_truth.size()) {
		cout << "error...estimations and ground truth must have same size" << endl;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		res = estimations[i] - ground_truth[i];
		ressq = res.array()*res.array();
		accum = accum + ressq;
	}

	mean = accum.array() / estimations.size();

	rmse = mean.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj_(3, 4);
	Hj_.setZero();

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = pow(px, 2) + pow(py, 2);
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (abs(c1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj_;
	}

	//compute the Jacobian matrix
	Hj_ << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;

	return Hj_;
}

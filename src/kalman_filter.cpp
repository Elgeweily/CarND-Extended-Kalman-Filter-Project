#define M_pi 3.1415926535897

#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
	MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

	VectorXd y = z - H_ * x_;
	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

	//new estimate
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

	// define measurement function
	VectorXd h_func = VectorXd(3);
	h_func(0) = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
	h_func(1) = atan2(x_(1), x_(0));
	h_func(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / sqrt(pow(x_(0), 2) + pow(x_(1), 2));

	VectorXd y = z - h_func;

	while (y(1) > M_pi) {
		y(1) -= 2 * M_pi;
	}
	while (y(1) < -M_pi) {
		y(1) += 2 * M_pi;
	}

	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

	//new estimate
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

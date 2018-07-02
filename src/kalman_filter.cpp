#define M_pi 3.1415926535897

#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	// Kalman filter equations
	VectorXd y = z - H_ * x_;
	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

	//new estimate
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// define measurement function
	VectorXd h_func = VectorXd(3);
	h_func(0) = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
	h_func(1) = atan2(x_(1), x_(0));
	h_func(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / sqrt(pow(x_(0), 2) + pow(x_(1), 2));

	// extended Kalman filter equations
	VectorXd y = z - h_func;

	// limit bearing angle between -pi & pi
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

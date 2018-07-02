#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices and vectors
  ekf_.x_ = VectorXd(4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  H_laser = MatrixXd(2, 4);
  R_laser = MatrixXd(2, 2);
  R_radar = MatrixXd(3, 3);
  ekf_.Q_ = MatrixXd(4, 4);

  // object covariance matrix
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;

  // measurement matrix - laser
  H_laser << 1, 0, 0, 0,
	  0, 1, 0, 0;

  // measurement covariance matrix - laser
  R_laser << 0.0225, 0,
	  0, 0.0225;

  // measurement covariance matrix - radar
  R_radar << 0.09, 0, 0,
	  0, 0.0009, 0,
	  0, 0, 0.09;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float range = measurement_pack.raw_measurements_(0);
      float bearing = measurement_pack.raw_measurements_(1);
        
      ekf_.x_ << range * cos(bearing), range * sin(bearing), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */  
	  ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }
	
	// update timestamp
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/


  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  // set the F matrix
  ekf_.F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  // calculate powers of dt  
  float dt_2 = pow(dt, 2);
  float dt_3 = pow(dt, 3);
  float dt_4 = pow(dt, 4);

  // noise variance values
  float noise_ax = 9;
  float noise_ay = 9;

  // set the process covariance matrix Q
  ekf_.Q_ << (dt_4 / 4) * noise_ax, 0, (dt_3 / 2) * noise_ax, 0,
	  0, (dt_4 / 4) * noise_ay, 0, (dt_3 / 2) * noise_ay,
	  (dt_3 / 2) * noise_ax, 0, dt_2 * noise_ax, 0,
	  0, (dt_3 / 2) * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  z = measurement_pack.raw_measurements_;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // radar updates
		  Hj = tools.CalculateJacobian(ekf_.x_);
		  if (!Hj.isZero(0)) {
			  ekf_.H_ = Hj;
			  ekf_.R_ = R_radar;
			  ekf_.UpdateEKF(z);
		  }
		  else {
			  // division by zero, skip update step
		  }

  } else {
      // laser updates
	  ekf_.H_ = H_laser;
	  ekf_.R_ = R_laser;
	  ekf_.Update(z);

	  // print the output
	  cout << "x_ = " << ekf_.x_ << endl;
	  cout << "P_ = " << ekf_.P_ << endl;
  }
}

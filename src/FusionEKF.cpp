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

  // initializing matrices
  x = VectorXd(4);
  P = MatrixXd(4, 4);
  F = MatrixXd(4, 4);
  H_laser = MatrixXd(2, 4);
  Hj = MatrixXd(3, 4);
  R_laser = MatrixXd(2, 2);
  R_radar = MatrixXd(3, 3);
  Q = MatrixXd(4, 4);

  //object covariance matrix
  P << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;

  //measurement matrix - laser
  H_laser << 1, 0, 0, 0,
	  0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser << 0.0225, 0,
	  0, 0.0225;

  //measurement covariance matrix - radar
  R_radar << 0.09, 0, 0,
	  0, 0.0009, 0,
	  0, 0, 0.09;

  noise_ax = 9;
  noise_ay = 9;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    //ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float range = measurement_pack.raw_measurements_(0);
      float bearing = measurement_pack.raw_measurements_(1);
        
      x << range * cos(bearing), range * sin(bearing), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */  
      x << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }

	ekf_.Init(x, P, F, H_laser, R_laser, Q);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //Set the F matrix
  F << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  //Set the process covariance matrix Q
  Q << ((pow(dt, 4) / 4.0) * noise_ax), 0, ((pow(dt, 3) / 2.0) * noise_ax), 0,
	  0, ((pow(dt, 4) / 4.0) * noise_ay), 0, ((pow(dt, 3) / 2.0) * noise_ay),
	  ((pow(dt, 3) / 2.0) * noise_ax), 0, (pow(dt, 2) * noise_ax), 0,
	  0, ((pow(dt, 3) / 2.0) * noise_ay), 0, (pow(dt, 2) * noise_ay);

  ekf_.F_ = F;
  ekf_.Q_ = Q;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  z = measurement_pack.raw_measurements_;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  Hj = tools_.CalculateJacobian(ekf_.x_);
	  if (!Hj.isZero(0)) {
		  ekf_.H_ = Hj;
		  ekf_.R_ = R_radar;
		  ekf_.UpdateEKF(z);
	  }

  } else {
    // Laser updates
	  ekf_.H_ = H_laser;
	  ekf_.R_ = R_laser;
	  ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

  cout << "x_test = " << ekf_.x_test << endl;
  cout << "P_test = " << ekf_.P_test << endl;
}

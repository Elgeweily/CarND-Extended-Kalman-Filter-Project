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
  x_ = VectorXd(4);
  P_ = MatrixXd(4, 4);
  F_ = MatrixXd(4, 4);
  H_laser = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  R_laser = MatrixXd(2, 2);
  R_radar = MatrixXd(3, 3);
  Q_ = MatrixXd(4, 4);

  //object covariance matrix
  P_ << 1, 0, 0, 0,
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
        
      x_ << range * cos(bearing), range * sin(bearing), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */  
      x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }

	ekf_.Init(x_, P_, F_, H_laser, R_laser, Q_);

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
  F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  //Set the process covariance matrix Q
  Q_ << ((pow(dt, 4) / 4.0) * noise_ax), 0, ((pow(dt, 3) / 2.0) * noise_ax), 0,
	  0, ((pow(dt, 4) / 4.0) * noise_ay), 0, ((pow(dt, 3) / 2.0) * noise_ay),
	  ((pow(dt, 3) / 2.0) * noise_ax), 0, (pow(dt, 2) * noise_ax), 0,
	  0, ((pow(dt, 3) / 2.0) * noise_ay), 0, (pow(dt, 2) * noise_ay);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  z_ = measurement_pack.raw_measurements_;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  Hj_ = tools_.CalculateJacobian(x_);
	  ekf_.Init(x_, P_, F_, Hj_, R_radar, Q_);
	  ekf_.UpdateEKF(z_);
  } else {
    // Laser updates
	  ekf_.Init(x_, P_, F_, H_laser, R_laser, Q_);
	  ekf_.Update(z_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

  /**
  * RMSE and Jacobian Equations.
  */
  Tools tools_;

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;
  float dt;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::VectorXd x_;
  Eigen::MatrixXd P_;
  Eigen::MatrixXd F_;
  Eigen::MatrixXd H_laser;
  Eigen::MatrixXd Hj_;
  Eigen::MatrixXd R_laser;
  Eigen::MatrixXd R_radar;
  Eigen::MatrixXd Q_;

  Eigen::VectorXd z_;

  float noise_ax;
  float noise_ay;
};

#endif /* FusionEKF_H_ */

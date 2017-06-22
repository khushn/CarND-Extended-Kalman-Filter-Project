#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  //create a 4D state vector, we don't know yet the values of the x state
  Eigen::VectorXd x = VectorXd(4);

  //state covariance matrix P
  Eigen::MatrixXd P = MatrixXd(4, 4);
  P << 1, 0, 0, 0,
      0, 1, 0, 0, 
      0, 0, 1000, 0, 
      0, 0, 0, 1000;

  //the initial transition matrix F_
  Eigen::MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  // We let the process noise be a dummy initially
  // It is computed in the ProcessMesurement() function using
  noise_ax=9;
  noise_ay=9;
  Eigen::MatrixXd Q=MatrixXd(4, 4);

  ekf_=KalmanFilter();
  ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

/**
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // RADAR Not coded!!!

    return;
  }
  **/

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
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      cout << "rho: " << rho << ", phi: " << phi << ", rho dot: " << rho_dot << endl;
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
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
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  cout << "dt: " << dt << endl;
  if (dt <.05) {
    cout << "Skipping prediction step as dt (small): " << dt << endl;
    // TODO: YOUR CODE HERE
  } else {
  //1. Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2)=dt;
    ekf_.F_(1, 3)=dt;
    
    //2. Set the process covariance matrix Q
    float dt2=dt*dt;
    float dt3=dt*dt2;
    float dt4=dt*dt3;
    ekf_.Q_=MatrixXd(4, 4);
    ekf_.Q_ << (dt4/4)*noise_ax, 0, (dt3/2)*noise_ax, 0, 
                0, (dt4/4)*noise_ay, 0, (dt3/2)*noise_ay,
                (dt3/2)*noise_ax, 0, dt2*noise_ax, 0, 
                0, (dt3/2)*noise_ay, 0, dt2*noise_ay;

    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    float rho = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float rho_dot = measurement_pack.raw_measurements_[2];
    cout << "rho: " << rho << ", phi: " << phi << ", rho dot: " << rho_dot << endl;
    float px = rho * cos(phi);
    float py = rho * sin(phi);
    float vx = rho_dot * cos(phi);
    float vy = rho_dot * sin(phi);
    cout << "px: " << px << ", py: " << py << ", vx: " << vx << ", vy: " << vy << endl;
    VectorXd cartesian = VectorXd(4);
    cartesian << px, py, vx, vy;
    // Radar updates
    Tools tools;
    Hj_ = tools.CalculateJacobian(cartesian);
    VectorXd z = VectorXd(3);
    z << rho, phi, rho_dot;
    ekf_.UpdateEKF(z, Hj_, R_radar_);
  } else {
    // Laser updates
    VectorXd z = VectorXd(2);
    z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

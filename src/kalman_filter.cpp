#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const Eigen::MatrixXd & Hj, const Eigen::MatrixXd & R_rdr) {
  VectorXd y = z - compute_hx();
  // IMPORTANT: We need to normalize the phi differntial
  // Else, gives problems if ground truth is in 2nd quadrant (e.g. 170 degree)
  // And the prediction is in 3rd quadrand (e.g. -170 degrees)
  // So it calculates the difference as 340 degrees 
  y[1] = normalize_angle(y[1]);

  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_rdr;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}

VectorXd KalmanFilter::compute_hx() {
  VectorXd h_x = VectorXd(3);
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  h_x[0] = sqrt(px*px + py*py);

  float phi = 0.0;
  if (abs(px) < .001) {
    if (py>0) {
      phi = M_PI/2;
    } else {
      phi = - M_PI/2;
    }
  } else {
    // NOTE: Using atan does not distinguish between y -ve and x -ve, i.e. between 
    // 2nd and 3rd quadrant. So Atan2 must me used
    phi = atan2(py, px);
  }

  
  h_x[1] = phi;
  h_x[2] = (px*vx + py*vy) / h_x[0];
  cout << "inside compute_hx --" << endl;
  cout << "x_: " << x_ << endl;
  cout << "h_x: " << h_x[0] << ", "<< h_x[1] << ", " << h_x[2] << endl;
  return h_x;
}

float KalmanFilter::normalize_angle(float angle){
  //normalize pi to bring it between -pi and +pi
  while (angle > M_PI) {
    angle -= 2 * M_PI;
  }
  while (angle < -M_PI) {
    angle += 2 * M_PI;
  }
  return angle;
}
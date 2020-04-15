#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <iostream>
#include "kalman_filter.hpp"


KalmanFilter::KalmanFilter() {;}

KalmanFilter::~KalmanFilter() {;}

void KalmanFilter::Init(const Eigen::VectorXd &x_in,
                        const Eigen::VectorXd &u_in,
                        const Eigen::MatrixXd &P,
                        const Eigen::MatrixXd &F,
                        const Eigen::MatrixXd &B,
                        const Eigen::MatrixXd &Q,
                        const Eigen::MatrixXd &H,
                        const Eigen::MatrixXd &R,
                        const double t) {
  assert(x_in.size() == P.rows());
  assert(x_in.size() == P.cols());
  assert(x_in.size() == F.rows());
  assert(x_in.size() == F.cols());
  assert(x_in.size() == B.rows());
  assert(u_in.size() == B.cols());
  assert(x_in.size() == Q.rows());
  assert(x_in.size() == Q.cols());
  assert(x_in.size() == H.cols());
  x_ = x_in;
  u_ = u_in;
  P_ = P;
  F_ = F;
  B_ = B;
  Q_ = Q;
  H_ = H;
  R_ = R;
  t_ = t;
}


void KalmanFilter::Predict() {
  Eigen::MatrixXd Ft = F_.transpose();
  x_ = F_ * x_ + B_ * u_;
  P_ = F_ * P_ * Ft + Q_;
}


void KalmanFilter::Update(const Eigen::VectorXd &z) {
  assert(z.size() == H_.rows());
  assert(z.size() == R_.rows());
  assert(z.size() == R_.cols());
  Eigen::MatrixXd Ht = H_.transpose();
  Eigen::VectorXd z_pred = H_ * x_;
  Eigen::VectorXd y = z - z_pred;
  Eigen::MatrixXd S = H_ * P_ * Ht + R_;
  Eigen::MatrixXd Si = S.inverse();
  Eigen::MatrixXd K = P_ * Ht * Si;

  // new estimate
  x_ = x_ + (K * y);
  int x_size = x_.size();
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


void KalmanFilter::set_x(const Eigen::VectorXd &x_in) {
  assert(x_in.size() == x_.size());
  x_ = x_in;
}


void KalmanFilter::set_t(const double t_new) {
  t_ = t_new;
}

void KalmanFilter::set_F(const Eigen::MatrixXd &F) {
  assert(F.rows() == F_.rows());
  assert(F.cols() == F_.cols());
  F_ = F;
}

void KalmanFilter::set_H(const Eigen::MatrixXd &H) {
  assert(H.rows() == H_.rows());
  assert(H.cols() == H_.cols());
  H_ = H;
}

void KalmanFilter::set_R(const Eigen::MatrixXd &R) {
  assert(R.rows() == R_.rows());
  assert(R.cols() == R_.cols());
  R_ = R;
}

Eigen::VectorXd KalmanFilter::get_x() {
  return x_;
}

Eigen::MatrixXd KalmanFilter::get_P() {
  return P_;
}

double KalmanFilter::get_t() {
  return t_;
}

#ifndef INCLUDE_KALMAN_FILTER_HPP_
#define INCLUDE_KALMAN_FILTER_HPP_
#include "Eigen/Dense"

inline double DegToRad(const double deg) {
  return deg * M_PI / 180.;
}

inline double RadToDeg(const double rad) {
  return rad * 180. / M_PI;
}

class KalmanFilter {
 public:
  KalmanFilter();

  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * x_in Initial state
   * u_in Initial control-input
   * P Initial state covariance
   * F Transition matrix
   * B Control-input matrix
   * Q Process covariance matrix
   * H Measurement matrix
   * R Measurement covariance matrix
   */
  void Init(const Eigen::VectorXd &x_in,
            const Eigen::VectorXd &u_in,
            const Eigen::MatrixXd &P,
            const Eigen::MatrixXd &F,
            const Eigen::MatrixXd &B,
            const Eigen::MatrixXd &Q,
            const Eigen::MatrixXd &H,
            const Eigen::MatrixXd &R,
            const double t);

  void Predict();

  // z measurement
  void Update(const Eigen::VectorXd &z);

  void set_x(const Eigen::VectorXd &x_in);
  void set_F(const Eigen::MatrixXd &F);
  void set_H(const Eigen::MatrixXd &H);
  void set_Q(const Eigen::MatrixXd &Q);
  void set_R(const Eigen::MatrixXd &R);
  void set_t(const double t_new);

  Eigen::VectorXd get_x();
  Eigen::MatrixXd get_P();
  double get_t();

 private:
  // state vector
  Eigen::VectorXd x_;
  // control input vector
  Eigen::VectorXd u_;
  // state covariance matrix
  Eigen::MatrixXd P_;
  // state transistion matrix
  Eigen::MatrixXd F_;
  // control input matrix
  Eigen::MatrixXd B_;
  // process covariance matrix
  Eigen::MatrixXd Q_;
  // measurement matrix
  Eigen::MatrixXd H_;
  // measurement covariance matrix
  Eigen::MatrixXd R_;
  // estimate time stamp
  double t_;
};
#endif  // INCLUDE_KALMAN_FILTER_HPP_

#include <unistd.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <kalman_filter.hpp>

typedef Eigen::Vector2d Vector2;
typedef Eigen::Vector3d Vector3;
typedef Vector3 Point;

struct measurement {
  double time_;
  double range_;
  double azimuth_;
  double rcs_;
  double range_rate_;
  Eigen::Vector2d v_;
};

void kf_init(KalmanFilter& kf, measurement& meas);


int main() {
  std::vector<measurement> raw_data, kf_data;
  // Read data
  std::fstream file;
  std::ofstream outfile;

  char *path = NULL;
  size_t size;
  path = getcwd(path, size);
  std::string file_name, outfile_name;
  file_name = std::string(path) + "/data.csv";
  outfile_name = std::string(path) + "/data_out.csv";

  file.open(file_name);
  outfile.open(outfile_name);

  std::vector<double> pair_data;
  std::string line;
  std::string data;
  if (file.good()) {
    while (getline(file, line, '\n')) {
      pair_data.clear();
      std::istringstream templine(line);
      while (getline(templine, data, ',')) {
        pair_data.push_back(std::atof(data.c_str()));
      }
      measurement d{pair_data[0],
                    pair_data[1],
                    pair_data[2],
                    pair_data[3],
                    pair_data[4],
                    Vector2{0, 0}};
      raw_data.push_back(d);
    }
  } else {
    std::cout << "Cannot load file " << file_name << std::endl;
    return 0;
  }
  std::cout << raw_data.size() << std::endl;

  KalmanFilter kf;
  for (int i = 0; i < raw_data.size(); i++) {
    measurement meas = raw_data[i];
    if (kf_data.size() == 0) {
      kf_data.push_back(meas);
      kf_init(kf, meas);
      continue;
    }
    if ((raw_data[i].time_ - raw_data[i-1].time_) > 0.25) {
      kf_data.push_back(meas);
      kf_init(kf, meas);
      continue;
    }


    double dt = meas.time_ - kf.get_t();
    Eigen::MatrixXd F(4, 4);
    F << 1, 0, dt, 0,
         0, 1, 0, dt,
         0, 0, 1, 0,
         0, 0, 0, 1;

    kf.set_F(F);

    Eigen::MatrixXd G(4, 2), Qv(2, 2);
    G << dt * dt / 2, 0,
         0, dt * dt / 2,
         dt, 0,
         0, dt;
    Qv << 4, 0,
          0, 4;
    kf.set_Q(G * Qv * G.transpose());

    // Predict
    kf.Predict();
    kf.set_t(meas.time_);
    Eigen::VectorXd pred = kf.get_x();

    Eigen::MatrixXd H(3, 4);
    double H20 = pred[1] * (pred[2]*pred[1]-pred[3]*pred[0]) / std::pow(pred[0] * pred[0] + pred[1] * pred[1], 3/2.);
    double H21 = pred[0] * (pred[3]*pred[0]-pred[2]*pred[1]) / std::pow(pred[0] * pred[0] + pred[1] * pred[1], 3/2.);
    double H22 = pred[0] / std::sqrt(pred[0] * pred[0] + pred[1] * pred[1]);
    double H23 = pred[1] / std::sqrt(pred[0] * pred[0] + pred[1] * pred[1]);
    H << 1, 0, 0, 0,
         0, 1, 0, 0,
         H20, H21, H22, H23;
    kf.set_H(H);

    // Haven't find an approach to map uncertainty from polar to cartesian
    // Here is an approximation method.
    Eigen::MatrixXd J(2, 2), J_inv(2, 2), R0(2, 2), R(3, 3);
    double range = meas.range_;
    double azimuth = DegToRad(meas.azimuth_);
    double cov_range = (0.25/3)*(0.25/3);
    double cov_azimuth = (DegToRad(1.)/3)*(DegToRad(1.)/3);
    double cov_rate = (0.12/3)*(0.12/3);

    J << cos(azimuth), -range * sin(azimuth),
         sin(azimuth), range * cos(azimuth);
    J_inv << cos(azimuth), sin(azimuth),
             -range * sin(azimuth), range * cos(azimuth);
    R0 << cov_range, 0,
          0, cov_azimuth;  // Need to be set
    R0 = J * R0 * J_inv;
    R << R0(0, 0), R0(0, 1), 0,
         R0(1, 0), R0(1, 1), 0,
         0, 0, cov_rate;
    kf.set_R(R);

    // Update
    Eigen::VectorXd z_in(3);
    double x = meas.range_ * std::cos(DegToRad(meas.azimuth_));
    double y = meas.range_ * std::sin(DegToRad(meas.azimuth_));
    z_in << x, y, meas.range_rate_;


    kf.Update(z_in);
    Eigen::VectorXd update = kf.get_x();

    std::cout << std::to_string(meas.time_) << "|"
              << "(x, y, vx, vy) = raw => ("
              << x << ", "
              << y << ", none, none) v.s. "
              << "pred => ("
              << pred[0] << ", "
              << pred[1] << ", "
              << pred[2] << ", "
              << pred[3] << ") v.s. "
              << "update => ("
              << update[0] << ", "
              << update[1] << ", "
              << update[2] << ", "
              << update[3] << ")" << std::endl;

    if (std::isnan(pred[0]) || std::isnan(pred[1])
        || std::isnan(pred[2]) || std::isnan(pred[3])
        || (Point(pred[0], pred[1], 0) - Point(x, y, 0)).norm() > 3) {
      kf_data.push_back(meas);
      kf_init(kf, meas);
      std::cout << "Prediction ERROR!!!" << std::endl;
      continue;
    }
    double range_update = Point(update[0], update[1], 0).norm();
    double azimuth_update = RadToDeg(std::atan2(update[1], update[0]));
    double rate_update = (update[0] * update[2] + update[1] * update[3])
                         /std::sqrt(update[0] * update[0] + update[1] * update[1]);
    measurement meas_update{meas.time_,
                            range_update,
                            azimuth_update,
                            meas.rcs_,
                            rate_update,
                            Vector2{update[2], update[3]}};
    kf_data.push_back(meas_update);

    outfile << std::to_string(meas.time_) << ", "
            << x << ", "
            << y << ", "
            << pred[0] << ", "
            << pred[1] << ", "
            << update[0] << ", "
            << update[1] << std::endl;
  }
  file.close();
  outfile.close();
  return 0;
}

void kf_init(KalmanFilter& kf, measurement& meas) {
  double x = meas.range_ * std::cos(DegToRad(meas.azimuth_));
  double y = meas.range_ * std::sin(DegToRad(meas.azimuth_));
  // EKF initialization
  double range = meas.range_;
  double azimuth = DegToRad(meas.azimuth_);
  double cov_range = (0.25/3)*(0.25/3);
  double cov_azimuth = (DegToRad(1.)/3)*(DegToRad(1.)/3);
  double cov_rate = (0.12/3)*(0.12/3);

  Eigen::VectorXd x_in(4), u_in(4);
  Eigen::MatrixXd P(4, 4), F(4, 4), B(4, 4),
                  Q(4, 4), H(3, 4), R(3, 3);

  x_in << x, y, meas.v_[0], meas.v_[1];
  u_in << 0, 0, 0, 0;
  P << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1, 0,
       0, 0, 0, 1;
  F = B = P;

  Eigen::MatrixXd G(4, 2), Qv(2, 2);
  G << 0.125, 0,
       0, 0.125,
       0.5, 0,
       0, 0.5;
  Qv << 1, 0,
        0, 1;
  Q = G * Qv * G.transpose();

  double H20 = x_in[1] * (x_in[2]*x_in[1]-x_in[3]*x_in[0]) / std::pow(x_in[0] * x_in[0] + x_in[1] * x_in[1], 3/2.);
  double H21 = x_in[0] * (x_in[3]*x_in[0]-x_in[2]*x_in[1]) / std::pow(x_in[0] * x_in[0] + x_in[1] * x_in[1], 3/2.);
  double H22 = x_in[0] / std::sqrt(x_in[0] * x_in[0] + x_in[1] * x_in[1]);
  double H23 = x_in[1] / std::sqrt(x_in[0] * x_in[0] + x_in[1] * x_in[1]);
  H << 1, 0, 0, 0,
       0, 1, 0, 0,
       H20, H21, H22, H23;

  // Haven't find an approach to map uncertainty from polar to cartesian
  // Here is an approximation method.
  Eigen::MatrixXd J(2, 2), J_inv(2, 2), R0(2, 2);
  J << cos(azimuth), -range * sin(azimuth),
       sin(azimuth), range * cos(azimuth);
  J_inv << cos(azimuth), sin(azimuth),
           -range * sin(azimuth), range * cos(azimuth);
  R0 << cov_range, 0,
        0, cov_azimuth;  // Need to be set
  R0 = J * R0 * J_inv;
  R << R0(0, 0), R0(0, 1), 0,
       R0(1, 0), R0(1, 1), 0,
       0, 0, cov_rate;

  kf.Init(x_in, u_in, P, F, B, Q, H, R, meas.time_);
}

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
  double azimuth_rate_;
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
                    0};
      raw_data.push_back(d);
    }
  } else {
    std::cout << "Cannot load file " << file_name << std::endl;
    return 0;
  }
  std::cout << raw_data.size() << " radar measurements loaded." << std::endl;

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
    Qv << 1, 0,
          0, 4;
    kf.set_Q(G * Qv * G.transpose());

    // Predict
    kf.Predict();
    kf.set_t(meas.time_);
    Eigen::VectorXd pred = kf.get_x();
    // Update
    Eigen::VectorXd z_in(3);
    z_in << meas.range_, meas.azimuth_, meas.range_rate_;


    kf.Update(z_in);
    Eigen::VectorXd update = kf.get_x();

    std::cout << std::to_string(meas.time_) << "|"
              << "(range, azimuth, range_rate, azimuth_rate) = raw => ("
              << meas.range_ << ", "
              << meas.azimuth_ << ", "
              << meas.range_rate_ << ", none) v.s. "
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
        || std::fabs(pred[0] - meas.range_) > 1
        || std::fabs(pred[1] - meas.azimuth_) > 4
        || std::fabs(pred[2] - meas.range_rate_) > 0.48) {
      kf_data.push_back(meas);
      kf_init(kf, meas);
      std::cout << "Prediction ERROR!!!" << std::endl;
      continue;
    }

    measurement meas_update{meas.time_,
                            update[0],
                            update[1],
                            meas.rcs_,
                            update[2],
                            update[3]};
    kf_data.push_back(meas_update);

    double x = meas.range_ * std::cos(DegToRad(meas.azimuth_));
    double y = meas.range_ * std::sin(DegToRad(meas.azimuth_));

    double x_pred = pred[0] * std::cos(DegToRad(pred[1]));
    double y_pred = pred[0] * std::sin(DegToRad(pred[1]));

    double x_update = update[0] * std::cos(DegToRad(update[1]));
    double y_update = update[0] * std::sin(DegToRad(update[1]));

    outfile << std::to_string(meas.time_) << ", "
            << x << ", "
            << y << ", "
            << x_pred << ", "
            << y_pred << ", "
            << x_update << ", "
            << y_update << std::endl;
  }
  file.close();
  outfile.close();
  return 0;
}


void kf_init(KalmanFilter& kf, measurement& meas) {
  // EKF initialization
  double cov_range = (0.25/3)*(0.25/3);
  double cov_azimuth = (1./3)*(1./3);
  double cov_rate = (0.12/3)*(0.12/3);

  Eigen::VectorXd x_in(4), u_in(4);
  Eigen::MatrixXd P(4, 4), F(4, 4), B(4, 4),
                  Q(4, 4), H(3, 4), R(3, 3);

  x_in << meas.range_, meas.azimuth_, meas.range_rate_, 0;
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
        0, 4;
  Q = G * Qv * G.transpose();

  H << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1, 0;
  R << cov_range, 0, 0,
       0, cov_azimuth, 0,
       0, 0, cov_rate;

  kf.Init(x_in, u_in, P, F, B, Q, H, R, meas.time_);
}


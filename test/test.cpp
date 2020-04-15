#include <unistd.h>
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

inline double DegToRad(const double deg) {
  return deg * M_PI / 180.;
}

inline double RadToDeg(const double rad) {
  return rad * 180. / M_PI;
}

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


    Eigen::MatrixXd J(2, 2), J_inv(2, 2), R0(2, 2);
    double range = meas.range_;
    double azimuth = DegToRad(meas.azimuth_);

    J << cos(azimuth), -range * sin(azimuth),
         sin(azimuth), range * cos(azimuth);
    J_inv << cos(azimuth), sin(azimuth),
             -range * sin(azimuth), range * cos(azimuth);
    R0 << pow(1.0/3, 2), 0,
          0, pow(0.25/3, 2);  // Need to be set
    R0 = J * R0 * J_inv;
    kf.set_R(R0);

    double dt = meas.time_ - kf.get_t();
    Eigen::MatrixXd F(4, 4);
    F << 1, 0, dt, 0,
         0, 1, 0, dt,
         0, 0, 1, 0,
         0, 0, 0, 1;

    kf.set_F(F);

    // Predict
    kf.Predict();
    kf.set_t(meas.time_);
    Eigen::VectorXd pred = kf.get_x();
    // Update
    Eigen::VectorXd z_in(2);
    double x = meas.range_ * std::cos(DegToRad(meas.azimuth_));
    double y = meas.range_ * std::sin(DegToRad(meas.azimuth_));
    z_in << x, y;
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
        || (Point(pred[0], pred[1], 0) - Point(x, y, 0)).norm() > 5) {
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
  double cov_azimuth = (1./3)*(1./3);
  double cov_range = (0.25/3)*(0.25/3);

  Eigen::VectorXd x_in(4), u_in(4);
  Eigen::MatrixXd P(4, 4), F(4, 4), B(4, 4),
                  Q(4, 4), H(2, 4), R(2, 2);

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

  H << 1, 0, 0, 0,
       0, 1, 0, 0;

  Eigen::MatrixXd J(2, 2), J_inv(2, 2), R0(2, 2);
  J << cos(meas.azimuth_), -meas.range_ * sin(meas.azimuth_),
       sin(meas.azimuth_), meas.range_ * cos(meas.azimuth_);
  J_inv << cos(meas.azimuth_), sin(meas.azimuth_),
           -meas.range_ * sin(meas.azimuth_), meas.range_ * cos(meas.azimuth_);
  R0 << cov_azimuth, 0,
        0, cov_range;  // Need to be set
  R = J * R0 * J_inv;

  kf.Init(x_in, u_in, P, F, B, Q, H, R, meas.time_);
}

#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.95;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // State dimension
  n_x_ = 5;

  // Augumented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  //Prediction sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Time when state is true
  time_us_ = 0;

  is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Radar components
  double rho;
  double phi;
  double rho_dot;

  // Laser components
  double px;
  double py;

  // Verify if is_initialized is set to true and if timestamp is equal to zero - used specifically for the "Reset" action 
  if (!is_initialized_ || meas_package.timestamp_ == 1477010443000000) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Convert radar from polar to cartesian coordinates.
    */
    // First measurement
    x_ << 0.2, 0.4, 0, 0, 0;
      
    // Create and initialize the covarinace matrix P
    P_ << 0.1, 0, 0, 0, 0,
          0, 0.1, 0, 0, 0,
          0, 0, 0.1, 0, 0,
          0, 0, 0, 0.1, 0,
          0, 0, 0, 0, 0.1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      rho_dot = meas_package.raw_measurements_[2];
      
      double cos_phi = cos(phi);
      double sin_phi = sin(phi);
        
      double px = rho * cos_phi;
      double py = rho * sin_phi;
        
      x_ << px, py, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        
      // Set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
      
    time_us_ = meas_package.timestamp_;

    // Done initializing, no need to predict or update
    is_initialized_ = true;
      
    return;
  }

  // Compute the time elapsed between the current and previous measurements
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
  } else {
      UpdateLidar(meas_package);
  }

  // Print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Generate sigma points
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Q = MatrixXd(2, 2);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  Q(0,0) = std_a_ * std_a_;
  Q(0,1) = 0;
  Q(1,0) = 0;
  Q(1,1) = std_yawdd_ * std_yawdd_;
  P_aug.bottomRightCorner(2,2) = Q;

  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;

  for(int i=0; i<n_aug_; ++i) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  // Predict sigma points
  double v_k, psi_k, psi_dot_k, nu_a_k, nu_psi_ddot_k;

  for(int i=0; i<2 * n_aug_ + 1; ++i) {
    VectorXd temp_vector_x = VectorXd(5);
    VectorXd temp_vector_noise = VectorXd(5);

    v_k = Xsig_aug.col(i)(2);
    psi_k = Xsig_aug.col(i)(3);
    psi_dot_k = Xsig_aug.col(i)(4);
    nu_a_k = Xsig_aug.col(i)(5);
    nu_psi_ddot_k = Xsig_aug.col(i)(6);

    temp_vector_noise << 0.5 * delta_t * delta_t * cos(psi_k) * nu_a_k,
                         0.5 * delta_t * delta_t * sin(psi_k) * nu_a_k,
                         delta_t * nu_a_k,
                         0.5 * delta_t * delta_t * nu_psi_ddot_k,
                         delta_t * nu_psi_ddot_k;

    if(fabs(psi_dot_k) < 0.001) {
      temp_vector_x << v_k * cos(psi_k) * delta_t, v_k * sin(psi_k) * delta_t, 0, 0, 0;
    } else {
      temp_vector_x << (v_k / psi_dot_k) * (sin(psi_k + psi_dot_k * delta_t) - sin(psi_k)),
                       (v_k / psi_dot_k) * (-cos(psi_k + psi_dot_k * delta_t) + cos(psi_k)),
                       0,
                       psi_dot_k * delta_t,
                       0;
    }

    Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_) + temp_vector_x + temp_vector_noise;
  }

  // Predict mean/covariance matrix using sigma points
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for(int i=1; i<2 * n_aug_ + 1; i++) {
    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }

  // Predict state mean
  VectorXd x_mean = VectorXd(n_x_);
  x_mean.fill(0.0);

  for(int i=0; i<2 * n_aug_ + 1; i++) {
    x_mean = x_mean + weights_(i) * Xsig_pred_.col(i);
  }

  x_ = x_mean;

  // Predict state covariance matrix
  MatrixXd P_cal = MatrixXd(n_x_, n_x_);
  P_cal.fill(0.0);

  for(int i=0; i<2 * n_aug_ + 1; i++) {
    // State difference
    MatrixXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);

    // Angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_cal = P_cal + weights_(i) * x_diff * x_diff.transpose();
  }

  P_ = P_cal;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Get Lidar components
  double px_z = meas_package.raw_measurements_[0];
  double py_z = meas_package.raw_measurements_[1];

  VectorXd z = VectorXd(2);

  z << px_z, py_z;

  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  double px, py, v, psi, psi_dot;

  for(int i=0; i<2 * n_aug_ + 1; ++i) {
    px = Xsig_pred_.col(i)(0);
    py = Xsig_pred_.col(i)(1);
    v = Xsig_pred_.col(i)(2);
    psi = Xsig_pred_.col(i)(3);
    psi_dot = Xsig_pred_.col(i)(4);

    VectorXd temp_vector = VectorXd(n_z);
    temp_vector(0) = px;
    temp_vector(1) = py;

    Zsig.col(i) = temp_vector;
  }

  Estimate(meas_package, z, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Get Radar components
  double rho = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double rho_dot = meas_package.raw_measurements_[2];

  VectorXd z = VectorXd(3);
  z << rho, phi, rho_dot;

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  double px, py, v, psi, psi_dot;

  for(int i=0; i<2*n_aug_+1; ++i) {
    px = Xsig_pred_.col(i)(0);
    py = Xsig_pred_.col(i)(1);
    v = Xsig_pred_.col(i)(2);
    psi = Xsig_pred_.col(i)(3);
    psi_dot = Xsig_pred_.col(i)(4);

    VectorXd temp_vector = VectorXd(n_z);

    temp_vector(0) = sqrt(px*px + py*py);
    temp_vector(1) = atan2(py, px);
    temp_vector(2) = (px * cos(psi) * v + py * sin(psi) * v) / sqrt(px*px + py*py);

    // Angle normalization
    while(temp_vector(1) > M_PI) temp_vector(1) -= 2.*M_PI;
    while(temp_vector(1) < -M_PI) temp_vector(1) += 2.*M_PI;

    Zsig.col(i) = temp_vector;
  }

  Estimate(meas_package, z, Zsig, n_z);  
}

/**
 * Calculates the estimations for either Radar or Lidar
 * @param meas_package The measurement at k+1
 * @param Zsig The matrix with predictions
 * @param n_z The number of components for Lidar and Radar readings
 */
void UKF::Estimate(MeasurementPackage meas_package, VectorXd z, const MatrixXd Zsig, int n_z) {
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // Calculate mean predicted measurement
  for(int i=0; i<2 * n_aug_ + 1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Calculate measurement covariance matrix S
  MatrixXd Component_R = MatrixXd(n_z, n_z);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
     Component_R << std_radr_ * std_radr_, 0, 0,
                    0, std_radphi_ * std_radphi_, 0,
                    0, 0, std_radrd_ * std_radrd_;
  } else {
      Component_R << std_laspx_ * std_laspx_, 0,
                     0, std_laspy_ * std_laspy_;
  }

  for(int i=0; i<2*n_aug_+1; ++i) {
    S = S + weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
  }

  S = S + Component_R;

  MatrixXd z_diff = z - z_pred;

  MatrixXd temp_vector2 = z_diff.transpose() * S.inverse();

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    double nis_radar_i = (temp_vector2 * z_diff).value();
    nis_radar_.push_back(nis_radar_i);
  } else {
    double nis_lidar_i = (temp_vector2 * z_diff).value();
    nis_lidar_.push_back(nis_lidar_i);
  }
  
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // Calculate cross correlation matrix
  for(int i=0; i<2*n_aug_+1;++i) {
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);

    Tc = Tc + weights_(i) * x_diff * (Zsig.col(i) - z_pred).transpose();
  }

  // Calculate Kalman Gain K
  MatrixXd K = Tc * S.inverse();

  // Update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
}

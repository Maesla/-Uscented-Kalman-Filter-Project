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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* State dimension
   n_x_ = 5;

   ///* Augmented state dimension
   n_aug_ = 7;

   ///* Sigma point spreading parameter
   lambda_ = 3 - n_aug_;

   ///* Weights of sigma points
    weights_ = VectorXd(2*n_aug_+1);

    P_ = MatrixXd::Identity(n_x_, n_x_);

    Q_ = MatrixXd(2,2);
    Q_(0,0) = std_a_*std_a_;
    Q_(1,1) = std_yawdd_*std_yawdd_;

    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ +1);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_)
  {
    time_us_ = meas_package.timestamp_;
    x_.fill(0);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      cout << "Init by Radar" << endl;

      float phi = meas_package.raw_measurements_[0];
      float rho = meas_package.raw_measurements_[1];
      float px = cos(rho)*phi;
      float py = sin(rho)*phi;
      x_(0) = px;
      x_(1) = py;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      cout << "Init by Laser" << endl;

      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    is_initialized_ = true;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  cout << "Count: " << count << endl;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else //LASER
  {

  }

  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_<< endl;


  count ++;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd augmentedSigmaPoints = GenerateAugmentedSigmaPoints();

  cout << "Augmented Sigma Points = " << endl << augmentedSigmaPoints << endl;

  PredictAugmentedSigmaPoints(augmentedSigmaPoints, delta_t);

  cout << "Predicted Sigma Points = " << endl << Xsig_pred_ << endl;

  PredictMeanAndCovariance();

  //Xsig_pred_ = predictedSigmaPoints;
}

MatrixXd UKF::GenerateAugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  int num_col = 2 * n_aug_ + 1;
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, num_col);

  x_aug.fill(0);
  x_aug.head(n_x_) = x_;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q_;

  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug.col(0)  = x_aug;
  MatrixXd offset = sqrt(lambda_+n_aug_)*A;

  for (int i = 0; i < n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + offset.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - offset.col(i);
  }

  return Xsig_aug;
}


void UKF::PredictAugmentedSigmaPoints(MatrixXd augmentedSigmaPoints, float dt)
{
  int num_col = 2 * n_aug_ + 1;

  //MatrixXd Xsig_pred = MatrixXd(n_x_, num_col);
  Xsig_pred_.fill(0.0);

  for(int i = 0; i < num_col; i++)
  {
      Xsig_pred_.col(i) = PredictAugmentedSigmaPoint(augmentedSigmaPoints.col(i), dt);
  }

  //return Xsig_pred;
}

VectorXd UKF::PredictAugmentedSigmaPoint(VectorXd xsig_pred_col, float dt)
{
    double px = xsig_pred_col[0];
    double py = xsig_pred_col[1];
    double v = xsig_pred_col[2];
    double yaw = xsig_pred_col[3];
    double yaw_rate = xsig_pred_col[4];
    double acceleration_noise = xsig_pred_col[5];
    double yaw_acceleration_noise = xsig_pred_col[6];

    bool isYawRateZero = fabs(yaw_rate) < 0.01;

    double px_predicted = 0;
    double py_predicted = 0;

    if(isYawRateZero)
    {
        px_predicted = v*cos(yaw)*dt;
        py_predicted = v*sin(yaw)*dt;
    }
    else
    {
        px_predicted = (v/yaw_rate)*(sin(yaw+yaw_rate*dt)-sin(yaw));
        py_predicted = (v/yaw_rate)*(-cos(yaw+yaw_rate*dt)+cos(yaw));
    }

    double v_predicted = 0;
    double yaw_predicted = yaw_rate*dt;
    double yaw_rate_predicted = 0;

    double dt_squared = dt*dt;
    double px_noise = 0.5*dt_squared*cos(yaw)*acceleration_noise;
    double py_noise = 0.5*dt_squared*sin(yaw)*acceleration_noise;
    double v_noise = dt*acceleration_noise;
    double yaw_noise = 0.5*dt_squared*yaw_acceleration_noise;
    double yaw_rate_noise = dt*yaw_acceleration_noise;

    VectorXd x_predicted = VectorXd(5);
    x_predicted(0) = px + px_predicted + px_noise;
    x_predicted(1) = py + py_predicted + py_noise;
    x_predicted(2) = v + v_predicted + v_noise;
    x_predicted(3) = yaw + yaw_predicted + yaw_noise;
    x_predicted(4) = yaw_rate + yaw_rate_predicted + yaw_rate_noise;


    return x_predicted;
}

void UKF::PredictMeanAndCovariance(){

    weights_.fill(1.0/(2.0*(lambda_+n_aug_)));
    weights_[0] = lambda_/(lambda_+n_aug_);

    x_.fill(0.0);
    //predict state mean
    int num_col = 2*n_aug_+1;
    for (int i = 0; i < num_col; i++)
    {
        x_ = x_ + weights_[i]*Xsig_pred_.col(i);
    }


    P_.fill(0.0);
    //predict state covariance matrix
    for (int i = 0; i < num_col; i++)
    {
      VectorXd aux = SubstractAndKeepAngleNormalized(Xsig_pred_.col(i), x_, 3);
      P_ = P_ + weights_[i]*aux*aux.transpose();
    }

}

VectorXd UKF::SubstractAndKeepAngleNormalized(VectorXd a, VectorXd b, int indexOfAngle){
  VectorXd aux = a - b;
  while (aux(indexOfAngle)> M_PI) aux(indexOfAngle)-=2.*M_PI;
  while (aux(indexOfAngle)<-M_PI) aux(indexOfAngle)+=2.*M_PI;

  return aux;
}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3;
  int num_col = 2 * n_aug_ + 1;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, num_col);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  weights_.fill(1.0/(2.0*(lambda_+n_aug_)));
  weights_[0] = lambda_/(lambda_+n_aug_);


  if (count > 60)
  {
    int a = 1;
  }
  for (int i = 0; i < num_col; i++)
  {
    double px = Xsig_pred_.col(i)[0];
    double py = Xsig_pred_.col(i)[1];
    double v = Xsig_pred_.col(i)[2];
    double yaw = Xsig_pred_.col(i)[3];

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    cout << "Px Py" << endl << px << " " << py << endl;

    double rho = sqrt(px*px + py*py);
    double phi = atan2(py,px);

    cout << "Rho " << endl << rho << endl;

    //double rho_rate = (px*v1 + py*v2 ) / sqrt(px*px + py*py);
    double rho_rate = 0.0;
    if (fabs(rho) > 0.1)
      rho_rate = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;

    //Zsig.col(i)[0] = rho;
    //Zsig.col(i)[1] = phi;
    //Zsig.col(i)[2] = rho_rate;

    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_rate;
  }

  cout << "Z sigma = " << endl << Zsig << endl;


  z_pred.fill(0.0);
  for(int i = 0; i < num_col; i++)
  {
    z_pred = z_pred + weights_[i]*Zsig.col(i);
  }

  S.fill(0.0);
  for(int i = 0; i < num_col; i++)
  {
    VectorXd aux = SubstractAndKeepAngleNormalized(Zsig.col(i), z_pred, 1);
    S = S + weights_[i]*aux*aux.transpose();
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  S = S + R;

  //Update Asigment

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < num_col; i++)
  {

    try {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;


      cout << "X diff = " << endl << x_diff << endl;
      cout << "Z diff = " << endl << z_diff << endl;


      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    catch (...)
    {
      cout << "Exception occurred";
    }

//    Tc = Tc + weights_[i] * x_diff * z_diff.transpose();
  }


  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //VectorXd z = meas_package.raw_measurements_;
  VectorXd z = VectorXd(n_z);
  z.fill(0.0);
  z(0) = meas_package.raw_measurements_[0];
  z(1) = meas_package.raw_measurements_[1];
  z(2) = meas_package.raw_measurements_[2];


  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();


}

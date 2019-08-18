#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  x_ << 1,1,1,1,1;
  P_ << 1, 0, 0, 0, 0,
	    0, 1, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1);
  //VectorXd weights = VectorXd(2*n_aug_+1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
	
	if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
	   float ro     = meas_package.raw_measurements_(0);
	  float phi    = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
	  x_(0) = ro     * cos(phi);
	  x_(1) = ro     * sin(phi);	
}
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		float px = meas_package.raw_measurements_(0);
		float py = meas_package.raw_measurements_(1);
		x_(0) = px;
		x_(1) = py;
	}
	//previous_timestamp_ =meas_package.timestamp_;
    is_initialized_ = true;
	return;
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  int n_sigma = 2 * n_x_ + 1;
  // augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
  //augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(7, 7);
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  //created augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_;i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*L.col(i);
  }
  // Predicted sigma points
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1);
  //predict sigma points
  for(int i=0;i<2*n_aug_+1;i++)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2,i);
	double yaw = Xsig_aug(3,i);
	double yawd = Xsig_aug(4,i);
	double nu_a = Xsig_aug(5,i);
	double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p,py_p;
	
	//avoiding division by zero
    if (fabs(yawd) > 0.001){
		px_p = p_x + v/yawd*(sin(yaw + yawd*delta_t) - sin(yaw));
		py_p = p_y + v/yawd*(cos(yaw) - cos(yaw+yawd*delta_t) );
	}
	else {
		px_p = p_x + v*delta_t*cos(yaw);
		py_p = p_y + v*delta_t*sin(yaw);
	}

	double v_p = v;
	double yaw_p = yaw + yawd*delta_t;
	double yawd_p = yawd;

	//add noise
	px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
	v_p = v_p + nu_a*delta_t;

	yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
	yawd_p = yawd_p + nu_yawdd*delta_t;
   
   //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //weights
  //VectorXd weights = VectorXd(2*n_aug_+1);
  // predicted state
  weights_(0) = lambda_/(lambda_+ n_aug_);
  for(int i=1;i<2*n_aug_+1;i++){
	double weight = 0.5/(lambda_ + n_aug_);
	weights_(i) = weight;
  } 
  x_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
  	x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }
  P_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
	VectorXd x_diff = Xsig_pred_.col(i) - x_;
    std::cout<<"Test "<<M_PI;
	while(x_diff(3)>M_PI){
		x_diff(3) = x_diff(3) - 2.*M_PI;
	}
	while(x_diff(3)<-M_PI){
		x_diff(3) += 2.*M_PI;
	}
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd H = MatrixXd(2,4);
  MatrixXd R = MatrixXd(2,2);
  H << 1, 0, 0, 0,
       0, 1, 0, 0;
  R << 0.0225, 0,
  	0,	0.0225;
  VectorXd z_pred = H*x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H*P_*Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_*Ht;
  MatrixXd K = PHt*Si;

  //new estimate
  x_ = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H) + P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z=3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for(int i=0; i<2*n_aug_+1;i++){
	double p_x = Xsig_pred_(0,i);
	double p_y = Xsig_pred_(1,i);
	double v = Xsig_pred_(2,i);
	double yaw = Xsig_pred_(3,i);
	
	double v1 = cos(yaw)*v;
	double v2 = sin(yaw)*v;

	//measurement model
	Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
	Zsig(0,i) = atan2(p_y,p_x);
	Zsig(2,i) = (p_x*v1 + p_y*v2)/Zsig(0,i);
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
  	z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
	VectorXd z_diff = Zsig.col(i) - z_pred;

	while(z_diff(1)>M_PI){
	  z_diff(1) -= 2.*M_PI;
	}
	while(z_diff(1)<M_PI){
	  z_diff(1) += 2.*M_PI;
	}

	S += weights_(i)*z_diff*z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R<<std_radr_*std_radr_,0,0,
	 0, std_radphi_*std_radphi_, 0,
	 0, 0, std_radrd_*std_radrd_;
  
  S = S + R;

  // modify it (creating example vetor for incoming radar measurement)
  VectorXd z = VectorXd(n_z);
  z<<5.9214,
	 0.2187,
	 2.0062;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){

  	VectorXd z_diff = Zsig.col(i) - z_pred;
	while(z_diff(1)>M_PI) z_diff(1) -= 2.*M_PI;
	while(z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	VectorXd x_diff = Xsig_pred_.col(i) - x_;
	while(x_diff(3)>M_PI) x_diff(3) -= 2.*M_PI;
	while(x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

	Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }
  //kalman gain k
  MatrixXd K = Tc*S.inverse();
  VectorXd z_diff = z - z_pred;

  while(z_diff(1)>M_PI)z_diff(1) -= 2.*M_PI;
  while(z_diff(1)<-M_PI)z_diff(1) += 2.*M_PI;

  // update state mean and covariance
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose(); 
}

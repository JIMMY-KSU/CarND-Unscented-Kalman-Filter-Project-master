#include "ukf.h"
#include "tools.h"
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
	n_x_ = 5;
	n_aug_ = 7;
	n_sig_ = 2*n_aug_+1;
	lambda_ = 3 - n_aug_;


  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_ << 1, 1, 0, 0, 0;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5; //0.172273; //1.5; // default was 30;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//0.527692; // 0.5; // default was 30;


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
  Xsig_pred_= MatrixXd(n_x_, 2 * n_aug_ + 1);

  time_us_ = 0.0;

  NIS_radar_ =0.0; // NIS is normalized innovation squared.

  NIS_laser_ =0.0;

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_ / (lambda_+n_aug_);
    for (int i =1;i<2*n_aug_+1;i++) {
  	  weights_(i) = 0.5/(lambda_+n_aug_);
    }

  is_initialized_ = false;


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
	if (!is_initialized_) {
		time_us_ = meas_package.timestamp_;
		if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			cout << "Laser Measurement" << endl;
			float px = meas_package.raw_measurements_[0];
			float py = meas_package.raw_measurements_[1];

			// CTRV Model State Vector
			x_ << px, py, 0, 0, 0;

		}
		else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
			cout << "Radar Measurement" << endl;
			float rho = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];
			float rho_dot = meas_package.raw_measurements_[2];

			float px = rho * cos(phi);
			float py = rho * sin(phi);

			// CTRV Model State Vector
			x_ << px, py, 0, 0, 0;
		}
		is_initialized_ = true;
	}
	else {
		  /*****************************************************************************
		   *  Prediction
		   ****************************************************************************/
		double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
		cout << "delta_t = " << delta_t << endl;
		time_us_ = meas_package.timestamp_;
		cout << "Start Prediction" << endl;
		Prediction(delta_t);
		cout << "Done Prediction" << endl;
		//Switch between lidar and radar measurements
		if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			cout << "Start UpdateLidar" << endl;
			cout << "meas px = " << meas_package.raw_measurements_[0] << "\t meas py = " << meas_package.raw_measurements_[1] << endl;
			UpdateLidar(meas_package);
			cout << "Done UpdateLidar" << endl;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			cout << "Start UpdateRadar" << endl;
			cout << "meas rho = " << meas_package.raw_measurements_[0] << "\t meas phi = " << meas_package.raw_measurements_[1] << "\t meas rhodot=" << meas_package.raw_measurements_[2]  <<  endl;
			UpdateRadar(meas_package);
			cout << "Done UpdateRadar" << endl;
		}
	}


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
/** 1. AugmentedSigmaPoints **/
	cout << "/** 1. AugmentedSigmaPoints **/" << endl;
	//create augmented mean state
	VectorXd x_aug_(n_aug_);//n_aug_ = 7
	x_aug_.fill(0);
	x_aug_.head(n_x_) = x_;

	//create augmented covariance matrix
	MatrixXd P_aug_(n_aug_,n_aug_);
	P_aug_.fill(0);
	P_aug_.topLeftCorner(n_x_,n_x_) = P_;
	P_aug_(5,5) = std_a_*std_a_; // pow(std_a_,2)
	P_aug_(6,6) = std_yawdd_*std_yawdd_; // pow(std_yawdd_,2)

	//create square root matrix
	MatrixXd L_ = P_aug_.llt().matrixL();

	//create augmented sigma points
	MatrixXd Xsig_aug_(n_aug_, 2*n_aug_+1);
	Xsig_aug_.col(0) = x_aug_;
	for (int i=0;i<n_aug_;i++) {
		Xsig_aug_.col(i+1) 				= x_aug_ + sqrt(lambda_ + n_aug_)* L_.col(i);
		Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L_.col(i);
	}

	 cout << "Xsig_aug_ \n"<< Xsig_aug_ <<endl;

/** 2. SigmaPointPrediction **/
	cout << "/** 2. SigmaPointPrediction **/" << endl;
	for (int i=0;i<2*n_aug_+1;i++) {
		//extract values for better readability
		double p_x = Xsig_aug_(0,i);
		double p_y = Xsig_aug_(1,i);
		double v = Xsig_aug_(2,i);
		double yaw = Xsig_aug_(3,i);
		double yawd = Xsig_aug_(4,i);
		double nu_a = Xsig_aug_(5,i);
		double nu_yawdd = Xsig_aug_(6,i);

		//predicted state values
		double px_p, py_p;
		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}
		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;
		//yaw_p = yaw_p - (2. * M_PI) * floor((yaw_p + M_PI)/(2. *M_PI));

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;


	}
	cout << "Xsig_pred_ \n" << Xsig_pred_ <<endl;

/** 3. PredictMeanAndCovariance **/
	cout << "/** 3. PredictMeanAndCovariance **/" << endl;
	//predict state mean
	/*for (int i = 0; i<2*n_aug_+1; i++)
	  {
	      x_ = x_+weights_(i)*Xsig_pred_.col(i);
	  }*/ //This is where I had issue from for two weeks!
	  x_ = Xsig_pred_ * weights_;

	//predict state covariance matrix
	P_.fill(0.0);
	for (int i =0;i<2*n_aug_+1;i++) {
		//cout << "i=" << i << endl;
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//cout << "before x_diff(3)=" << x_diff(3) << endl;
		//while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		//while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;//this will take too long time hence propose to use ceil function
		x_diff(3) = x_diff(3) - (2. * M_PI) * floor((x_diff(3) + M_PI)/(2. *M_PI));
		//x_diff(3) = x_diff(3)-ceil((x_diff(3)-M_PI)/(2.*M_PI))*2.*M_PI;
		//cout << "after x_diff(3)=" << x_diff(3) << endl;
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

	}
	cout << "P_  \n" << P_  <<endl;



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
	//transform sigma points into measurement space
	cout << "Started Laser Update!" << endl;
	if (x_(0) == 0) return;
	MatrixXd Zsig_ = MatrixXd(2, 2*n_aug_+1);
	VectorXd z_pred_ = VectorXd(2);

	//measurement model
	Zsig_.row(0) = Xsig_pred_.row(0); //px
	Zsig_.row(1) = Xsig_pred_.row(1); //py


	//z_pred_ = Zsig_ * weights_;
	z_pred_.fill(0.0);
	  for (int i=0; i < 2*n_aug_+1; i++) {
	      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
	  }
	//create measurement covariance S
	MatrixXd S_(2,2);
	S_.fill(0);

	//create matrix for cross correlation Tc
	MatrixXd Tc_(n_x_,2);
	Tc_.fill(0);

	for (int i=0;i<2*n_aug_+1 ;i++){
		//measurement residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;
		//cout << "before z_diff(1)=" << z_diff(1) << endl;
		//z_diff(1) = z_diff(1) - (2. * M_PI) * floor((z_diff(1) + M_PI)/(2. *M_PI));
		//angle normalization
		//		while (z_diff(1)>M_PI) z_diff(1) -= 2*M_PI;
		//		while (z_diff(1)<-M_PI) z_diff(1) += 2*M_PI;
		//cout << "after z_diff(1)=" << z_diff(1) << endl;
		S_ = S_ + weights_(i) * z_diff * z_diff.transpose();

		//State residual
		VectorXd tx_diff = Xsig_pred_.col(i) - x_;
		cout << "before tx_diff(3)=" << tx_diff(3) << endl;
		tx_diff(3) = tx_diff(3)  - (2. * M_PI) * floor((tx_diff(3) + M_PI)/(2. *M_PI));
		//angle normalization
		//		while (tx_diff(3)>M_PI) tx_diff(3)-= 2.*M_PI;
		//		while (tx_diff(3)<-M_PI) tx_diff(3)+= 2.*M_PI;
		cout << "after tx_diff(3)=" << tx_diff(3) << endl;
		Tc_ = Tc_ + weights_(i) * tx_diff * z_diff.transpose();

	}

	//Laser measurement noise covariance matrix
	MatrixXd R_laser_(2,2);
	R_laser_ << std_laspx_*std_laspx_,0,
						0,std_laspy_*std_laspy_;

	S_ = S_ + R_laser_ ;
	cout << "S_  \n" << S_<< endl;
	//calculate unscented kalman filter gain K_
	cout << "Tc_  \n" << Tc_<< endl;
	cout << "S_.inverse() \n" << S_.inverse()<< endl;
	MatrixXd K_ = Tc_ * S_.inverse();
	cout << "K_ \n" << K_<< endl;
	//Actual measurement - predicted measurement
	VectorXd A_ = meas_package.raw_measurements_ - z_pred_;
	//cout << "before A_(1)=" << A_(1) << endl;
	//A_(1) = A_(1) - (2. * M_PI) * floor((A_(1) + M_PI)/(2. *M_PI));
	//angle normalization
	//			while (A_(1)>M_PI) A_(1) -= 2*M_PI;
	//			while (A_(1)<-M_PI) A_(1) += 2*M_PI;
	//cout << "after A_(1) =" << A_(1)<< endl;
	//cout << "A_ \n" << A_<< endl;
	//Update state x
	x_ = x_ + K_*A_;
	cout << "x_ \n" << x_<< endl;
	//Update state covariance
	P_ = P_ - K_ * S_ * K_.transpose();

	//Update NIS
	NIS_laser_ = A_.transpose() * S_.inverse() * A_;



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
	cout << "Started Radar Update!" << endl;
	if (x_(0) == 0) return;
	//transform sigma points into measurement space
	MatrixXd Zsig_ = MatrixXd(3, 2*n_aug_+1);
	VectorXd z_pred_ = VectorXd(3);

	//measurement model
	//transform sigma points into measurement space
	  for (int i = 0;i<2*n_aug_+1;i++) {
		  double p_x = Xsig_pred_(0,i);
		  double p_y = Xsig_pred_(1,i);
		  double v = Xsig_pred_(2,i);
		  double yaw = Xsig_pred_(3,i);
		  double yawd = Xsig_pred_(4,i);

		  double v1 = cos(yaw)*v;
		  double v2 = sin(yaw)*v;

		  Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y); //rho
		  Zsig_(1,i) = atan2(p_y,p_x);			//phi
		  Zsig_(1,i) = Zsig_(1,i)  - (2. * M_PI) * floor((Zsig_(1,i) + M_PI)/(2. *M_PI));
		  if (Zsig_(0,i) > 0.001) {
			  Zsig_(2,i) = (p_x*v1+p_y*v2)/Zsig_(0,i); //rho_dot
		  }
		  else {
			  Zsig_(2,i) = 0; //rho_dot
		  }

	  }

	z_pred_ = Zsig_ * weights_;
	z_pred_(1) = z_pred_(1)  - (2. * M_PI) * floor((z_pred_(1) + M_PI)/(2. *M_PI));
	//create measurement covariance S
	MatrixXd S_(3,3);
	S_.fill(0);


	//create matrix for cross correlation Tc
	MatrixXd Tc_(n_x_,3);
	Tc_.fill(0);

	for (int i=0;i<2*n_aug_+1 ;i++){
		//measurement residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;
		cout << "before z_diff(1)=" << z_diff(1) << endl;
		//angle normalization
		//while (z_diff(1)>M_PI) z_diff(1) -= 2*M_PI;
		//while (z_diff(1)<-M_PI) z_diff(1) += 2*M_PI;
		z_diff(1) = z_diff(1)  - (2. * M_PI) * floor((z_diff(1) + M_PI)/(2. *M_PI));
		cout << "after z_diff(1)=" << z_diff(1) << endl;

		S_ = S_ + weights_(i) * z_diff * z_diff.transpose();

		//State residual
		VectorXd tx_diff = Xsig_pred_.col(i) - x_;
		cout << "before tx_diff(3)=" << tx_diff(3) << endl;
		//angle normalization
		//while (tx_diff(3)>M_PI) tx_diff(3)-= 2.*M_PI;
		//while (tx_diff(3)<-M_PI) tx_diff(3)+= 2.*M_PI;
		tx_diff(3) = tx_diff(3)  - (2. * M_PI) * floor((tx_diff(3) + M_PI)/(2. *M_PI));
		cout << "after tx_diff(3)=" << tx_diff(3) << endl;
		//cout<< "tx_diff \n" <<tx_diff << endl;
		Tc_ = Tc_ + weights_(i) * tx_diff * z_diff.transpose();

	}

	//Laser measurement noise covariance matrix
	MatrixXd R_radar_(3,3);
	R_radar_ << std_radr_*std_radr_,0,0,
						0,std_radphi_*std_radphi_,0,
						0,0,std_radrd_*std_radrd_;

	S_ = S_ + R_radar_ ;
	cout << "S_  \n" << S_<< endl;

	//calculate unscented kalman filter gain K_
	cout << "Tc_  \n" << Tc_<< endl;
		cout << "S_.inverse() \n" << S_.inverse()<< endl;
	MatrixXd K_ = Tc_ * S_.inverse();
	cout << "K_ \n" << K_<< endl;
	//Actual measurement - predicted measurement
	VectorXd A_ = meas_package.raw_measurements_ - z_pred_;
	cout << "before A_(1)=" << A_(1) << endl;
	A_(1) = A_(1) - (2. * M_PI) * floor((A_(1) + M_PI)/(2. *M_PI));
	//angle normalization
	//		while (A_(1)>M_PI) A_(1) -= 2*M_PI;
	//		while (A_(1)<-M_PI) A_(1) += 2*M_PI;
	cout << "after A_(1) =" << A_(1)<< endl;
	//Update state x
	x_ = x_ + K_*A_;
	cout << "x_" << x_<< endl;
	//Update state covariance
	P_ = P_ - K_ * S_ * K_.transpose();

	//Update NIS
	NIS_radar_ = A_.transpose() * S_.inverse() * A_;




}

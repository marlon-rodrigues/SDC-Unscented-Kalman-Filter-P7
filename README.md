# SDC-Unscented-Kalman-Filter-P7

## Unscented Kalman Filter Project
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project required obtaining RMSE values that are lower that the tolerance outlined in the [project reburic](https://review.udacity.com/#!/rubrics/783/view). 

Data File
---
The data file information is provided by the simulator.

Each line in the data file represents a sensor measurement where the first column tells you if the measurement comes from radar (R) or lidar (L).

For a row containing radar data, the columns are: sensor_type, rho_measured, phi_measured, rhodot_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth.

For a row containing lidar data, the columns are: sensor_type, x_measured, y_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth.

Whereas radar has three measurements (rho, phi, rhodot), lidar has two measurements (x, y).


File Structure
---
The project files are located in the "src" folder. The files are the following:

*main.cpp - communicates with the Term 2 Simulator receiving data measurements, calls a function to run the Unscented Kalman filter, calls a function to calculate RMSE. Reads in the data and sends a sensor measurement to ufk.cpp.

*ukf.cpp - initializes the filter, defines the predict function, the update function for lidar, and the update function for radar.

*tools.cpp- function to calculate RMSE.

Video with Results
---
A video of the results named "Unscented-KF-P7.mov" is located in the repository. Lidar measurements are red circles, radar measurements are blue circles with an arrow pointing in the direction of the observed angle, and estimation markers are green triangles. The video shows what the simulator looks like when a c++ script is using its Unscented Kalman filter to track the object. The simulator provides the script the measured data (either lidar or radar), and the script feeds back the measured estimation marker, and RMSE values from its Kalman filter.

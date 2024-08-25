# Unscented Kalman Filter (UKF) Implementation

<img src="media/ukf_highway_tracked.gif" width="700" height="400" />


This project implements an Unscented Kalman Filter (UKF) in C++ for sensor fusion, utilizing data from noisy radar and lidar measurements to estimate the state of multiple cars on a highway. The UKF is designed to estimate the state of a moving object (e.g., a car) by processing noisy sensor data and predicting the object's position, velocity, and orientation.

## Overview

The Unscented Kalman Filter is a powerful algorithm that improves on the limitations of the Extended Kalman Filter (EKF) by providing a better approximation of nonlinear functions. The UKF leverages a set of sigma points to capture the mean and covariance more accurately, making it suitable for complex scenarios involving highly nonlinear dynamics and measurement models.

This implementation is structured to handle real-time data from laser and radar sensors, allowing for accurate tracking and state estimation even in the presence of significant noise.

## Key Features

- **Sensor Fusion:** Integrates data from both laser and radar sensors to provide a more accurate state estimation.
- **Nonlinear State Prediction:** Utilizes the unscented transform to handle nonlinear motion and measurement models effectively.
- **Noise Handling:** Configurable noise parameters for both process noise and measurement noise, enabling robust performance in various environments.
- **State Initialization:** Automatically initializes the state vector based on the first incoming measurement, whether from a laser or radar sensor.
- **Angle Normalization:** Ensures angles are normalized within the range of \(-\pi\) to \(\pi\), avoiding issues related to angle wrapping.



The main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ukf_highway


<img src="media/ukf_highway.png" width="700" height="400" />

`main.cpp` is using `highway.h` to create a straight 3 lane highway environment with 3 traffic cars and the main ego car at the center. 
The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the 
other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car's has
it's own UKF object generated for it, and will update each indidual one during every time step. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.

---

## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
 * PCL 1.2

# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program


---
## Results
[UKF-xy.png]: ./Results/UKF-xy.png
![alt text][UKF-xy.png]

Input file: obj_pose-laser-radar-synthetic-input.txt
Accuracy - RMSE:
[0.0690372, 0.0824925, 0.336984, 0.219354] is less than required [.09, .10, .40, .30]. 

AS below figure for NIS (Normalized Innovation Squared) show the NIS for lidar is below 95% line (5.991) for chi-square with degree of dreedom of 2. 
UKF-NIS for Lidar: 
[UKF-NIS-Lidar.png]: ./Results/UKF-NIS-Lidar.png
![alt_text][UKF-NIS-Lidar.png]

AS below figure for NIS (Normalized Innovation Squared) show the NIS for radar is below 95% line (7.815) for chi-square with degree of dreedom of 3. 
UKF-NIS for Radar: 
[UKF-NIS_Radar.png]: ./Results/UKF-NIS-Radar.png
![alt_text][UKF-NIS-Radar.png]


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`

## For more derails regarding Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/c3eb3583-17b2-4d83-abf7-d852ae1b9fff/concepts/f437b8b0-f2d8-43b0-9662-72ac4e4029c1)
for instructions and the project rubric.

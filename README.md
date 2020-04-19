# Implementation of Kalman Filter.
This package is to implement Kalman Filter.  

## Requirements
1. Eigen3

## Content
1. Implementation of Kalman Filter basic functions  
2. Three model cases implementation of Extended Kalman Filter (EKF)  
3. Visualization script

## 1. Implementation of Kalman Filter basic functions
The functions are defined in **~/kalman_filter/include/kalman.hpp** and **~/kalman_filter/src/kalman.cpp**  

## 2. Three model cases implementation of Extended Kalman Filter
Here we use radar data for example.  
Radar provides the range **r**, the azimuth **a**, the RCS value, and the range rate **r'** for each target measurements.  
In **~/kalman_filter/data.csv**, the information of one target measurement is provided with the time stamp.  
The test cases will use the provided data to implement EKF, and are all consided as constant velocity motion.  

* Compile the code
```
cd ~/kalman_filter
cmake ../kalman_filter
make
```

### Case I:
In case I, the state model is defined as **[x, y, vx, vy]**, where **(x, y)** is the position of the measurement and **(vx, vy)** is the velocity.  
The observation model is defined as **[x, y]**, where **(x, y)** is calculated from **rcos(a)** and **rsin(a)**.  
```
cd ~/kalman_filter
./devel/lib/kalman_filter/kalman_1
```

* Result:
<img src="https://github.com/tom13133/kalman_filter/blob/master/images/kalman_1.svg" width="1000">

### Case II:
In case II, the state model is defined as **[r, a, r', a']**, where **(r, a)** are the range and the azimuth. **(r', a')** are range rate and azimuth rate.  
The observation model is defined as **[r, a, r']**.  
```
cd ~/kalman_filter
./devel/lib/kalman_filter/kalman_2
```
* Result:
<img src="https://github.com/tom13133/kalman_filter/blob/master/images/kalman_2.svg" width="1000">

### Case III:
In case I, the state model is defined as **[x, y, vx, vy]**, where **(x, y)** is the position of the measurement and **(vx, vy)** is the velocity.  
The observation model is defined as **[x, y, r']**, where **(x, y)** is calculated from **rcos(a)** and **rsin(a)**, **r'** is range rate.  
```
cd ~/kalman_filter
./devel/lib/kalman_filter/kalman_3
```
* Result:
<img src="https://github.com/tom13133/kalman_filter/blob/master/images/kalman_3.svg" width="1000">

### Conclusion:
The possible reasons to affect the performance:  
1. The type of Motion model (i.e. constant velocity)  
2. The definition of state and observation model.  
3. The uncertainty **Q**, **R** of the model. (**This is important.**)  


## 3. Visualization script
To visualize the performance, **~/kalman_filter/result.py** is used.  
The script draw the radar measurements before and after EKF is applied (saved in **~/kalman_filter/data_out.csv**) and compare them with ground truth (saved in **~/kalman_filter/ground_truth.csv**).  
The lidar measurements are used as ground truth, and the transformation between radar and lidar has been defined in the script **result.py**.  
After running the script, the above figures are shown on the screen.  

* In **~/kalman_filter/data**, there is another data for testing.  





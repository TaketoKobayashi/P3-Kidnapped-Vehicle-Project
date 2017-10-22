## **Kidnapped Vehicle Project**

The goals / steps of this project are the following:

* Your robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

In this project I will implement a 2 dimensional particle filter in C++. Your particle filter will be given a map and some initial localization information (analogous to what a GPS would provide). At each time step your filter will also get observation and control data..

[//]: # (Image References)
[image1]: ./outputs/1.png

Final result:

![alt text][image1]

---

## Particle Filter Process


#### 1. Read the range measurement data into program.

This is provided in main.cpp file.

#### 2. Particle Filter Process:

##### 1. initialize our position from GPS input.

I initialized with 144 particles, take consideration: standard variance is between -0.3 ~ 0.3, so I divided the error space 0.6*0.6 into 12*12 grid, each grid's resolution is 0.05*0.05, so the particles number is 144. Use `normal_distribution` to generate particles.

##### 2. Predict the particles state, add control input(yaw rate and velocity).

In this step, we should consider the 2 situations, yaw_rate is zero and not zero.

##### 3. Update the particles weights.

During the update step, we update our particle weights using map landmark positions and feature measurements, and we did data association and set association.

##### 4. Resample the particles weights.

During resampling, we resample M times (M is range of 0 to length_of_particleArray) drawing a particle i (i is the particle index) proportional to its weight, and according to Sebastian's course to implement the resampling wheel algorithm.

The new set of particles represents the Bayes filter posterior probability. We now have a refined estimate of the vehicles position based on input evidence.


## Test Implementation:

RMSE screenshot as below:

![alt text][image1]





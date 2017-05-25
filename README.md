# MPC Control
## Self-Driving Car Engineer Nanodegree Program

Controlling steering and throttle input via MPC Control to drive a vehicle as fast as safely possible around a track in the Udacity Simulator.
The vehicle has a path along the center of the track as input.

## Result (max speed run): [YouTube video: MPC Control](https://www.youtube.com/watch?v=5UmWq_d3WW4)

## Model

The model state is:

* px: X-position of the vehicle in the forward direction
* py: Y-position of the vehicle in the lateral direction
* psi: Orientation of the vehicle
* v: Velocity of the vehicle

The actuators of the vehicle are:

* deltaPsi: Steering angle
* a: Acceleration

The update equations for our model used to predict future states are:
* px(t+1)  = px(t) + v(t) \* cos(psi(t)) \* dt
* py(t+1)  = py(t) + v(t) \* sin(psi(t)) \* dt
* psi(t+1) = psi(t) + v(t) / Lf \* deltaPsi \* dt
* v(t+1)   = v(t) + a \* dt;

Where `dt` is the timestep between predictions and `Lf` is the distance between the front and the center of gravity of the vehicle, which determines its turning radius.  
Within the simulator, the slip angle of the vehicle is zero as long as it is on the paved road, so cornering forces are not taken into account.

## Timesteps and Frequency

The final timestep value for the controller is 15 with a frequency of 0.12 seconds.

I started out with 8 steps at 0.1 seconds which provided decent results at low speeds, but quickly failed when the velocity was increased.
The car reacted too slowly to directional changes in the road and would run out of bounds on sharp turns. Additionally, the vehicle would tend to oscillate on fast straights and not deal well with latency (see below).

If either parameter was set too high, the vehicle would tend to drive too conservatively as it was planning its path down the road too far in advance to still maintain a high speed.

## Cost Function

In order to find the path with the lowest associated cost while maintaining a fairly fast speed, the following parameters where penalized with different factors
* Cross Track Error^2 * `2000`
* Orientation Error^2 * `1500`
* Deviation from reference velocity^2 * `1`. Reference velocity is 95 mp/h for safe laps around the track and 100 mp/h for "race trim".
* Use of steering actuator^2 * `20000`
* Use of acceleration actuator^2 * `1`
* Difference of sequential actuations for steering^2 * `2`
* Difference of sequential actuations for acceleration^2 * `1`

Giving considerably more weight to the cross track and orientation errors kept the vehicle on the road well at low speeds.
However, with increasing speed, oscillation would become worse and eventually kick the vehicle off the road. Severly penalizing the use of the steering actuator resulted in overall more smooth steering inputs and a stable drive.
A slightly larger penalty for the difference in sequential steering actuations was applied to keep constant cycles of full acceleration and braking at bay if the vehicle's reference speed is set below 100.

## Polynomial Fitting and MPC Preprocessing

The provided waypoints are transformed into vehicle space and then fitted to a polynomial. I use a 2-dimensional polynomial to fit the path since it gave more robust results than a 3-dimensional fit. With the planned path in vehicle space, the initial x/y position and orientation of the vehicle can be set to zero and fed into the MPC solver as part of the state vector.

## Dealing with Latency

In order to simulate a system that is closer to real-life, latency of 0.1 seconds between a cycle of the MPC controller and the actual actuation was artificially introduced.  
Predictably, this had two distinct effects:
* it affected driving at high speeds towards tight turns the most - the vehicle would react too late to recognizing a sharp turn and run off-road
* any existing oscillation was amplified


 To account for this, I made the model drive more conservatively by increasing the amount of track that it looks ahead via the timestep parameter (see above). Additionally, I further increased the penalty to the use of the steering actuator.

## Dependencies

* cmake >= 3.5
* make >= 4.1
* gcc/g++ >= 5.4
* uWebSockets == 0.14, but the master branch will probably work just fine
* Ipopt
* CppAD
* Eigen. This is already part of the repo so you shouldn't have to worry about it.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`

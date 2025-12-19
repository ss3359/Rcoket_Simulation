# Rcoket_Simulation
This is a 3D simulation of a rocket ascending to the moon using C++ code and the Runge-Kutta algorithm for updating the position and velocity. 
Below is the steps needed to run the code and the explanations on how the coode workd 


## Features
- Newtonian point-mass rocket model
- Thrust, gravity, drag, and wind forces
- Exponential atmosphere model
- Explicit RK4 time stepping
- Modular force decomposition

## Equations of Motion
dr/dt = v  
dv/dt = F(r, v, t) / m (Acceleration) 

## Build
g++ main.cpp MONTEREY_JACK.cpp -o rocket

## Run
./rocket



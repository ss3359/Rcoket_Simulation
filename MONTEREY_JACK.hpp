//  MONTEREY_JACK.hpp
//  CHEESE
//
//  Created by Owner on 12/12/25.
//

#ifndef MONTEREY_JACK_hpp
#define MONTEREY_JACK_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <vector>
#include <string>
using namespace std;




/*
    We are going to design a rocket simulation using
    a Runge-Kutta stepping method
*/

class Rocket{
    private:
        double G = 6.6743e-11; // m^3 kg^-1 s^-2 Gravitational Constant
        double m = 100; //Initial mass of the rocket (kg)
        double M = 5.972e24; // Mass of the Earth (kg)
        double Cd = 0.2; // Drag Coefficient For The Rocket
        double A= 0.0314; // Cross Sectional Area, πr^2 which r=0.1 meters
        double dt = 0.01; //Time Step

    public:
        vector<double> position=vector<double>(3);
        vector<double> velocity=vector<double>(3);
        
        Rocket(double px, double py, double pz,double vx, double vy, double vz){
            position[0]=px; position[1]=py; position[2]=pz;
            velocity[0]=vx; velocity[1]=vy; velocity[2]=vz;
        }

        //Calculate the length of vector
        double vector_length( const vector<double> &r){
            double sum=0;
            for(double r_c: r){
                sum += r_c*r_c;
            }
            return sqrt(sum);
        }

        //Calculate The Forces Acting On The Rocket As It is moving.
    
        double AirDensity(const vector<double> &pos, const vector<double> &vel);
    
        vector<double> Wind(const vector<double> &pos,const vector<double> &vel);
    
        vector<double> Weight(const vector<double> &pos,const vector<double> &vel);
        vector<double> RelativeVel(const vector<double> &pos,const vector<double> &vel);
    
    
        vector<double> Lift( const vector<double>& pos, const vector<double> &vel);
    
        vector<double> Thrust(const vector<double> &pos,const vector<double> &vel,double t);
    
        vector<double> Drag(const vector<double> &pos,const vector<double> &vel);
    
    vector<double> AddVectors(const vector<double> &a,const vector<double> &b,const vector<double> &c,const vector<double> &d);
        
        vector<double> Acceleration(const vector<double>& pos, const vector<double>& vel, double t);
        void UpdateVelocityAndPosition(vector<double> &pos, vector<double> &vel, double t);
};
double Fuel_Burn(double t);
vector<double> operator+(const vector<double>& a,const vector<double> &b);
vector<double> operator-(const vector<double> &a,const vector<double> &b);
vector<double> operator*(double c,const vector<double> &a);
vector<double>& operator+=(vector<double>&a,const vector<double>&b);

vector<double> AddVectors(const vector<double> &a, const vector<double> &b,const vector<double> &c,const vector<double> &d);

#endif /* MONTEREY_JACK_hpp */

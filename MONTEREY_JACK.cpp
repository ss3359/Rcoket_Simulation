//
//  MONTEREY_JACK.cpp
//  CHEESE
//
//  Created by Owner on 12/12/25.
//

#include "MONTEREY_JACK.hpp"

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <vector>
#include <string>
using namespace std;


vector<double> Rocket:: Weight(const vector<double>& pos, const vector<double> &vel){
    double g = 9.8; //gravitatinoal acceleration ms^-2
    return {0,0,-m*g};

};

vector<double> Rocket:: Wind(const vector<double> & pos,const vector<double>& vel){
    
    vector<double> result(pos.size(),0);
    
    if(pos[2] >=0 && pos[2]<= 1000)
        result[0] = 8*vel[2]/1000;
    else if(pos[2]>1000 && pos[2]<=5000)
        result[0] = 12*((vel[2]-1000)/4000) + 8;
    else
        result[0]=20;
    result[1]=0;
    
    return result;
    
}
vector<double> Rocket:: RelativeVel(const vector<double> & pos, const vector<double> &vel){
    
    vector<double> wind=Wind(pos,vel);
    return vel-wind;
}

vector<double> Rocket:: Drag(const vector<double> &pos,const vector<double>& vel){
    vector<double> result(pos.size());
    
    
    vector<double> rel_vel=RelativeVel(pos,vel);
    double scalar=-(0.5)*AirDensity(pos,vel)*Cd*A*vector_length(rel_vel);
    
    return scalar*rel_vel;
};
vector<double> Rocket:: Thrust(const vector<double> &pos, const vector<double>& vel,double t){
    vector<double> result(vel.size());
    
    double scalar= Fuel_Burn(t);
    vector<double> n = {0,0,1};
    result=scalar*n;
    
    return result;
};

vector<double> Rocket::Acceleration(const vector<double> & pos, const vector<double> & vel, double t){
    vector<double> result;
    vector<double> weight=Weight(pos,vel), drag=Drag(pos,vel);
    vector<double> thrust=Thrust(pos,vel,t);
        
    vector<double> DT=drag+thrust;
    result=weight+DT;
    
    return (1/m)*result;
    
}

void Rocket:: UpdateVelocityAndPosition( vector<double>& pos, vector<double> &vel, double t){
    
    /*
     The RK4 method will correspond to these equations
     drdt(r, v)= v
     dvdt(r, v)= F/m where F is the sum of the forces

    State:
    r' = v
    v' = a(r, v, t)

    RK4 scheme:
        J = velocity stages
        K = acceleration stages

     */

// TODO:
// - Add variable mass during burn
// - Add gravity turn / pitch program
// - Output trajectory to CSV for plotting

    int iter = 100;
    
    for(int i =1; i <=iter; i++){
        vector<double> K1 = Acceleration(pos,vel,t);
        vector<double> J1 = vel;
        
        vector<double> step_one_J=(dt/2)*K1;
        step_one_J+=pos;
        
        vector<double> step_one_K=(dt/2)*J1;
        step_one_K+=vel;
        
        
        
        vector<double> K2=Acceleration(step_one_J,step_one_K,t+(dt/2));
        vector<double> J2 =step_one_J;
        
        vector<double> step_two_J=(dt/2)*K2;
        step_two_J+=pos;
        
        vector<double> step_two_K=(dt/2)*J2;
        step_two_K+=vel;
        
        vector<double> K3 = Acceleration(step_two_J,step_two_K,t+(dt/2));
        vector<double> J3 = step_two_J;
        
        
        vector<double> step_three_J=(dt)*K3;
        step_three_J+=pos;
        
        vector<double> step_three_K=(dt)*J3;
        step_three_K+=vel;
        
        vector<double> K4 = Acceleration(step_three_J,step_three_K,t+dt);
        vector<double> J4= step_three_K;
        

        
        pos+=(dt/6)*(J1+2*J2+2*J3+J4);
        vel+=(dt/6)*(K1+2*K2+2*K3+K4);
        
        cout<<"New Postion: <"<<pos[0]<<","<<pos[1]<<","<<pos[2]<<">"<<endl;
        cout<<"New Velocity: <"<<vel[0]<<","<<vel[1]<<","<<vel[2]<<">"<<endl;
        cout<<endl;

    }
    
}


double Rocket:: AirDensity(const vector<double> &pos,const vector<double> &vel){
    
    return 1.225*exp(-pos[2]/8500);
}

vector<double> operator+(const vector<double> &a,const vector<double> &b){
    vector<double> result(a.size());
    
    for(int i=0; i<a.size(); i++){
        result[i]=a[i] + b[i];
    }
    return result;
}
vector<double> operator-(const vector<double> &a,const vector<double> &b){
    vector<double> result(a.size());
    
    for(int i=0; i<a.size(); i++){
        result[i]=a[i] - b[i];
    }
    return result;
}
vector<double> operator*(double c,const vector<double> &a){
    vector<double> result(a.size());
    
    for(int i=0; i<a.size(); i++){
        result[i]=c*a[i];
    }
    return result;
}



vector<double>& operator+=(vector<double>&a,const vector<double>&b){
    
    for(int i=0; i<a.size(); i++){
        a[i]+=b[i];
    }
    return a;
}



double Fuel_Burn(double t){
    if(0<=t && t<=10){
        return 3000;
    }
    else
        return 0;
}

vector<double> Rocket:: AddVectors(const vector<double> &a, const vector<double> &b,const vector<double> &c,const vector<double> &d){
    vector<double> result1, result2, sum1, sum2, total;
    
    
    result1=(2.0)*b; result2=(2.0)*c;
    sum1=a+result1; sum2=result2+d;
    
    total=sum1+sum2;
    
    return (dt/6)*total;
}

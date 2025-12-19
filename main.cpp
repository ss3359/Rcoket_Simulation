//
//  main.cpp
//  CHEESE
//
//  Created by Owner on 11/10/25.
//

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <thread>

#include "MONTEREY_JACK.hpp"

using namespace std;


int main(){
    Rocket r(0, 0, 0, 0, 0, 0);
    vector<double> P = r.position, V=r.velocity;
    double t = 0.0;
    r.UpdateVelocityAndPosition(P, V, t);
    return 0;
}


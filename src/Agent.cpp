//
//  Agent.cpp
//  mySketch
//
//  Created by Max Jappert on 07/05/2024.
//

#include "Agent.hpp"
#include <vector>

// Constructor definition
Agent::Agent(std::vector<float> loc_, std::vector<float> dir_) {
    // Initialize member variables if needed
    loc = loc_;
    dir = dir_;
}

// Constructor definition
Agent::Agent() {
}

// Destructor definition
Agent::~Agent() {
    // Clean up resources if needed
}

void Agent::rotate(float theta) {
    float cosTheta = cos(theta);
    float sinTheta = sin(theta);
    std::vector<float> rotatedVec(2);
    rotatedVec[0] = dir[0] * cosTheta - dir[1] * sinTheta; // x' = x*cos(theta) - y*sin(theta)
    rotatedVec[1] = dir[0] * sinTheta + dir[1] * cosTheta; // y' = x*sin(theta) + y*cos(theta)
    
    dir = rotatedVec;
}


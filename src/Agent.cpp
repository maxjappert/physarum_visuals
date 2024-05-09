//
//  Agent.cpp
//  mySketch
//
//  Created by Max Jappert on 07/05/2024.
//

#include "Agent.hpp"
#include <vector>
#include <glm/glm.hpp>  // Include the main GLM header
#include <glm/trigonometric.hpp>  // Include this for trigonometric functions

// Constructor definition
Agent::Agent(glm::vec2 loc_, glm::vec2 dir_) {
    loc = loc_;
    dir = dir_;
}

Agent::Agent(float locx, float locy, float dirx, float diry) {
    loc.x = locx;
    loc.y = locy;
    dir.x = dirx;
    dir.y = diry;
    //glm::vec2 loc(locx, locy);
    //glm::vec2 dir(dirx, diry);
}

// Constructor definition
//Agent::Agent() {
//}

// Destructor definition
Agent::~Agent() {
    // Clean up resources if needed
}

void Agent::rotate(float theta) {
    float cosTheta = glm::cos(theta);
    float sinTheta = glm::sin(theta);
    
    glm::vec2 rotatedVec;
    
    rotatedVec.x = dir.x * cosTheta - dir.y * sinTheta;
    rotatedVec.y = dir.x * sinTheta + dir.y * cosTheta;
    
    dir = rotatedVec;
}


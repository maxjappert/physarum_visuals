//
//  Agent.hpp
//  mySketch
//
//  Created by Max Jappert on 07/05/2024.
//

#ifndef Agent_hpp
#define Agent_hpp

#include <stdio.h>
#include <vector>
#include <glm/glm.hpp>  // Include the main GLM header

class Agent {
public:
    // Constructor
    //Agent();
    Agent(glm::vec2 loc_, glm::vec2 dir_);
    Agent(float locx, float locy, float dirx, float diry);

    // Destructor
    ~Agent();

    // Example method
    glm::vec2 loc;
    glm::vec2 dir;
    
    void rotate(float theta);
    

private:
    // Member variable
};

#endif /* Agent_hpp */

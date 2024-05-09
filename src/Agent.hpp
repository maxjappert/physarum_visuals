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

class Agent {
public:
    // Constructor
    Agent();
    Agent(std::vector<float> loc_, std::vector<float> dir_);

    // Destructor
    ~Agent();

    // Example method
    std::vector<float> loc;
    std::vector<float> dir;
    
    void rotate(float theta);
    

private:
    // Member variable
};

#endif /* Agent_hpp */

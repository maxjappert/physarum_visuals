#include "ofApp.h"
#include "Agent.hpp"
#include <glm/glm.hpp>
#include <glm/trigonometric.hpp>
#include <glm/vector_relational.hpp>
#include <glm/gtc/random.hpp>  // For glm::linearRand

const int numAgents = 1000000;
Agent** agents = new Agent*[numAgents];
float** trailMap;
float stepSize = 1;
int width;
int height;

float decayT = 2;
float sensoryAngle = 3.141/8;
float rotationAngle = 3.141/4;
float sensorOffset = 9;
float deposit = 5;
float maxDeposit = 1000;

GLuint vbo;

int kernelSize = 3;
float sigma = 0.5;

ofImage image;

ofPixels pixels;

ofTexture trailMapTexture;


// Helper function to create a Gaussian kernel
std::vector<std::vector<float>> createGaussianKernel(int kernelSize, float sigma) {
    std::vector<std::vector<float>> kernel(kernelSize, std::vector<float>(kernelSize));
    float sum = 0.0;
    int halfSize = kernelSize / 2;
    for (int i = -halfSize; i <= halfSize; i++) {
        for (int j = -halfSize; j <= halfSize; j++) {
            float value = exp(-(i * i + j * j) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
            kernel[i + halfSize][j + halfSize] = value;
            sum += value;
        }
    }

    // Normalize the kernel
    for (int i = 0; i < kernelSize; i++) {
        for (int j = 0; j < kernelSize; j++) {
            kernel[i][j] /= sum;
        }
    }
    return kernel;
}

// Function to apply Gaussian blur to a 2D array
void applyGaussianBlur(float** array, int width, int height, int kernelSize, float sigma) {
    std::vector<std::vector<float>> kernel = createGaussianKernel(kernelSize, sigma);
    int halfSize = kernelSize / 2;

    // Create a temporary array to store the blurred values
    float** tempArray = new float*[width];
    for (int i = 0; i < width; ++i) {
        tempArray[i] = new float[height];
        std::copy(array[i], array[i] + height, tempArray[i]);
    }

    for (int x = halfSize; x < width - halfSize; ++x) {
        for (int y = halfSize; y < height - halfSize; ++y) {
            float sum = 0.0;
            for (int i = -halfSize; i <= halfSize; ++i) {
                for (int j = -halfSize; j <= halfSize; ++j) {
                    sum += tempArray[x + i][y + j] * kernel[halfSize + i][halfSize + j];
                }
            }
            array[x][y] = sum;
        }
    }

    // Delete the temporary array
    for (int i = 0; i < width; ++i) {
        delete[] tempArray[i];
    }
    delete[] tempArray;
}

// Function to create a normalized 2D vector using float
glm::vec2 createRandomNormalized2DVector() {
    float angle = glm::linearRand(0.0f, 2 * glm::pi<float>());
    return glm::vec2(glm::cos(angle), glm::sin(angle));
}

//--------------------------------------------------------------
void ofApp::setup(){
    //ofSetVerticalSync(true);
    ofBackground(20);
    glColor4f(1, 1, 1, 1);
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glPointSize(0.1f);
    
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, numAgents * sizeof(float) * 2, NULL, GL_STATIC_DRAW); // Use GL_DYNAMIC_DRAW if data changes often
    glVertexPointer(2, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);

    width = (int)ofGetWidth();
    height = (int)ofGetHeight();
    
    for (int i = 0; i < numAgents; i++) {
        glm::vec2 newLoc(glm::linearRand(0, width-1), glm::linearRand(0, height-1));
        glm::vec2 newDir = createRandomNormalized2DVector();
        //agents[i] = new Agent(newLoc, newDir);
        agents[i] = new Agent(newLoc.x, newLoc.y, newDir.x, newDir.y);
    }
    
    trailMap = new float*[width];
    
    
    for (int i = 0; i < width; i++) {
        trailMap[i] = new float[height];
        for (int j = 0; j < height; j++) {
            trailMap[i][j] = 0;
        }
    }
    
    //pixels.allocate(width, height, OF_PIXELS_RGB);
    //trailMapTexture.allocate(width, height, GL_R32F);  // Allocate a texture with floating-point precision
}

bool isLegalCoords(int x, int y) {
    return x >= 0 && x < width && y >= 0 && y < height;
}

glm::vec2 getRotatedVec(glm::vec2 input, float theta) {

    float cosTheta = glm::cos(theta);
    float sinTheta = glm::sin(theta);
    
    glm::vec2 rotatedVec;
    
    rotatedVec.x = input.x * cosTheta - input.y * sinTheta;
    rotatedVec.y = input.x * sinTheta + input.y * cosTheta;
    
    return rotatedVec;
}

glm::vec2 cutVec(glm::vec2 vec) {
    if (vec.x < 0) {
        vec.x = 0;
    } else if (vec.x >= width-1) {
        vec.x = width-1;
    }
    
    if (vec.y < 0) {
        vec.y = 0;
    } else if (vec.y > height-1) {
        vec.y = height-1;
    }
    
    return vec;
}

float ofApp::map(float value, float inputMin, float inputMax, float outputMin, float outputMax, bool clamp) {
    if (fabs(inputMax - inputMin) < FLT_EPSILON){
        return outputMin; // avoid zero division
    } else {
        float outVal = ((value - inputMin) / (inputMax - inputMin) * (outputMax - outputMin) + outputMin);

        if(clamp){
            if(outputMax < outputMin){
                if(outVal < outputMax) outVal = outputMax;
                else if(outVal > outputMin) outVal = outputMin;
            }else{
                if(outVal > outputMax) outVal = outputMax;
                else if(outVal < outputMin) outVal = outputMin;
            }
        }
        return outVal;
    }
}


//--------------------------------------------------------------
void ofApp::update(){
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            if (trailMap[i][j] - decayT >= 0) {
                trailMap[i][j] -= decayT;
            } else {
                trailMap[i][j] = 0;
            }
        }
    }
        
    applyGaussianBlur(trailMap, width, height, kernelSize, sigma);
    
    //image.setFromPixels(pixels);
    //image.draw(0, 0);
    
    for (int i = 0; i < numAgents; i++) {
        // Motor stage
        agents[i]->loc += stepSize*agents[i]->dir;
        
        int newLocX = (int)round(agents[i]->loc.x);
        int newLocY = (int)round(agents[i]->loc.y);
        if (isLegalCoords(newLocX, newLocY)) {
            if (trailMap[newLocX][newLocY] + deposit < maxDeposit) {
                trailMap[newLocX][newLocY] += deposit;
            } else {
                trailMap[newLocX][newLocY] = maxDeposit;
            }
        } else {
            agents[i]->rotate(glm::linearRand(0.0f, 2*glm::pi<float>()));
        }
        
        // Sensory stage
        Agent agent = *agents[i];
        glm::vec2 sensorVecF = agents[i]->dir * sensorOffset;
        glm::vec2 F_vec = cutVec(agents[i]->loc + sensorVecF);
        
        float F = trailMap[(int)round(F_vec.x)][(int)round(F_vec.y)];
        
        glm::vec2 sensorVecFL = getRotatedVec(sensorVecF, -sensoryAngle);
        glm::vec2 FL_vec = cutVec(agents[i]->loc + sensorVecFL);
        
        //std::cout << (int)round(FL_vec.x) << std::endl;
        //std::cout << (int)round(FL_vec.y) << std::endl;
        
        float FL = trailMap[(int)round(FL_vec.x)][(int)round(FL_vec.y)];
        
        glm::vec2 sensorVecFR = getRotatedVec(sensorVecF, sensoryAngle);
        glm::vec2 FR_vec = cutVec(agents[i]->loc + sensorVecFR);
        
        float FR = trailMap[(int)round(FR_vec.x)][(int)round(FR_vec.y)];
        
        if (F > FL && F > FR) {
            continue;
        } else if (F < FL && F < FR) {
            if (glm::linearRand(0.0f, 1.0f) < 0.5) {
                agents[i]->rotate(-rotationAngle);
            } else {
                agents[i]->rotate(rotationAngle);
            }
        } else if (FL < FR) {
            agents[i]->rotate(rotationAngle);
        } else if (FL > FR) {
            agents[i]->rotate(-rotationAngle);
        }
    }
}

//--------------------------------------------------------------
void ofApp::draw(){
    
    // Update VBO data if necessary
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    float* vboData = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    for(int i = 0; i < numAgents; i++) {
        vboData[2*i] = agents[i]->loc.x;
        vboData[2*i+1] = agents[i]->loc.y;
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    // Draw the points
    glDrawArrays(GL_POINTS, 0, numAgents);

    // Cleanup if needed, though not necessary every frame
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    
    //ofSaveFrame();
}

//--------------------------------------------------------------
void ofApp::exit(){

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseScrolled(int x, int y, float scrollX, float scrollY){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}

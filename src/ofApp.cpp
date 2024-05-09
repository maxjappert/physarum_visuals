#include "ofApp.h"
#include "Agent.hpp"

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
std::vector<float> createRandomNormalized2DVector() {
    // Random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd;
    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd());
    // Normal distribution with mean 0 and standard deviation 1, specialized for float
    std::normal_distribution<float> dist(0.0f, 1.0f);

    // Generate two independent normally distributed random numbers (floats)
    float x = dist(gen);
    float y = dist(gen);

    // Calculate the magnitude of the vector (float precision)
    float magnitude = std::sqrtf(x*x + y*y);

    // Normalize the vector components to unit length
    x /= magnitude;
    y /= magnitude;

    // Create a vector to hold the normalized components (floats)
    std::vector<float> normalizedVector = {x, y};

    return normalizedVector;
}

//--------------------------------------------------------------
void ofApp::setup(){
    //ofSetVerticalSync(true);
    ofBackground(20);
    
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::uniform_real_distribution<> disDouble(0, ofGetWidth());
    
    for (int i = 0; i < numAgents; i++) {
        std::vector<float> randomVec = createRandomNormalized2DVector();
        agents[i] = new Agent({(float)disDouble(gen), (float)disDouble(gen)}, createRandomNormalized2DVector());
        //agents[i] = new Agent({100, 100}, createRandomNormalized2DVector());
    }
    
    width = (int)ofGetWidth();
    height = (int)ofGetHeight();
    trailMap = new float*[width];
    
    for (int i = 0; i < width; ++i) {
        trailMap[i] = new float[height];
        for (int j = 0; j < height; j++) {
            trailMap[i][j] = 0;
        }
    }
    
    pixels.allocate(width, height, OF_PIXELS_RGB);
}

bool isLegalCoords(int x, int y) {
    return x >= 0 && x < (int)ofGetWidth() && y >= 0 && y < (int)ofGetHeight();
}

// Function to generate a random float between a and b
float getRandomFloat(float a, float b) {
    // Static variables for the random engine and distribution
    static std::random_device rd;  // Obtain a seed from the system entropy device, or whatever is available
    static std::mt19937 gen(rd()); // Seed the generator
    static std::uniform_real_distribution<> dis(a, b); // Define the range [a, b]

    return dis(gen); // Generate and return the random float
}

std::vector<float> getRotatedVec(std::vector<float> input, float theta) {

    float cosTheta = cos(theta);
    float sinTheta = sin(theta);
    std::vector<float> rotatedVec(2);
    rotatedVec[0] = input[0] * cosTheta - input[1] * sinTheta; // x' = x*cos(theta) - y*sin(theta)
    rotatedVec[1] = input[0] * sinTheta + input[1] * cosTheta; // y' = x*sin(theta) + y*cos(theta)

    return rotatedVec;
}

std::vector<float> mult2dVector(std::vector<float> vec, float factor) {
    return {vec[0]*factor, vec[1]*factor};
}

int cutX(int x) {
    if (x < 0) {
        x = 0;
    } else if (x >= ofGetWidth()) {
        x = ofGetWidth()-1;
    }
    
    return x;
}

int cutY(int y) {
    if (y < 0) {
        y = 0;
    } else if (y >= ofGetHeight()) {
        y = ofGetHeight()-1;
    }
    
    return y;
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
    
}

//--------------------------------------------------------------
void ofApp::draw(){
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glPointSize(0.1f);
    glBegin(GL_POINTS);
    for(int i = 0; i < numAgents; i++) {
        glColor4f(1, 1, 1, 1);
        glVertex2f(agents[i]->loc[0], agents[i]->loc[1]);
    }
    glEnd();
        
    //for (int i = 0; i < numAgents; i++) {
    //    ofSetColor(255, 255, 255);
    //    ofDrawCircle(agents[i]->loc[0], agents[i]->loc[1], 0.1);
    //}
    
    #pragma omp parallel for
    for (int i = 0; i < width; i++) {
        #pragma omp parallel for
        for (int j = 0; j < height; j++) {
            int mappedValue = round(map(trailMap[i][j], 0, maxDeposit, 0, 255, true));
            //ofColor color(mappedValue, mappedValue, mappedValue); // Red color
            //pixels.setColor(i, j, color);
            //ofSetColor(mappedValue, mappedValue, mappedValue);
            //ofDrawRectangle(i, j, 1, 1);
            
            //std::cout << mappedValue << std::endl;
            //std::cout << trailMap[i][j] << std::endl;
            
            if (trailMap[i][j] >= 0) {
                trailMap[i][j] -= decayT;
                //ofSetColor(trailMap[i][j], trailMap[i][j], trailMap[i][j]);
                //ofDrawCircle(j, i, 1);
                //ofSetColor(255, 255, 255);
            }
        }
    }
        
    applyGaussianBlur(trailMap, width, height, kernelSize, sigma);
    
    image.setFromPixels(pixels);
    //image.draw(0, 0);
    
    #pragma omp parallel for
    for (int i = 0; i < numAgents; i++) {
        // Motor stage
        //agents[i]->rotate(0.01);
        agents[i]->loc[0] += stepSize*agents[i]->dir[0];
        agents[i]->loc[1] += stepSize*agents[i]->dir[1];
        int newLocX = round(agents[i]->loc[0]);
        int newLocY = round(agents[i]->loc[1]);
        if (isLegalCoords(newLocX, newLocY)) {
            if (trailMap[newLocX][newLocY] + deposit < maxDeposit) {
                trailMap[newLocX][newLocY] += deposit;
            } else {
                trailMap[newLocX][newLocY] = maxDeposit;
            }
        } else {
            agents[i]->rotate(getRandomFloat(0, 2*3.141));
        }
        
        // Sensory stage
        std::vector<float> sensorVecF = mult2dVector(agents[i]->dir, sensorOffset);
        int F_x = cutX(round(agents[i]->loc[0] + sensorVecF[0]));
        int F_y = cutY(round(agents[i]->loc[1] + sensorVecF[1]));
        
        float F = trailMap[F_x][F_y];
        
        std::vector<float> sensorVecFL = getRotatedVec(sensorVecF, -sensoryAngle);
        int FL_x = cutX(round(agents[i]->loc[0] + sensorVecFL[0]));
        int FL_y = cutY(round(agents[i]->loc[1] + sensorVecFL[1]));
        
        float FL = trailMap[FL_x][FL_y];
        
        std::vector<float> sensorVecFR = getRotatedVec(sensorVecF, sensoryAngle);
        int FR_x = cutX(round(agents[i]->loc[0] + sensorVecFR[0]));
        int FR_y = cutY(round(agents[i]->loc[1] + sensorVecFR[1]));
        
        float FR = trailMap[FR_x][FR_y];
        
        if (F > FL && F > FR) {
            continue;
        } else if (F < FL && F < FR) {
            if (getRandomFloat(0, 1) < 0.5) {
                agents[i]->rotate(-rotationAngle);
            } else {
                agents[i]->rotate(rotationAngle);
            }
        } else if (FL < FR) {
            agents[i]->rotate(rotationAngle);
        } else if (FL > FR) {
            agents[i]->rotate(-rotationAngle);
        }
        
        //ofSetColor(255, 0, 0);
        //ofDrawCircle(FR_x, FR_y, 1);
        //ofDrawCircle(FL_x, FL_y, 1);
        //ofSetColor(255, 255, 255);
    }
    
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

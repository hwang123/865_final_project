
#ifndef __ROI__h
#define __ROI__h


#include "Image.h"
// #include "basicImageManipulation.h"
#include <iostream>
#include <cmath>

using namespace std;

class ROI {
public:
    ROI(
        Image img, 
        Image mask,
        float score,
        int x1,
        int x2,
        int y1,
        int y2
    );

    int getWidth()             const { return width; } 
    int getHeight()            const { return height; }
    int getChannels()          const { return channels; } 
    float getPredictionScore() const { return predictionScore; }
    Image& getMask()             { return &mask; }
    Image& getImg()              { return &img; }

    void writeMask();
    void writeImg();

    // vector<int> getBoundingBox() const {
    //     vector<int> ret{x1, x2, y1, y2};
    //     return ret;
    // };

private:
    // Images
    Image img;
    Image mask;
    // Dimensions
    int width;
    int height;
    int channels;
    // Bounding box coordinates
    int x1;
    int x2;
    int y1;
    int y2;
    // Prediction score
    float predictionScore;

};

 
#endif

/* -----------------------------------------------------------------
 * File:    a5_main.cpp
 * Author:  Michael Gharbi <gharbi@mit.edu>
 * Created: 2015-09-30
 * -----------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------*/


#include "Image.h"
#include "basicImageManipulation.h"
#include "morphing.h"
#include "ROI.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

vector<string> split(string const &input) { 
    istringstream buffer(input);
    vector<string> ret;

    copy(istream_iterator<string>(buffer), 
        istream_iterator<string>(),
        back_inserter(ret));
    return ret;
}

vector<ROI> readFromFile(string fileName) {
    vector<ROI> rois; 

    ifstream inFile;
    inFile.open("../res/" + fileName);

    string line;

    // Get the width and height of the image
    getline(inFile, line);
    vector<string> dimensions = split(line);

    int width = stoi(dimensions[0]);
    int height = stoi(dimensions[1]);
    cout << "Image has size of " << width << ", " << height << endl;

    // Get number of ROI
    getline(inFile, line);
    int numROI = stoi(line);
    cout << "Image has " << numROI << " regions of interest." << endl;

    Image reference("./Input/dogs.png");

    for (int i = 0; i < numROI; ++i) {
        // Get the class ID and prediction score
        getline(inFile, line);
        vector<string> roiInfo = split(line);
        int classId = stoi(roiInfo[0]);
        float predictionScore = stof(roiInfo[1]);

        // Get the ROI bounding box
        getline(inFile, line);
        vector<string> boundingBox = split(line);

        int x1 = stoi(boundingBox[1]);
        int x2 = stoi(boundingBox[3]);
        int y1 = stoi(boundingBox[0]);
        int y2 = stoi(boundingBox[2]);

        int roiWidth = x2-x1+1;
        int roiHeight = y2-y1+1;

        // Create an Image to store the silhouette of the ROI
        Image mask = Image(roiWidth, roiHeight, 1);
        // Create an Image to store the ROI content itself
        Image img = Image(roiWidth, roiHeight, 3);

        for (int y = 0; y < roiHeight; ++y) {
            // Get the first line of the mask, which represents a row in the ROI
            getline(inFile, line);
            vector<string> roiRow = split(line);

            for (int x = 0; x < roiWidth; ++x) {
                // Get the actual pixel value
                float pixelVal = stof(roiRow[x]);
                mask(x, y) = pixelVal;

                for (int z = 0; z < img.channels(); ++z) {
                    img(x, y, z) = reference(x1 + x, y1 + y, z);
                }
            }   
        }
        ROI roi(img, mask, predictionScore, x1, x2, y1, y2);
        roi.writeImg("Output/dog" + to_string(i) + ".png");
        rois.push_back(roi);
    }
    inFile.close();

    return rois;
}

int main() {
    // Parse file and get regions of interest
    vector<ROI> rois = readFromFile("coco_res.txt");
}


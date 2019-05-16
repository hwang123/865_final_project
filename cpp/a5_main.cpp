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

vector<ROI> readFromFile(string fileName, Image &reference, int paddingRadius) {
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

    for (int i = 0; i < numROI; ++i) {
        // Get the class ID and prediction score
        getline(inFile, line);
        vector<string> roiInfo = split(line);
        int classId = stoi(roiInfo[0]);
        float predictionScore = stof(roiInfo[1]);

        // Get the ROI bounding box
        getline(inFile, line);
        vector<string> boundingBox = split(line);

        int x1 = stoi(boundingBox[1]) + paddingRadius;
        int y1 = stoi(boundingBox[0]) + paddingRadius;

        int x2 = stoi(boundingBox[3]) + paddingRadius;
        int y2 = stoi(boundingBox[2]) + paddingRadius;

        int roiWidth = x2-x1+1;
        int roiHeight = y2-y1+1;

        cout << "ROI dimensions: " << roiWidth << "," << roiHeight << endl;

        int paddedWidth = roiWidth + 2*paddingRadius;
        int paddedHeight = roiHeight + 2*paddingRadius;

        // Create an Image to store the silhouette of the ROI
        Image mask = Image(paddedWidth, paddedHeight, 1);
        // Create an Image to store the ROI content itself
        Image img = Image(paddedWidth, paddedHeight, 3);

        // for (int y = 0; y < roiHeight; ++y) {
        //     // Get the first line of the mask, which represents a row in the ROI
        //     getline(inFile, line);
        //     vector<string> roiRow = split(line);

        //     for (int x = 0; x < roiWidth; ++x) {
        //         // Get the actual pixel value
        //         float pixelVal = stof(roiRow[x]);
        //         mask(x, y) = pixelVal;

        //         int refIdxX = x1 - paddingRadius + x;
        //         int refIdxY = y1 - paddingRadius + y;

        //         for (int z = 0; z < img.channels(); ++z) {
        //             img(x, y, z) = reference(refIdxX, refIdxY, z);
        //         }
        //     }   
        // }

        for (int y = 0; y < paddedHeight; ++y) {
            // Get the first line of the mask, which represents a row in the ROI
            bool heightInMask = y >= paddingRadius && y < paddedHeight - paddingRadius;
            vector<string> roiRow;
            if (heightInMask) {
                getline(inFile, line);
                roiRow = split(line);       
            }

            for (int x = 0; x < paddedWidth; ++x) {
                float maskVal = 0.0f;
                bool widthInMask = x >= paddingRadius && x < paddedWidth - paddingRadius;

                if (heightInMask && widthInMask) {
                    maskVal = stof(roiRow[x-paddingRadius]);
                }
                mask(x, y) = maskVal;

                int refIdxX = x1 - 2*paddingRadius + x;
                int refIdxY = y1 - 2*paddingRadius + y;

                for (int z = 0; z < img.channels(); ++z) {
                    img(x, y, z) = reference(refIdxX, refIdxY, z);
                }
            }   
        }

        ROI roi(img, mask, predictionScore, x1 - 2*paddingRadius, x2, y1 - 2*paddingRadius, y2);
        roi.writeImg("Output/zebra/roi.png");

        rois.push_back(roi);
    }
    inFile.close();

    return rois;
}


Image compose(Image &base, vector<ROI> &rois, vector<Image> &styledImgs) {
    Image overlayed = base;

    for (int i = 0; i < styledImgs.size(); ++i){
        Image styled = styledImgs[i];
        ROI roi = rois[i];

        // Scale the styled image to fit the ROI
        float factor = float(roi.getWidth()) / styled.width();  
        cout << factor << endl;
        styled.write("./Output/highway/debug/styled.png");
        Image styledScaled = scaleBicubic(styled, factor, 3, 3);
        styledScaled.write("./Output/highway/debug/styledScaled.png");

        overlayed = overlayStyleTransferBoundary(overlayed, styledScaled, roi);
    }

    return overlayed;
}

void transferZebra(string styleName) {
    Image reference("./Input/zebra/zebras.png");
    vector<ROI> rois = readFromFile("zebra_data.txt", reference, 0);

    vector<Image> styledImgs;
    styledImgs.push_back(Image("./Input/zebra/zebra_" + styleName + ".png"));

    Image composed = compose(reference, rois, styledImgs);
    composed.write("./Output/zebra/composed_" + styleName + ".png");
}

void transferHighway(string styleName1, string styleName2, string styleName3) {
    Image reference("./Input/highway/highway.png");
    vector<ROI> rois = readFromFile("highway_data.txt", reference, 25);

    rois[0].getMask().write("./Output/highway/debug/mask0.png");
    vector<Image> styledImgs;
    styledImgs.push_back(Image("./Input/highway/" + styleName1 + "/roi0.png"));
    styledImgs.push_back(Image("./Input/highway/" + styleName2 + "/roi1.png"));
    styledImgs.push_back(Image("./Input/highway/" + styleName3 + "/roi2.png"));

    Image composed = compose(reference, rois, styledImgs);
    composed.write("./Output/highway/composed_" + styleName1 + ".png");    
}

int main() {
    // -----ZEBRA-----
    // transferZebra("popart");
    // transferZebra("sketch");
    // transferZebra("popart_bw");
    transferHighway("special", "special", "special");
}


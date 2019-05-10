#include "ROI.h"

using namespace std;

ROI::ROI(
    Image img, 
    Image mask,
    float score,
    int x1,
    int x2,
    int y1,
    int y2
): img(img), mask(mask){
    // img = img;
    // mask = mask;

    width = img.width();
    height = img.height();
    channels = img.channels();
    predictionScore = score;
    x1 = x1;
    x2 = x2;
    y1 = y1;
    y2 = y2;
};

void ROI::writeMask() {

}

void ROI::writeImg() {
	
}
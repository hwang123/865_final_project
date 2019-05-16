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
): img(img), mask(mask), x1(x1), x2(x2), y1(y1), y2(y2) {
    // img = img;
    // mask = mask;

    width = img.width();
    height = img.height();
    channels = img.channels();
    predictionScore = score;
};

void ROI::writeMask(string filename) {
	mask.write(filename);
}

void ROI::writeImg(string filename) {
	img.write(filename);
}

vector<int> ROI::getBoundingBox() {
	vector<int> bb{x1, x2, y1, y2};
	return bb;
}
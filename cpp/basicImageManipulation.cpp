/* --------------------------------------------------------------------------
 * File:    basicImageManipulation.cpp
 * Created: 2015-09-23
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/


#include "basicImageManipulation.h"
#include <assert.h>
using namespace std;
// FINAL PROJECT ---------------------------------------

// Generates two images of the same size of the mask input.
// The first image is a binary mask outlining the boundary of the mask.
// The second image is a 2-channel image denoting the normal direction to the boundary at each of
//  the boundary points. x-direction in the first channel, y-direction in the second.
vector<Image> boundaryTrace(Image &mask) {
    // First find starting pixel
    vector<Image> outputs;
    Image boundary = Image(mask.width(), mask.height(), mask.channels());
    Image slopes = Image(mask.width(), mask.height(), 2);
    
    vector<pair<int,int>> directions{
        make_pair(-1, 0),
        make_pair(-1, -1),
        make_pair(0, -1),
        make_pair(1, -1),
        make_pair(1, 0),
        make_pair(1, 1),
        make_pair(0, 1),
        make_pair(-1, 1)
    };

    int xo, yo;
    bool breakLoop = false;
    for (int y = 0; y < mask.height(); ++y) {
        if (breakLoop) break;
        for (int x = 0; x < mask.width(); ++x) {
            float pv = mask(x, y);
            if (pv > 0) {
                xo = x;
                yo = y;
                breakLoop = true;
                break;
            }
        }
    }

    boundary(xo, yo) = 1.0f;

    int directionIdx = 0;
    bool keepSearching = true;

    while (keepSearching) {
        for (int k = directionIdx + 1; k < directionIdx + 9; ++k) {
            int idx = k % 8;

            pair<int, int> direction = directions[idx];
            int newX = xo + direction.first;
            int newY = yo + direction.second;

            float boundaryVal = boundary.smartAccessor(newX, newY, 0);
            float maskVal = mask.smartAccessor(newX, newY, 0);

            if (newX == xo && newY == yo) {
                keepSearching = false;
                break;
            } else if (boundaryVal > 0) {
                keepSearching = false;
                break;
            } else if (maskVal > 0) {
                boundary(newX, newY) = 1.0f;

                int slopeIdx = (idx + 6) % 8;
                pair<int, int> slope = directions[slopeIdx];
                slopes(newX, newY, 0) = slope.first;
                slopes(newX, newY, 1) = slope.second;

                directionIdx = idx + 4;
                xo = newX;
                yo = newY;
                break;
            }
        }
    }

    outputs.push_back(boundary);
    outputs.push_back(slopes);
    return outputs;
}

Image padBoundary(vector<Image> &boundaries, int radius) {
    Image boundary = boundaries[0];
    Image directions = boundaries[1];
    Image padded = Image(boundary.width(), boundary.height(), boundary.channels());

    for (int x = 0; x < boundary.width(); ++x) {
        for (int y = 0; y < boundary.height(); ++y) {
            float boundaryVal = boundary(x, y);

            if (boundaryVal > 0) {
                int dx = directions(x, y, 0);
                int dy = directions(x, y, 1);
                
                for (int i = -3; i < 4; ++i) {
                    for (int j = -3; j < 4; ++j) {
                        int xIdx = x + i;
                        int yIdx = y + j;
                        for (int delta = 0; delta < radius; ++delta) {
                            int finalXIdx = xIdx + dx * delta;
                            int finalYIdx = yIdx + dy * delta;
                            if (finalXIdx < 0 || finalYIdx < 0 || finalXIdx > padded.width() - 1 || finalYIdx > padded.height() - 1) {
                                break;
                            }
                            float previousVal = padded(finalXIdx, finalYIdx);
                            float gaussian = exp(-1.0f * pow(delta, 2) / 200.0f);
                            // float gaussian = 1.0f - (1.0f / radius) * delta;

                            padded(finalXIdx, finalYIdx) = max(previousVal, gaussian);
                        }
                    }
                }
            }
        }
    }

    return padded;
}

Image combineMaskBoundary(Image &mask, Image &boundary) {
    Image combined = Image(boundary.width(), boundary.height(), boundary.channels());

    for (int x = 0; x < combined.width(); ++x) {
        for (int y = 0; y < combined.height(); ++y) {
            combined(x, y) = max(mask.smartAccessor(x, y,0), boundary(x, y));
        }
    }

    return combined;
}

Image overlayStyleTransferBoundary(Image& base, Image& styled, ROI& roi, int blurRadius) {
    Image overlayed = base;
    Image mask = roi.getMask();
    vector<Image> boundaries = boundaryTrace(mask);

    Image paddedBoundary = padBoundary(boundaries, floor(roi.getWidth()/8.0f));
    Image combinedBoundary = combineMaskBoundary(mask, paddedBoundary);
    combinedBoundary.write("./Output/mask.png");
    vector<int> bb = roi.getBoundingBox();
    int x1 = bb[0], x2 = bb[1], y1 = bb[2], y2 = bb[3];

    for (int i = 0; i < x2-x1-1; ++i) {
        for (int j = 0; j < y2-y1-1; ++j) {
            
            float alpha = combinedBoundary(i,j);
            // If the pixel is in the semantic mask
            if (alpha > 0) {
                int xBaseIdx = x1 + i;
                int yBaseIdx = y1 + j;

                for (int z = 0; z < base.channels(); ++z) {
                    float overlayedValue = alpha*styled(i, j, z) + (1.0f-alpha)*base(xBaseIdx, yBaseIdx, z);
                    overlayed(xBaseIdx, yBaseIdx, z) = overlayedValue;
                }
            } 
        }
    }

    return overlayed;
}

Image overlayStyleTransferNaive(Image& base, Image& styled, ROI& roi) {
    Image overlayed = base;

    vector<int> bb = roi.getBoundingBox();
    int x1 = bb[0], x2 = bb[1], y1 = bb[2], y2 = bb[3];

    for (int i = 0; i < x2-x1; ++i) {
        for (int j = 0; j < y2-y1; ++j) {
            int xBaseIdx = x1 + i;
            int yBaseIdx = y1 + j;
            for (int z = 0; z < 3; ++z) {
                overlayed(xBaseIdx, yBaseIdx, z) = styled(i, j, z);
            }
        }
    }

    return overlayed;
}

Image overlayStyleTransferAlpha(Image& base, Image& styled, ROI& roi, int radius) {
    Image overlayed = base;

    vector<int> bb = roi.getBoundingBox();
    int x1 = bb[0], x2 = bb[1], y1 = bb[2], y2 = bb[3];

    int roiWidth = x2-x1;
    int roiHeight = y2-y1;

    // First fully transfer the style to a tighter region
    for (int i = radius; i < roiWidth - radius; ++i) {
        for (int j = radius; j < roiHeight - radius; ++j) {
            int xBaseIdx = x1 + i;
            int yBaseIdx = y1 + j;
            for (int z = 0; z < 3; ++z) {
                overlayed(xBaseIdx, yBaseIdx, z) = styled(i, j, z);
            }
        }
    }

    // Then alpha-blend the outer band of the ROI
    // Can be thought of as 4 separate for-loops over each of the rectangular ROI edges
    int bandWidth = x2 - x1 - 2*radius;
    int bandHeight = y2 - y1 - 2*radius;

    float dx = 1.0f/radius;
    for (int i = 0; i < roiWidth; ++i) {
        for (int j = 0; j < radius; ++j) {
            int xIdx = x1 + i;
            int yIdx = y1 + j;

            float alpha = dx*j;
            for (int z = 0; z < base.channels(); ++z) {
                overlayed(xIdx, yIdx, z) = alpha*styled(i, j, z) + (1.0f-alpha)*base(xIdx, yIdx, z);
            }
        }
    }

    for (int i = 0; i < roiWidth; ++i) {
        for (int j = roiHeight - radius; j < roiHeight; ++j) {
            int xIdx = x1 + i;
            int yIdx = y1 + j;

            float alpha = 1.0f - dx*(j - roiHeight + radius);
            for (int z = 0; z < base.channels(); ++z) {
                overlayed(xIdx, yIdx, z) = alpha*styled(i, j, z) + (1.0f-alpha)*base(xIdx, yIdx, z);
            }
        }
    }

    for (int i = 0; i < radius; ++i) {
        for (int j = 0; j < roiHeight; ++j) {
            int xIdx = x1 + i;
            int yIdx = y1 + j;

            float alpha = dx * i;
            for (int z = 0; z < base.channels(); ++z) {
                overlayed(xIdx, yIdx, z) = alpha*styled(i, j, z) + (1.0f-alpha)*base(xIdx, yIdx, z);
            }
        }
    }

    for (int i = roiWidth - radius; i < roiWidth; ++i) {
        for (int j = 0; j < roiHeight; ++j) {
            int xIdx = x1 + i;
            int yIdx = y1 + j;

            float alpha = 1.0f - dx*(i - roiWidth + radius);
            for (int z = 0; z < base.channels(); ++z) {
                overlayed(xIdx, yIdx, z) = alpha*styled(i, j, z) + (1.0f-alpha)*base(xIdx, yIdx, z);
            }
        }
    }

    return overlayed;    
}

Image overlayStyleTransferMaskNaive(Image& base, Image& styled, ROI& roi) {
    Image overlayed = base;

    Image mask = roi.getMask();
    vector<int> bb = roi.getBoundingBox();
    int x1 = bb[0], x2 = bb[1], y1 = bb[2], y2 = bb[3];

    for (int i = 0; i < x2-x1; ++i) {
        for (int j = 0; j < y2-y1; ++j) {
            int maskVal = mask(i,j);

            // If the pixel is in the semantic mask
            if (maskVal > 0) {
                int xBaseIdx = x1 + i;
                int yBaseIdx = y1 + j;
                for (int z = 0; z < 3; ++z) {
                    overlayed(xBaseIdx, yBaseIdx, z) = styled(i, j, z);
                }
            } 
        }
    }

    return overlayed;
}


// -----------------------------------------------------

// --------- HANDOUT PS05 ------------------------------
// -----------------------------------------------------
//
Image scaleNN(const Image &im, float factor) {
    // --------- HANDOUT  PS05 ------------------------------
    // create a new image that is factor times bigger than the input by using
    // nearest neighbor interpolation.
    assert(factor > 0);
    if (factor == 1.0f) return im;

    float scalingFactor = factor > 1 ? factor : 1.0/factor;
    int scaledWidth = im.width() * scalingFactor;
    int scaledHeight = im.height() * scalingFactor;

    Image scaledImg = Image(scaledWidth, scaledHeight, im.channels());

    for (int x = 0; x < scaledWidth; ++x) {
        for (int y = 0; y < scaledHeight; ++y) {
            for (int z = 0; z < im.channels(); ++z) {
                int xIdx = round(x/scalingFactor);
                int yIdx = round(y/scalingFactor);
                
                scaledImg(x, y, z) = im.smartAccessor(xIdx, yIdx, z, true);
            }
        }
    }

    return scaledImg;
}

float interpolateLin(const Image &im, float x, float y, int z, bool clamp) {
    // --------- HANDOUT  PS05 ------------------------------
    // bilinear interpolation samples the value of a non-integral
    // position (x,y) from its four "on-grid" neighboring pixels.
    //  |           |
    // -1-----------2-
    //  |           |  *: my coordinates (x,y) are not integral
    //  |  *        |     since I am not on the pixel grid :(
    //  |           |  1: top-left
    //  |           |  2: top-right
    //  |           |  3: bottom-right
    // -4-----------3- 4: bottom-left, what are our coordinates?
    //  |           |    We are willing to share some color
    //                   information with * ! Of course, the pixel
    //                   closest to * should influence it more.
    int leftX = floor(x), rightX = ceil(x);
    int topY = floor(y), botY = ceil(y);

    float alphaXR = x - floor(x);
    float alphaXL = 1.0 - alphaXR;

    float alphaYB = y - floor(y);
    float alphaYT = 1.0 - alphaYB;

    float topLeftVal = im.smartAccessor(leftX, topY, z, clamp);
    float topRightVal = im.smartAccessor(rightX, topY, z, clamp);
    float botLeftVal = im.smartAccessor(leftX, botY, z, clamp);
    float botRightVal = im.smartAccessor(rightX, botY, z, clamp);

    float topInterpolatedVal = topLeftVal*alphaXL + topRightVal*alphaXR;
    float botInterpolatedVal = botLeftVal*alphaXL + botRightVal*alphaXR;

    return topInterpolatedVal*alphaYT + botInterpolatedVal*alphaYB;
}

Image scaleLin(const Image &im, float factor) {
    assert(factor > 0);
    if (factor == 1.0f) return im;

    float scalingFactor = factor > 1 ? factor : 1.0/factor;
    int scaledWidth = im.width() * scalingFactor;
    int scaledHeight = im.height() * scalingFactor;

    Image scaledImg = Image(scaledWidth, scaledHeight, im.channels());

    for (int x = 0; x < scaledWidth; ++x) {
        for (int y = 0; y < scaledHeight; ++y) {
            for (int z = 0; z < im.channels(); ++z) {
                float xIdx = x/scalingFactor;
                float yIdx = y/scalingFactor;
                
                scaledImg(x, y, z) = interpolateLin(im, xIdx, yIdx, z, true);
            }
        }
    }

    return scaledImg;
}

float getBicubicValue(float dx, float B, float C) {
    if (dx < 1) {
        return (12.0 - 9.0*B - 6.0*C)*pow(dx, 3.0) + (-18.0 + 12.0*B + 6.0*C)*pow(dx, 2) + (6.0 - 2.0*B);
    } else {
        return (-1.0*B - 6.0*C)*pow(dx, 3.0) + (6.0*B + 30.0*C)*pow(dx, 2.0) + (-12.0*B - 48.0*C)*dx + (8.0*B + 24.0*C);
    }
}

float interpolateBicubic(const Image &im, float x, float y, int z, float B, float C, bool clamp) {
    int startingX = floor(x) - 1, startingY = floor(y) - 1; 
    int endingX = ceil(x) + 1, endingY = ceil(y) + 1;

    float total = 0;
    for (int i = startingX; i < endingX + 1; ++i) {
        float dx = abs(x-i);
        float kx = (1.0/6.0) * getBicubicValue(dx, B, C);
        for (int j = startingY; j < endingY + 1; ++j) {
            float dy = abs(y-j);
            float ky = (1.0/6.0) * getBicubicValue(dy, B, C);

            float pixel = im.smartAccessor(i, j, z);
            total+= kx*ky*pixel;
        }
    }

    return total;
}

Image scaleBicubic(const Image &im, float factor, float B, float C) {
    // --------- HANDOUT  PS05 ------------------------------
    // create a new image that is factor times bigger than the input by using
    // a bicubic filter kernel with Mitchell and Netravali's parametrization
    // see "Reconstruction filters in computer graphics", Mitchell and Netravali 1988
    // or http://entropymine.com/imageworsener/bicubic/
    assert(factor > 0);
    if (factor == 1.0f) return im;

    float scalingFactor = factor > 1 ? factor : 1.0/factor;
    int scaledWidth = im.width() * scalingFactor;
    int scaledHeight = im.height() * scalingFactor;

    Image scaledImg = Image(scaledWidth, scaledHeight, im.channels());

    for (int x = 0; x < scaledWidth; ++x) {
        for (int y = 0; y < scaledHeight; ++y) {
            for (int z = 0; z < im.channels(); ++z) {
                float xIdx = x/scalingFactor;
                float yIdx = y/scalingFactor;
                
                scaledImg(x, y, z) = interpolateBicubic(im, xIdx, yIdx, z, B, C, true);
            }
        }
    }
    return scaledImg;
}

float evalSinc(float x) {
    if (x == 0.0) {
        return 1.0;
    }
    return sin(M_PI*x)/(M_PI*x);
}

float interpolateLanczos(const Image &im, float x, float y, int z, float a) {
    int startingX = floor(x) - (a-1), startingY = floor(y) - (a-1); 
    int endingX = ceil(x) + (a-1), endingY = ceil(y) + (a-1);

    float total = 0;
    for (int i = startingX; i < endingX + 1; ++i) {
        float dx = abs(x-i);
        float kx = evalSinc(dx)*evalSinc(dx/a);
        for (int j = startingY; j < endingY + 1; ++j) {
            float dy = abs(y-j);
            float ky = evalSinc(dy)*evalSinc(dy/a);

            float pixel = im.smartAccessor(i, j, z);
            total+= kx*ky*pixel;
        }
    }

    return total;
}

Image scaleLanczos(const Image &im, float factor, float a) {
    // --------- HANDOUT  PS05 ------------------------------
    // create a new image that is factor times bigger than the input by using
    // a Lanczos filter kernel
    assert(factor > 0);
    if (factor == 1.0f) return im;

    float scalingFactor = factor > 1 ? factor : 1.0/factor;
    int scaledWidth = im.width() * scalingFactor;
    int scaledHeight = im.height() * scalingFactor;

    Image scaledImg = Image(scaledWidth, scaledHeight, im.channels());

    for (int x = 0; x < scaledWidth; ++x) {
        for (int y = 0; y < scaledHeight; ++y) {
            for (int z = 0; z < im.channels(); ++z) {
                float xIdx = x/scalingFactor;
                float yIdx = y/scalingFactor;
                
                scaledImg(x, y, z) = interpolateLanczos(im, xIdx, yIdx, z, a);
            }
        }
    }
    return scaledImg;
}

Image rotate(const Image &im, float theta) {
    // --------- HANDOUT  PS05 ------------------------------
    // rotate an image around its center by theta

	// // center around which to rotate
    // float centerX = (im.width()-1.0)/2.0;
    // float centerY = (im.height()-1.0)/2.0;

    return im; // changeme

}

// -----------------------------------------------------
// --------- END --- PS05 ------------------------------


// --------- HANDOUT PS01 ------------------------------
// -----------------------------------------------------

// Change the brightness of the image
// const Image & means a reference to im will get passed to the function,
// but the compiler won't let you modify it within the function.
// So you will return a new image
Image brightness(const Image &im, float factor) {
    // // --------- HANDOUT  PS01 ------------------------------
	// // Image output(im.width(), im.height(), im.channels());
	// // Modify image brightness
	// // return output;
	// return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    return im * factor;
}

Image contrast(const Image &im, float factor, float midpoint) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Image output(im.width(), im.height(), im.channels());
    // // Modify image contrast
    // // return output;
	// return Image(1,1,1); //Change this

    // --------- SOLUTION PS01 ------------------------------
    return (im - midpoint) * factor + midpoint;
}

Image color2gray(const Image &im, const std::vector<float> &weights) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Image output(im.width(), im.height(), 1);
    // // Convert to grayscale
	// return Image(1,1,1); //Change this

    // --------- SOLUTION PS01 ------------------------------
    Image output(im.width(), im.height(), 1);
    for (int i = 0 ; i < im.width(); i++ ) {
        for (int j = 0 ; j < im.height(); j++ ) {
            output(i,j,0) = im(i,j,0) * weights[0] + im(i,j,1) * weights[1] + im(i,j,2) *weights[2];
        }
    }
    return output;
}

// For this function, we want two outputs, a single channel luminance image
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<Image> lumiChromi(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create the luminance image
    // // Create the chrominance image
    // // Create the output vector as (luminance, chrominance)
	// return std::vector<Image>(); //Change this

    // --------- SOLUTION PS01 ------------------------------

    // Create the luminance
    Image im_luminance = color2gray(im);

    // Create chrominance images
    // We copy the input as starting point for the chrominance
    Image im_chrominance = im;
    for (int c = 0 ; c < im.channels(); c++ ) {
        for (int y = 0 ; y < im.height(); y++) {
            for (int x = 0 ; x < im.width(); x++) {
                im_chrominance(x,y,c) = im_chrominance(x,y,c) / im_luminance(x,y);
            }
        }
    }

    // Stack luminance and chrominance in the output vector, luminance first
    return std::vector<Image>{im_luminance, im_chrominance};
}

Image lumiChromi2rgb(const vector<Image> & lc) {
    // luminance is lc[0]
    // chrominance is lc[1]

    // Create chrominance images
    // We copy the input as starting point for the chrominance
    Image im = Image(lc[1].width(), lc[1].height(), lc[1].channels());
    for (int c = 0 ; c < im.channels(); c++ ) {
      for (int y = 0 ; y < im.height(); y++) {
        for (int x = 0 ; x < im.width(); x++) {
            im(x,y,c) = lc[1](x,y,c) * lc[0](x,y);
        }
      }
    }
    return im;
}


// Modify brightness then contrast
Image brightnessContrastLumi(const Image &im, float brightF, float contrastF, float midpoint) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Modify brightness, then contrast of luminance image
    // return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    // Separate luminance and chrominance
    std::vector<Image> lumi_chromi = lumiChromi(im);
    Image im_luminance             = lumi_chromi[0];
    Image im_chrominance           = lumi_chromi[1];

    // Process the luminance channel
    im_luminance = brightness(im_luminance, brightF);
    im_luminance = contrast(im_luminance, contrastF, midpoint);

    // Multiply the chrominance with the new luminance to get the final image
    for (int i = 0 ; i < im.width(); i++ ){
        for (int j = 0 ; j < im.height(); j++) {
            for (int c = 0; c < im.channels(); c++) {
                im_chrominance(i,j,c) = im_chrominance(i,j,c) * im_luminance(i,j);
            }
        }
    }
    // At this point, im_chrominance olds the complete processed image
    return im_chrominance;
}


Image rgb2yuv(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create output image of appropriate size
    // // Change colorspace
    // return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    Image output(im.width(), im.height(), im.channels());
    for (int j = 0 ; j < im.height(); j++) {
        for (int i = 0 ; i < im.width(); i++) {
            output(i,j,0) =   0.299 * im(i,j,0) + 0.587 * im(i,j,1) + 0.114 * im(i,j,2);
            output(i,j,1) = - 0.147 * im(i,j,0) - 0.289 * im(i,j,1) + 0.436 * im(i,j,2);
            output(i,j,2) =   0.615 * im(i,j,0) - 0.515 * im(i,j,1) - 0.100 * im(i,j,2);
        }
    }
    return output;
}


Image yuv2rgb(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create output image of appropriate size
    // // Change colorspace
    // return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    Image output(im.width(), im.height(), im.channels());
    for (int j = 0 ; j < im.height(); j++) {
        for (int i = 0; i < im.width(); i++)
        {
            output(i,j,0) =  im(i,j,0) + 0     * im(i,j,1) + 1.14  * im(i,j,2);
            output(i,j,1) =  im(i,j,0) - 0.395 * im(i,j,1) - 0.581 * im(i,j,2);
            output(i,j,2) =  im(i,j,0) + 2.032 * im(i,j,1) + 0     * im(i,j,2);
        }
    }
    return output;
}


Image saturate(const Image &im, float factor) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create output image of appropriate size
    // // Saturate image
    // // return output;
    // return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    Image output = rgb2yuv(im); // Change colorspace
    for (int i = 0 ; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            output(i,j,1) = output(i,j,1) * factor;
            output(i,j,2) = output(i,j,2) * factor;
        }
    }
    output = yuv2rgb(output); // Back to RGB
    return output;
}


// Return two images in a C++ vector
std::vector<Image> spanish(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Remember to create the output images and the output vector
    // // Push the images onto the vector
    // // Do all the required processing
    // // Return the vector, color image first
	// return std::vector<Image>(); //Change this

    // --------- SOLUTION PS01 ------------------------------
    // Extract the luminance
    Image output_L = color2gray(im);

    // Convert to YUV for manipulation
    Image output_C = rgb2yuv(im);

    for (int j = 0; j < im.height(); j++) {
        for (int i = 0; i < im.width(); i++) {
            output_C(i,j,0) = 0.5;              // constant luminance
            output_C(i,j,1) = -output_C(i,j,1); // opposite chrominance
            output_C(i,j,2) = -output_C(i,j,2); // opposite chrominance
        }
    }
    // Convert back to RGB
    output_C = yuv2rgb(output_C);

    // Location of the black dot
    int bdot_x = floor(im.width()/2);
    int bdot_y = floor(im.height()/2);

    // Add the black dot to Luminance, and Chrominance images
    output_L(bdot_x, bdot_y,0) = 0.0f;
    output_C(bdot_x, bdot_y,0) = 0.0f; // black is 0
    output_C(bdot_x, bdot_y,1) = 0.0f;
    output_C(bdot_x, bdot_y,2) = 0.0f;

    // Pack the images in a vector, chrominance first
    return std::vector<Image>{output_C, output_L} ;
}


// White balances an image using the gray world assumption
Image grayworld(const Image & im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // Implement automatic white balance by multiplying each channel
    // of the input by a factor such that the three channel of the output image
    // have the same mean value. The mean value of the green channel
    // is taken as reference.
    // return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    // Compute the mean per channel
    float mean_r = 0, mean_g = 0, mean_b = 0;
    float N = im.width()*im.height();
    for (int j = 0 ; j < im.height(); j++) {
        for (int i = 0 ; i < im.width(); i++) {
            mean_r += im(i,j,0);
            mean_g += im(i,j,1);
            mean_b += im(i,j,2);
        }
    }
    mean_r /= N;
    mean_g /= N;
    mean_b /= N;

    Image output = im;
    for (int j = 0 ; j < im.height();j ++) {
        for (int i = 0 ; i < im.width(); i++) {
            output(i,j,0) = output(i,j,0)/mean_r*mean_g;
            // dont process output(i,j,1), since the mean of
            // the green channel is already at the right value
            output(i,j,2) = output(i,j,2)/mean_b*mean_g;
        }
    }
    return output;
}


Image gamma_code(const Image &im, float gamma) {
    // // --------- HANDOUT  PS01 ------------------------------
    // Image output(im.width(), im.height(), im.channels());
    // Gamma encodes the image
    // return output;

    // --------- SOLUTION PS01 ------------------------------
    Image output = Image(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.number_of_elements(); ++i){
        output(i) = pow(im(i), (1/gamma));
    }
    return output;
}

// -----------------------------------------------------
// --------- END --- PS01 ------------------------------

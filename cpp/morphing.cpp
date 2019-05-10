/* -----------------------------------------------------------------
 * File:    morphing.cpp
 * Created: 2015-09-25
 * -----------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------*/




#include <cassert>
#include "morphing.h"

using namespace std;

Vec2f operator+(const Vec2f & a, const Vec2f & b) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return the vector sum of a an b
    return Vec2f(a.x + b.x, a.y + b.y); // change me
}


Vec2f operator-(const Vec2f & a, const Vec2f & b) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return a-b
    return Vec2f(a.x - b.x, a.y - b.y); // change me
}


Vec2f operator*(const Vec2f & a, float f) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return a*f
    return Vec2f(f * a.x, f * a.y); // change me
}

Vec2f operator/(const Vec2f & a, float f) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return a/f
    return Vec2f(a.x / f, a.y / f); // change me
}

float dot(const Vec2f & a, const Vec2f & b) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return the dot product of a and b.
    return a.x*b.x + a.y*b.y; // change me
}

float length(const Vec2f & a) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return the length of a.
    return pow(pow(a.x, 2.0) + pow(a.y, 2.0), 1.0/2.0); // change me
}

Vec2f perpendicular(const Vec2f & a) {
    // --------- HANDOUT  PS05 ------------------------------
    // Return a vector that is perpendicular to a.
    // Either direction is fine.
    return Vec2f(-a.y, a.x);
}

// The Segment constructor takes in 2 points P(x1,y1) and Q(x2,y2) corresponding to
// the ends of a segment and initialize the local reference frame e1,e2.
Segment::Segment(Vec2f P_, Vec2f Q_) : P(P_), Q(Q_) {
    // // --------- HANDOUT  PS05 ------------------------------
    // // The initializer list above ": P(P_), Q(Q_)" already copies P_
    // // and Q_, so you don't have to do it in the body of the constructor.
    // You should:
    // * Initialize the local frame e1,e2 (see header file)
    lPQ = length(P-Q);
    e1 = (P-Q) / length(P-Q);
    Vec2f ortho = perpendicular(e1);
    e2 = ortho / length(ortho);
}


Vec2f Segment::XtoUV(Vec2f X) const {
    // --------- HANDOUT  PS05 ------------------------------
    // Compute the u,v coordinates of a point X with
    // respect to the local frame of the segment.
    // e2 ^
    //    |
    // v  +  * X
    //    | /
    //    |/
    //    *--+------>-----*
    //    P  u     e1     Q
    //                    u=1
    //
    // * Be careful with the different normalization for u and v
    float U = dot((X-P), (Q-P));
    U = U / pow(length(Q-P), 2.0);

    float V = dot((X-P), perpendicular(Q-P));
    V = V / length(Q-P);

    return Vec2f(U, V); 
}


Vec2f Segment::UVtoX(Vec2f uv) const {
    // --------- HANDOUT  PS05 ------------------------------
    // compute the (x, y) position of a point given by the (u,v)
    // location relative to this segment.
    // * Be careful with the different normalization for u and v
    float u = uv.x;
    float v = uv.y;
    Vec2f left = P + ((Q - P) * u);
    Vec2f right = (perpendicular(Q-P) * v) / length(Q-P);
    return left + right;
}


float Segment::distance(Vec2f X) const {
    // --------- HANDOUT  PS05 ------------------------------
    // Implement distance from a point X(x,y) to the segment. Remember the 3
    // cases from class.
    Vec2f UV = XtoUV(X);
    if (UV.x < 0) {
        return length(X-P);
    } else if (UV.x > 1) {
        return length(X-Q);
    }
    return abs(UV.y);
}


Image warpBy1(const Image &im, const Segment &segBefore, const Segment &segAfter){
    // --------- HANDOUT  PS05 ------------------------------
    // Warp an entire image according to a pair of segments.
    Image warped = Image(im.width(), im.height(), im.channels());
    
    for (int i = 0; i < warped.width(); ++i) {
        for (int j = 0; j < warped.height(); ++j) {
            Vec2f X(i, j);
            Vec2f UV = segAfter.XtoUV(X);
            Vec2f X_src = segBefore.UVtoX(UV);

            for (int z = 0; z < warped.channels(); ++z) {
                warped(i, j, z) = interpolateLin(im, X_src.x, X_src.y, z, true);
            }
        }
    }
    return warped;
}


float Segment::weight(Vec2f X, float a, float b, float p) const {
    // --------- HANDOUT  PS05 ------------------------------
    // compute the weight of a segment to a point X(x,y) given the weight
    // parameters a,b, and p (see paper for details).
    float numerator = pow(length(X), p);
    float denominator = a + distance(X);
    return pow(numerator/denominator, b);
}


Image warp(const Image &im, const vector<Segment> &src_segs,
        const vector<Segment> &dst_segs, float a, float b, float p) {

    Image warped = Image(im.width(), im.height(), im.channels());

    for (int i = 0; i < warped.width(); ++i) {
        for (int j = 0; j < warped.height(); ++j) {
            Vec2f X(i, j);
            for (int z = 0; z < warped.channels(); ++z) {
                float pixelValue = 0;
                float total_weight = 0.0;
                Vec2f finalIdx(0, 0);
                for (int k = 0; k < src_segs.size(); ++k) {
                    Segment segBefore = src_segs[k];
                    Segment segAfter = dst_segs[k];

                    Vec2f UV = segAfter.XtoUV(X);
                    Vec2f X_src = segBefore.UVtoX(UV);

                    float segWeight = segAfter.weight(X, a, b, p);
                    total_weight += segWeight;
                    finalIdx = finalIdx + X_src*segWeight;
                }

                finalIdx = finalIdx / total_weight;
                warped(i, j, z) = interpolateLin(im, finalIdx.x, finalIdx.y, z, true);
            }
        }
    }
    return warped;
}


vector<Image> morph(const Image &im_before, const Image &im_after,
        const vector<Segment> &segs_before, const vector<Segment> &segs_after,
        int N, float a, float b, float p) {
    // --------- HANDOUT  PS05 ------------------------------
    // return a vector of N+2 images: the two inputs plus N images that morphs
    // between im_before and im_after for the corresponding segments. im_before should be the first image, im_after the last.
    float dt = 1.0/(N+1.0);
    vector<Image> morphed_imgs;
    morphed_imgs.push_back(im_before);
    for (float t = dt; t < 1.0; t += dt){
        vector<Segment> segs_middle;

        for (int i = 0; i < segs_before.size(); ++i) {
            Segment segBefore = segs_before[i];
            Segment segAfter = segs_after[i];

            Vec2f P_before = segBefore.getP(), Q_before = segBefore.getQ();
            Vec2f P_after = segAfter.getP(), Q_after = segAfter.getQ();

            Vec2f diffP = P_after - P_before;
            Vec2f diffQ = Q_after - Q_before;

            Vec2f newP = P_before + (diffP*t);
            Vec2f newQ = Q_before + (diffQ*t);
            Segment interpolated_seg(newP, newQ);
            // cout << "-----------------" << endl;
            // cout << "P before: " << P_before.x << "," << P_before.y << endl;
            // cout << "Q before: " << Q_before.x << "," << Q_before.y << endl;

            // cout << "P after: " << P_after.x << "," << P_after.y << endl;
            // cout << "Q after: " << Q_after.x << "," << Q_after.y << endl;       

            // cout << "P middle: " << newP.x << "," << newP.y << endl;
            // cout << "Q middle: " << newQ.x << "," << newQ.y << endl;            

            segs_middle.push_back(interpolated_seg);
        }

        Image warped_im_before = warp(im_before, segs_before, segs_middle, a, b, p);
        Image warped_im_after = warp(im_after, segs_after, segs_middle, a, b, p);

        Image interpolated_warp = warped_im_before*(1-t) + warped_im_after*t;
        // Image interpolated_warp = warped_im_before;

        morphed_imgs.push_back(interpolated_warp);
    }

    morphed_imgs.push_back(im_after);
    return morphed_imgs;
}

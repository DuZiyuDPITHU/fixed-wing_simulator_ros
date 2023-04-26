#ifndef _UNIFORM_BSPLINE_H_
#define _UNIFORM_BSPLINE_H_
#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>
#include "uniform_bspline.h"

using namespace std;

class BsplineOpt
{
private:
    double lambda1_;               // smoothness weight
    double lambda2_;               // distance weight
    double lambda3_;               // feasibility weight
    double lambda4_;               // curvature weight
    double lambda5_;               // fitting weight

    UniformBspline* bspline_ptr;
public:
    BsplineOpt() {}
    ~BsplineOpt() {}

};

#endif
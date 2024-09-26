#ifndef _FITGLE_H_
#define _FITGLE_H_

#include <cstdlib>
#include <cstdio>
#include <vector>
#include "Frame.h"
#include <memory>
#include <gsl/gsl_bspline.h>

namespace FITGLE_NS {

struct InputParameters
{
    double start;  //start distance r0
    double end;    //end distance r1
    double start_density;  //start density d0
    double end_density;    //end density d1
    int    splineOrder;
    int    numSplines;
    int    numSplinesDensity;
    double outputPrecision;
    double boxlox; 
    double boxloy; 
    double boxloz; 
    double boxhix; 
    double boxhiy; 
    double boxhiz; 
    int    steps;
    int    maxIter;
    FILE*  fileTraj;
};  //Structure to store input parameters

class FitGLE
{
public:
    FitGLE(int argc, char** argv);
    ~FitGLE();
    void exec();

    //helper functions
    void accumulatePairNormalEquation();
    void accumulateDensityNormalEquation();
    void leastSquareSolver(bool isPair);
    void output(int iter);

    //helper functions
    void clearNormalMatrix();
    double distance(std::vector<double> &, std::vector<double> &);
    double localDensity(int i);
    std::vector<double> parallelVelocity(int i, int j);
    std::vector<double> parallelUnitVector(int i, int j);

private:
    std::shared_ptr<class Frame> trajFrame;
    std::shared_ptr<struct InputParameters> info; 
    std::vector<double> divPoints;  // divide points of the b-spline radial ranges
    
    int totalSplineCount;
    std::vector<std::vector<double> > normalMatrix;
    std::vector<double> normalVector;
    std::vector<double> splineCoefficientsPair;
    std::vector<double> splineCoefficientsDensity;

    //gsl members for pairwise force b-splines
    gsl_bspline_workspace *bwPair;
    gsl_vector *splineValuePair;

    //gsl members for density b-splines
    gsl_bspline_workspace *bwDensity;
    gsl_vector *splineValueDensity;

};

}

#endif

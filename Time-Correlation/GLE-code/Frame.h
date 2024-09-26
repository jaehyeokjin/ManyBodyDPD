#ifndef _FRAME_H_
#define _FRAME_H_

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "FitGLE.h"

namespace FITGLE_NS {

class Frame
{
public:
    friend class FitGLE;
    Frame(int n, double lx, double ly, double lz, double hx, double hy, double hz, char* filename);
    ~Frame();
    void readFrame();
    void jumpBackHead();
    int get();
    double distance(int i, int j);
    void genLocalDensity();
private:
    FILE* trajectory;
    int numParticles;
    double boxlox;
    double boxloy;
    double boxloz;
    double boxhix;
    double boxhiy;
    double boxhiz;
    std::vector<std::vector<double> > positions;
    std::vector<std::vector<double> > residualForces;
    std::vector<std::vector<double> > velocities;
    std::vector<double> density;
};
}

#endif

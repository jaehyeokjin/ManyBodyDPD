#include "Frame.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
using namespace FITGLE_NS;

Frame::Frame(int n, double lx, double ly, double lz, double hx, double hy, double hz, char* fileName)
{
    trajectory = fopen(fileName, "r");
    numParticles = n;
    boxlox = lx;
    boxloy = ly;
    boxloz = lz;
    boxhix = hx;
    boxhiy = hy;
    boxhiz = hz;
    positions.resize(numParticles);
    residualForces.resize(numParticles);
    velocities.resize(numParticles);
    density.resize(numParticles);
    for (int i=0; i<numParticles; i++)
    {
      positions[i].resize(3);
      residualForces[i].resize(3);
      velocities[i].resize(3);
    }
    printf("finishing initializing trajectory frames, %d Particles\n", numParticles);
}

Frame::~Frame()
{
    fclose(trajectory);
    printf("cleaning up Frame Information\n");
}

int Frame::get()
{
  return numParticles;
}

void Frame::readFrame()
{
    ssize_t read;
    size_t len;
    char*  line = NULL;
    int    lineID;
    for (int iline = 0; iline < numParticles+9; iline++)
    {
        read = getline(&line, &len, trajectory);
        if (iline >= 9)
        {
           char* pch = strtok(line, " \t");
           int atomID = atoi(pch) - 1;
           pch = strtok(NULL, " \t");
           int type = atoi(pch);
           pch = strtok(NULL, " \t");
           double atomMass = atof(pch);
           pch = strtok(NULL, " \t");
           positions[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           positions[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           positions[atomID][2] = atof(pch);
           pch = strtok(NULL, " \t");
           velocities[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           velocities[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           velocities[atomID][2] = atof(pch);
           pch = strtok(NULL, " \t");
           residualForces[atomID][0] = atof(pch);
           pch = strtok(NULL, " \t");
           residualForces[atomID][1] = atof(pch);
           pch = strtok(NULL, " \t");
           residualForces[atomID][2] = atof(pch);
        }
    }
    genLocalDensity();
}

void Frame::jumpBackHead()
{
    // Move the file pointer to the start.
    fseek(trajectory, 0, SEEK_SET);
}

double Frame::distance(int i, int j)
{
    double dx = positions[i][0] - positions[j][0];
    double boxLengthx = boxhix - boxlox;
    if (dx > 0.5 * boxLengthx) dx = dx - boxLengthx;
    if (dx < -0.5 * boxLengthx) dx = dx + boxLengthx;

    double dy = positions[i][1] - positions[j][1];
    double boxLengthy = boxhiy - boxloy;
    if (dy > 0.5 * boxLengthy) dy = dy - boxLengthy;
    if (dy < -0.5 * boxLengthy) dy = dy + boxLengthy;

    double dz = positions[i][2] - positions[j][2];
    double boxLengthz = boxhiz - boxloz;
    if (dz > 0.5 * boxLengthz) dz = dz - boxLengthz;
    if (dz < -0.5 * boxLengthz) dz = dz + boxLengthz;
  
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void Frame::genLocalDensity()
{
    double h = 6.0;
    for (int i=0; i<numParticles; i++) density[i] = 0.0;

    for (int i=0; i < numParticles - 1; i++)
    {
        for (int j=i+1; j < numParticles; j++)
        {
            double rij = distance(i, j);
            //double weight = 0.5 - 0.5 * tanh((rij - h) / (0.1*h) );  // Change the function form here : YINING NOTE : Need to improve
            double weight = 0.0;
	    if (rij < h) 
            {
                weight = (1.0 + 3.0 * rij / h) * (1.0 - rij / h) * (1.0 - rij / h) * (1.0 - rij / h);
                //weight = (1.0 - rij / h) * (1.0 - rij / h);
            }
	    density[i] += weight;
	    density[j] += weight;
        }
    }
}

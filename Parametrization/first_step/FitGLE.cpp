#include "FitGLE.h"
#include "comm.h"
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_bspline.h>
#include <lapacke.h>

using namespace FITGLE_NS;

FitGLE::FitGLE(int argc, char** argv)
{
    if (argc != 3)
    {
        printf("./FitGLE.x [trajectory Filename] [number of Particles]\n");
    }

    assert(argc == 3);
    printf("Initializing FitGLE parameters...\n");

    // parsing the configuration parameters
    info = std::make_shared<InputParameters>();
    VAR_BEGIN
      GET_REAL(info->start)
      GET_REAL(info->end)
      GET_REAL(info->start_density)
      GET_REAL(info->end_density)
      GET_REAL(info->boxlox)
      GET_REAL(info->boxloy)
      GET_REAL(info->boxloz)
      GET_REAL(info->boxhix)
      GET_REAL(info->boxhiy)
      GET_REAL(info->boxhiz)
      GET_REAL(info->outputPrecision)
      GET_INT(info->splineOrder)
      GET_INT(info->numSplines)
      GET_INT(info->numSplinesDensity)
      GET_INT(info->steps)
      GET_INT(info->maxIter)
    VAR_END
           
    printf("set up trajectory files\n");
    trajFrame = std::make_shared<Frame>(atoi(argv[2]), info->boxlox, info->boxloy, info->boxloz, info->boxhix, info->boxhiy, info->boxhiz, argv[1]);

    // Initialize the Normal Equation matrix and vectors
    // Set up the size of splines according to order and numbers
    printf("set up b-spline data structures\n");
    int numBreaks = info->numSplines + 2 - info->splineOrder;
    int numBreaksDensity = info->numSplinesDensity + 2 - info->splineOrder;

    //totalSplineCount = info->numSplines + info->numSplinesDensity;
    totalSplineCount = 1;
    normalVector.resize(totalSplineCount);

    splineCoefficientsPair.resize(info->numSplines);
    splineCoefficientsDensity.resize(info->numSplinesDensity);

    normalMatrix.resize(totalSplineCount);
    printf("set up containers\n");
    for (auto&& i : normalMatrix)
    {
        i.resize(totalSplineCount);
    }

    // Initialize the spline set up
    bwPair = gsl_bspline_alloc(info->splineOrder, numBreaks);
    splineValuePair = gsl_vector_alloc(info->numSplines);
    gsl_bspline_knots_uniform(info->start, info->end, bwPair);

    bwDensity = gsl_bspline_alloc(info->splineOrder, numBreaksDensity);
    splineValueDensity = gsl_vector_alloc(info->numSplinesDensity);
    gsl_bspline_knots_uniform(info->start_density, info->end_density, bwDensity);

    // Initially assume a uniform density dependent form
    //for (int m=0; m<totalSplineCount; m++)
    //{
    //    splineCoefficientsDensity[m] = 1.0;
    //}
    printf("finishing configuration, entering normal equation accumulation\n");
}

FitGLE::~FitGLE()
{
    gsl_bspline_free(bwPair);
    gsl_vector_free(splineValuePair);
    gsl_bspline_free(bwDensity);
    gsl_vector_free(splineValueDensity);
    printf("Exiting the Regression process...\n");
}

inline double FitGLE::distance(std::vector<double> & A, std::vector<double> & B)
{
    double dx = A[0] - B[0];
    double boxLengthx = info->boxhix - info->boxlox;
    if (dx > 0.5 * boxLengthx) dx = dx - boxLengthx;
    if (dx < -0.5 * boxLengthx) dx = dx + boxLengthx;

    double dy = A[1] - B[1];
    double boxLengthy = info->boxhiy - info->boxloy;
    if (dy > 0.5 * boxLengthy) dy = dy - boxLengthy;
    if (dy < -0.5 * boxLengthy) dy = dy + boxLengthy;

    double dz = A[2] - B[2];
    double boxLengthz = info->boxhiz - info->boxloz;
    if (dz > 0.5 * boxLengthz) dz = dz - boxLengthz;
    if (dz < -0.5 * boxLengthz) dz = dz + boxLengthz;
  
    return sqrt(dx*dx + dy*dy + dz*dz);
}

inline std::vector<double> FitGLE::parallelVelocity(int i, int j)
{
    double dx = trajFrame->positions[i][0] - trajFrame->positions[j][0];
    double boxLengthx = info->boxhix - info->boxlox;
    if (dx > 0.5 * boxLengthx) dx = dx - boxLengthx;
    if (dx < -0.5 * boxLengthx) dx = dx + boxLengthx;

    double dy = trajFrame->positions[i][1] - trajFrame->positions[j][1];
    double boxLengthy = info->boxhiy - info->boxloy;
    if (dy > 0.5 * boxLengthy) dy = dy - boxLengthy;
    if (dy < -0.5 * boxLengthy) dy = dy + boxLengthy;

    double dz = trajFrame->positions[i][2] - trajFrame->positions[j][2];
    double boxLengthz = info->boxhiz - info->boxloz;
    if (dz > 0.5 * boxLengthz) dz = dz - boxLengthz;
    if (dz < -0.5 * boxLengthz) dz = dz + boxLengthz;

    double rij = sqrt(dx*dx + dy*dy + dz*dz);
    double eij[] = {dx/rij, dy/rij, dz/rij};
    std::vector<double> vij;
    std::transform(trajFrame->velocities[i].begin(), trajFrame->velocities[i].end(), trajFrame->velocities[j].begin(), std::back_inserter(vij), std::minus<double>());
     
    double projection = vij[0] * eij[0] + vij[1] * eij[1] + vij[2] * eij[2];
    vij[0] = projection * eij[0];
    vij[1] = projection * eij[1];
    vij[2] = projection * eij[2];
    return vij;
}

inline std::vector<double> FitGLE::parallelUnitVector(int i, int j)
{
    double dx = trajFrame->positions[i][0] - trajFrame->positions[j][0];
    double boxLengthx = info->boxhix - info->boxlox;
    if (dx > 0.5 * boxLengthx) dx = dx - boxLengthx;
    if (dx < -0.5 * boxLengthx) dx = dx + boxLengthx;

    double dy = trajFrame->positions[i][1] - trajFrame->positions[j][1];
    double boxLengthy = info->boxhiy - info->boxloy;
    if (dy > 0.5 * boxLengthy) dy = dy - boxLengthy;
    if (dy < -0.5 * boxLengthy) dy = dy + boxLengthy;

    double dz = trajFrame->positions[i][2] - trajFrame->positions[j][2];
    double boxLengthz = info->boxhiz - info->boxloz;
    if (dz > 0.5 * boxLengthz) dz = dz - boxLengthz;
    if (dz < -0.5 * boxLengthz) dz = dz + boxLengthz;

    double rij = sqrt(dx*dx + dy*dy + dz*dz);
    std::vector<double> eij;
    eij.push_back(dx / rij);
    eij.push_back(dy / rij);
    eij.push_back(dz / rij);
    //printf("%lf %lf %lf d %lf %lf %lf\n", dx, dy, dz, trajFrame->positions[i][0], trajFrame->positions[j][0], info->boxLength);
    return eij;
}

void FitGLE::clearNormalMatrix()
{
    for (int i=0; i<totalSplineCount; i++)
    {
        for (int j=0; j<totalSplineCount; j++)
        {
            normalMatrix[i][j] = 0.0;
        }
    }
}

// Accumulate the normal equation for pairwise interaction
void FitGLE::accumulatePairNormalEquation()
{
    int nall = trajFrame->numParticles;
    int nSplines = totalSplineCount;
    std::vector<std::vector<double> > frameMatrix(nSplines, std::vector<double>(3*nall));
    double normalFactor = 1.0 / info->steps;
   
    // Computing Matrix F_km 
    for (int i=0; i<nall; i++)
    {
        printf("[P%d] rho = %lf\n", i+1, trajFrame->density[i]);
        for (int j = i+1; j<nall; j++)
        {
            double rij = distance(trajFrame->positions[i], trajFrame->positions[j]);
	    double rhoi = trajFrame->density[i];
	    double rhoj = trajFrame->density[j];

            std::vector<double> dv = parallelUnitVector(i, j);
            double dwdr = 0.0;
            double h = 6.0;
	    if (rij < h)
            {
                //weight = (1.0 - rij / h) * (1.0 - rij / h);
                dwdr = -12.0 * (1.0 - rij / h) * (1.0 - rij / h) * rij / h;
            }

            frameMatrix[0][3*i]     += -dwdr * (rhoi + rhoj) * dv[0];
            frameMatrix[0][3*i + 1] += -dwdr * (rhoi + rhoj) * dv[1];
            frameMatrix[0][3*i + 2] += -dwdr * (rhoi + rhoj) * dv[2];
            frameMatrix[0][3*j]     -= -dwdr * (rhoi + rhoj) * dv[0];
            frameMatrix[0][3*j + 1] -= -dwdr * (rhoi + rhoj) * dv[1];
            frameMatrix[0][3*j + 2] -= -dwdr * (rhoi + rhoj) * dv[2];
            
	    // Calculate the (F(rho_i)+F(rho_j)) term
	    //printf("rhoi = %lf rhoj = %lf\n", rhoi, rhoj);
	    //gsl_bspline_eval(rhoi, splineValueDensity, bwDensity);
            //for (int m=0; m<info->numSplinesDensity; m++)
            //{
            //   double phim = gsl_vector_get(splineValueDensity, m);
            //   frameMatrix[m + info->numSplines][3*i]     += -phim * dv[0] * dwdr;
            //   frameMatrix[m + info->numSplines][3*i + 1] += -phim * dv[1] * dwdr;
            //   frameMatrix[m + info->numSplines][3*i + 2] += -phim * dv[2] * dwdr;
            //   frameMatrix[m + info->numSplines][3*j]     -= -phim * dv[0] * dwdr;
            //   frameMatrix[m + info->numSplines][3*j + 1] -= -phim * dv[1] * dwdr;
            //   frameMatrix[m + info->numSplines][3*j + 2] -= -phim * dv[2] * dwdr;
            //}

	    //gsl_bspline_eval(rhoj, splineValueDensity, bwDensity);
            ////printf("here\n");
            //for (int m=0; m<info->numSplinesDensity; m++)
            //{
            //   double phim = gsl_vector_get(splineValueDensity, m);
            //   frameMatrix[m + info->numSplines][3*i]     += -phim * dv[0] * dwdr;
            //   frameMatrix[m + info->numSplines][3*i + 1] += -phim * dv[1] * dwdr;
            //   frameMatrix[m + info->numSplines][3*i + 2] += -phim * dv[2] * dwdr;
            //   frameMatrix[m + info->numSplines][3*j]     -= -phim * dv[0] * dwdr;
            //   frameMatrix[m + info->numSplines][3*j + 1] -= -phim * dv[1] * dwdr;
            //   frameMatrix[m + info->numSplines][3*j + 2] -= -phim * dv[2] * dwdr;
            //}

	    //if (rij < info->end)
            //{
            //    gsl_bspline_eval(rij, splineValuePair, bwPair);
            //    //size_t istart, iend;
            //    //gsl_bspline_eval_nonezero(rij, Bk, &istart, &iend, bw);
            //    //printf("rij = %lf, %d %d\n", rij, i, j);    
            //    //std::vector<double> dv = parallelVelocity(i, j);
            //    //Check if force matching works fine
	    //    //printf("rhoi = %lf rhoj = %lf\n", f_rhoi, f_rhoj);

            //    for (int m=0; m<info->numSplines; m++)
            //    {
            //         double phim = gsl_vector_get(splineValuePair, m);
            //         if (phim < 1e-20)
	    //             continue;
	    //    
	    //         // Evaluate the Local Density Modulating Terms

            //         // For all three dimensions
            //         frameMatrix[m][3*i]     += phim * dv[0];
            //         frameMatrix[m][3*i + 1] += phim * dv[1];
            //         frameMatrix[m][3*i + 2] += phim * dv[2];
            //         frameMatrix[m][3*j]     -= phim * dv[0];
            //         frameMatrix[m][3*j + 1] -= phim * dv[1];
            //         frameMatrix[m][3*j + 2] -= phim * dv[2];
            //    }
            //}  
        }
        //printf("done with %d\n", i+1);
    }
        
    // Constructing the normal Matrix and normal Vector
    for (int m=0; m<nSplines; m++)
    {
        for (int n=0; n<nSplines; n++)
        {
            double sum = 0.0;
            for (int k=0; k<3 * nall; k++)
                sum += frameMatrix[m][k] * frameMatrix[n][k];
            normalMatrix[m][n] += sum * normalFactor;
        }
   
        double sum_b = 0.0; 
        for (int k=0; k<3 * nall; k++)
            sum_b += frameMatrix[m][k] * trajFrame->residualForces[k/3][k%3];
        normalVector[m] += sum_b * normalFactor;
    }
}

// Accumulate the normal equation for pairwise interaction
void FitGLE::accumulateDensityNormalEquation()
{
    int nall = trajFrame->numParticles;
    int nSplines = info->numSplines;
    std::vector<std::vector<double> > frameMatrix(nSplines, std::vector<double>(3*nall));
    double normalFactor = 1.0 / info->steps;
   
    // Computing Matrix F_km 
    for (int i=0; i<nall; i++)
    {
        for (int j = i+1; j<nall; j++)
        {
            double rij = distance(trajFrame->positions[i], trajFrame->positions[j]);

	    if (rij < info->end)
            {
	        double rhoi = trajFrame->density[i];
	        double rhoj = trajFrame->density[j];
	            
	        // Calculate the f(rij) term
	        gsl_bspline_eval(rij, splineValuePair, bwPair);
	        double f_rij = 0;
                for (int m=0; m<info->numSplines; m++)
                {
                   f_rij += splineCoefficientsPair[m] * gsl_vector_get(splineValuePair, m);
                }
                //if (rij < 3.0) printf("rij = %lf fij = %lf\n", rij, f_rij);

                std::vector<double> dv = parallelUnitVector(i, j);

                gsl_bspline_eval(rhoi, splineValueDensity, bwDensity);
                for (int m=0; m<nSplines; m++)
                {
                     double phim = gsl_vector_get(splineValueDensity, m);
                     if (phim < 1e-20)
		         continue;
		
                     // For all three dimensions
                     frameMatrix[m][3*i]     += phim * dv[0] * f_rij * 0.5;
                     frameMatrix[m][3*i + 1] += phim * dv[1] * f_rij * 0.5;
                     frameMatrix[m][3*i + 2] += phim * dv[2] * f_rij * 0.5;
                     frameMatrix[m][3*j]     -= phim * dv[0] * f_rij * 0.5;
                     frameMatrix[m][3*j + 1] -= phim * dv[1] * f_rij * 0.5;
                     frameMatrix[m][3*j + 2] -= phim * dv[2] * f_rij * 0.5;
                }
                gsl_bspline_eval(rhoj, splineValueDensity, bwDensity);
                for (int m=0; m<nSplines; m++)
                {
                     double phim = gsl_vector_get(splineValueDensity, m);
                     if (phim < 1e-20)
		         continue;
		
                     // For all three dimensions
                     frameMatrix[m][3*i]     += phim * dv[0] * f_rij * 0.5;
                     frameMatrix[m][3*i + 1] += phim * dv[1] * f_rij * 0.5;
                     frameMatrix[m][3*i + 2] += phim * dv[2] * f_rij * 0.5;
                     frameMatrix[m][3*j]     -= phim * dv[0] * f_rij * 0.5;
                     frameMatrix[m][3*j + 1] -= phim * dv[1] * f_rij * 0.5;
                     frameMatrix[m][3*j + 2] -= phim * dv[2] * f_rij * 0.5;
                }
            }  
        }
    }
        
    // Constructing the normal Matrix and normal Vector
    for (int m=0; m<nSplines; m++)
    {
        for (int n=0; n<nSplines; n++)
        {
            double sum = 0.0;
            for (int k=0; k<3 * nall; k++)
                sum += frameMatrix[m][k] * frameMatrix[n][k];
            normalMatrix[m][n] += sum * normalFactor;
        }
   
        double sum_b = 0.0; 
        for (int k=0; k<3 * nall; k++)
            sum_b += frameMatrix[m][k] * trajFrame->residualForces[k/3][k%3];
        normalVector[m] += sum_b * normalFactor;
    }
}

void FitGLE::leastSquareSolver(bool isPair)
{
    // Solving the least square normal equation G*phi = b
    int nSplines = totalSplineCount;
    double* G = new double[nSplines * nSplines];
    double* b = new double[nSplines];

    
    // Preconditioning the Normal Matrix
    /*std::vector<double> h(info->numSplines, 0.0);
    for (int i = 0; i < info->numSplines; i++) {
        for (int j = 0; j < info->numSplines; j++) {
            h[j] = h[j] + normalMatrix[i][j] * normalMatrix[i][j];  //mat->dense_fm_matrix[j * mat->accumulation_matrix_rows + i] * mat->dense_fm_matrix[j * mat->accumulation_matrix_rows + i];
        }
    }
    for (int i = 0; i < info->numSplines; i++) {
        if (h[i] < 1.0E-20) h[i] = 1.0;
        else h[i] = 1.0 / sqrt(h[i]);
    }
    for (int i =0; i < info->numSplines; i++)
    {
        for (int j=0; j < info->numSplines; j++)
           normalMatrix[i][j] *= h[j];
    }*/


    // Store the normalMatrix in container 
    for (int m=0; m<nSplines; m++)
    {
        for (int n=0; n<nSplines; n++)
        {
            G[m*info->numSplines + n] = normalMatrix[m][n];
        }
        b[m] = normalVector[m];
        printf("m %d %lf\n", m, b[m]);
    }
    

    // Solving the linear system using SVD

    int m = nSplines;
    int n = nSplines;
    int nrhs = 1;
    int lda = nSplines;
    int ldb = 1;
    double rcond = -1.0;
    int irank;
    double* singularValue = new double[nSplines];
    int solverInfo = LAPACKE_dgelss(LAPACK_ROW_MAJOR, m, n, nrhs, G, lda, b, ldb, singularValue, rcond, &irank); 
   
    printf("LSQ Solver Info: %d\n", solverInfo);

    if (isPair)
    {
        for (int m=0; m<info->numSplines; m++)
        {
            splineCoefficientsPair[m] = b[m];
            splineCoefficientsPair[m] = b[m + info->numSplines];
        }
    }
    else
    {
        for (int m=0; m<info->numSplines; m++)
        {
            splineCoefficientsPair[m] = b[m];
            splineCoefficientsPair[m] = b[m + info->numSplines];
        }
    }

    delete[] G;
    delete[] b;
}

// Dump out the spline coefficients as well as fitted tables at each iteration
void FitGLE::output(int iter)
{
    printf("output ====== Iteration < %d >\n", iter);
    double start = info->start;
    double end = info->end;
    double precision = info->outputPrecision;

    // Note: Only Works for Unix Operating System! Change to Boost Filesystem might be able to correct this
    char cmd[50];
    sprintf(cmd, "mkdir iter_%d", iter);
    const int dir= system(cmd);
    if (dir< 0)
    {
        return;
    }

    char buffer_pairSpline[50];
    char buffer_pairTable[50];
    sprintf(buffer_pairSpline, "./iter_%d/pairsplinecoeff.dat", iter);
    sprintf(buffer_pairTable, "./iter_%d/pair_table.dat", iter);
    printf("%s\t%s\n", buffer_pairSpline, buffer_pairTable);

    FILE* fb = fopen(buffer_pairSpline, "w");
    for (int m=0; m<info->numSplines; m++)
    {
	printf("size = %lu, %lf\n", splineCoefficientsPair.size(), splineCoefficientsPair[m]);
        fprintf(fb, "%lf\n", splineCoefficientsPair[m]);
    }
    fclose(fb);

    FILE* fp = fopen(buffer_pairTable, "w");

    while (start < end)
    {
        double gamma_r = 0.0;
        gsl_bspline_eval(start, splineValuePair, bwPair);
        for (int m=0; m<info->numSplines; m++)
        {
           gamma_r += splineCoefficientsPair[m] * gsl_vector_get(splineValuePair, m);
        }
        fprintf(fp, "%lf\t%lf\n", start, gamma_r);
        start = start + precision;
    }
    
    fclose(fp);

    printf("finishing pair\n");

    // Deal with the density
    start = info->start_density;
    end = info->end_density;

    char buffer_densitySpline[50];
    char buffer_densityTable[50];
    sprintf(buffer_densitySpline, "./iter_%d/densitysplinecoeff.dat", iter);
    sprintf(buffer_densityTable, "./iter_%d/density_table.dat", iter);

    FILE* fdb = fopen(buffer_densitySpline, "w");
    for (int m=0; m<info->numSplinesDensity; m++)
    {
        fprintf(fdb, "%lf\n", splineCoefficientsDensity[m]);
    }
    fclose(fdb);

    FILE* fdp = fopen(buffer_densityTable, "w");

    while (start < end)
    {
        double gamma_r = 0.0;
        gsl_bspline_eval(start, splineValueDensity, bwDensity);
        for (int m=0; m<info->numSplinesDensity; m++)
        {
           gamma_r += splineCoefficientsDensity[m] * gsl_vector_get(splineValueDensity, m);
        }
        fprintf(fdp, "%lf\t%lf\n", start, gamma_r);
        start = start + precision;
    }
    
    fclose(fdp);
}

// Execution Process
void FitGLE::exec()
{
        for (int i=0; i<info->steps; i++)
        {
            trajFrame->readFrame();
            printf("finish reading frame\n");
            accumulatePairNormalEquation();
            printf("finishing step %d (total %d)\n", i+1, info->steps);
        }
        //leastSquareSolver(true);
	//output(0);
	printf("###########Answer#############\n");
        printf("%lf\n", normalVector[0]/normalMatrix[0][0]);
	printf("###########Byebye#############\n");
        clearNormalMatrix(); // Clear up the normal matrix for density fitting
}

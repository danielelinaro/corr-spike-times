// Generation of correlated spike trains
// R. Brette - Nov 2007

#ifndef CORRELATEDSPIKES_H_
#define CORRELATEDSPIKES_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct {
	int n; // neuron number
	double t; // spike timing
} event;

extern gsl_rng *rnd;

// -------------------------------------------------------
// Statistics of rectified Gaussian processes (appendix A)
// -------------------------------------------------------
void rectifiedGaussian(double mu,double sigma, double *mur,double *sigmar);
void invRectifiedGaussian(double mur,double sigmar, double *mu,double *sigma);

// -------------
// Cox Processes
// -------------
// Homogeneous pool
void homogeneousPoolCoxCD1(double r, double c, double tauc, double dt, int ndt, int n, int *spike);
void homogeneousPoolCoxCD2(double r, double c, double tauc, double dt, int ndt, int n, int *spike);
int homogeneousPoolCoxED(double r, double c, double tauc, double T, int n, event *spike, double *rate);

// General correlation structures
void decomposeCox(gsl_matrix *C,gsl_matrix *L);
void generalCoxCD(gsl_matrix *C, double tauc, double dt, int ndt, int n, int *spike, double *Xstore, double *Ystore);
int generalCoxED(gsl_matrix *C, double tauc, double T, int n, event *spike);

// -----------------
// Mixture Processes
// -----------------
// Calculation of mixture matrices
int globalSyncMixture(double *R,int n,double c,double *nu,gsl_matrix *P);
void topographicSyncMixture(double r,int n,double cmax, double alpha,double *nu,gsl_matrix *P);
int generalFixedColumnSumsMixture(double *r,int n,gsl_matrix *C,double *nu,gsl_matrix *P);
int generalOptimizationMixture(double *r,int n,gsl_matrix *C,double *nu,gsl_matrix *P);
// Random shifts
double randomLatency();
// Simulation of mixture processes
int offlineMixture(gsl_matrix *P,double *nu,int M, int N, double tauc, double T, event *spike);
int homogeneousPoolMixtureA(double r, double c, double tauc, double T, int n, event *spike);
int homogeneousPoolMixtureB(double r, double c, double tauc, double T, int n, event *spike);

#endif /*CORRELATEDSPIKES_H_*/

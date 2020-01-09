//============================================================================
// Name        : Correlations.cpp
// Author      : Romain Brette
// Version     :
// Copyright   : Your copyright notice
// Description : Draft of algorithms to generate correlated spike trains
//============================================================================

#include <cstdio>
#include <ctime>
#include <libgen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "erfinv.h"
#include "correlatedspikes.h"

double fWNIF(double x,void *params) {
	double y;
	
	y=exp(x*x)*(1.+gsl_sf_erf(x));
	return y;
}

// Firing rate of an integrate-and-fire model with white noise input current
double rateIFWhiteNoise(double refrac,double tau,double theta_eff,double H_eff) {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	
	double result, error;
	     
	gsl_function F;
	F.function = &fWNIF;
	F.params = NULL;
	     
	gsl_integration_qag(&F, H_eff, theta_eff, 0, 1e-7, 10000, GSL_INTEG_GAUSS15,w, &result, &error); 
	
	gsl_integration_workspace_free (w);
	
	return 1./(refrac+result*tau*sqrt(M_PI));
}

// Calculate time differences between the events of two spike trains (empirical)
int CCVF(event *spike,int k,double *timediff,double maxtime,int nr1,int nr2) {
	// spike: ordered list of events (with 2 spike trains nr1 & nr2)
	// k: number of events
	// timediff: array of time differences between the two spike trains
	// maxtime: maximum time difference stored
	// nr1: index of first neuron
	// nr2: index of second neuron
	// returns the number of time differences
	const int nstore=5; // number of previous events stored
	double previous0[nstore],previous1[nstore]; // circular array of stored events
	int index0,index1;
	int ndiff=0;
	
	index0=0;
	index1=0;
	for(int i=0;i<nstore;i++) {
		previous0[i]=previous1[i]=0.;
	}
	
	// First fill the circular arrays
	int n0=0; // number of total events in the first spike train
	int n1=0;
	int n; // index of processed event
	for(n=0;n<k;n++) {
		if (spike[n].n==nr1) {
			previous0[index0]=spike[n].t;
			index0++;
			if (index0==nstore)
				index0=0;
			n0++;
			// Check first if the circular array is filled
			if ((n0>=nstore) && (n1>=nstore)) {
				// Write time differences
				for(int j=0;j<nstore;j++) {
					timediff[ndiff]=spike[n].t-previous1[j];
					if (timediff[ndiff]<=maxtime)
						ndiff++;
				}
			}
		} else if (spike[n].n==nr2) {
			previous1[index1]=spike[n].t;
			index1++;
			if (index1==nstore)
				index1=0;
			n1++;
			// Check first if the circular array is filled
			if ((n0>=nstore) && (n1>=nstore)) {
				// Write time differences
				for(int j=0;j<nstore;j++) {
					timediff[ndiff]=-(spike[n].t-previous0[j]);
					if (-timediff[ndiff]<=maxtime)
						ndiff++;
				}
			}
		}
	}
	
	return ndiff;
}

void generate_spike_trains(int N, double tauc, double dur, double rate_range[2]) {
	gsl_matrix *C=gsl_matrix_calloc(N,N);
	double r[N];
	
	// Correlation and rate matrix
	for(int i=0;i<N;i++) {
		//r[i]=gsl_rng_uniform(rnd)*10.+5.; // 5..15 Hz
		r[i] = gsl_rng_uniform(rnd)*(rate_range[1]+rate_range[0])/2 + rate_range[0];
		for(int j=0;j<N;j++) {
			double a;
			if (j<i) {
				a=(1.+gsl_rng_uniform(rnd))*2.; // 10..20
				gsl_matrix_set(C,i,j,a);
				gsl_matrix_set(C,j,i,a);
			}
		}
	}
	
	// Decomposition
	double nu[2*N];
	gsl_matrix *P=gsl_matrix_calloc(N,2*N);
	if (generalOptimizationMixture(r,N,C,nu,P))
		fprintf(stderr, "oups\n");
	
	// Simulation
	event *spike=new event[1000000];
	int nspikes=offlineMixture(P,nu,2*N,N,tauc,dur,spike);
        fprintf(stderr, "%d spikes generated.\n", nspikes);
	
	// Spike trains
	for(int i=0;i<nspikes;i++)
		fprintf(stdout, "%.3f %d\n", spike[i].t, spike[i].n);
	
	delete [] spike;	
	
	gsl_matrix_free(C);
	gsl_matrix_free(P);
}

int main(int argc, char *argv[]) {
	// Initialization of random numbers
	const gsl_rng_type * T;
        int N;
        double tauc, dur, rate_range[2];
	
	gsl_rng_env_setup();
        gsl_rng_default_seed = time(NULL);
	T = gsl_rng_default;
	rnd = gsl_rng_alloc (T);

        if (argc != 6) {
                fprintf(stderr, "usage: %s N tauc duration rate_min rate_max\n", basename(argv[0]));
                fprintf(stderr, "\n");
                fprintf(stderr, "    N              number of presynaptic sources\n");
                fprintf(stderr, "    tauc           correlation time constant (ms)\n");
                fprintf(stderr, "    duration       interval of time over which spikes are generated (ms)\n");
                fprintf(stderr, "    rate_min       minimum presynaptic firing rate (Hz)\n");
                fprintf(stderr, "    rate_max       maximum presynaptic firing rate (Hz)\n");
                fprintf(stderr, "\n");
                exit(0);
        }
        
        N = atoi(argv[1]);
        tauc = atof(argv[2]);
        dur = atof(argv[3]);
        rate_range[0] = atof(argv[4]);
        rate_range[1] = atof(argv[5]);

        generate_spike_trains(N,tauc,dur,rate_range);

	gsl_rng_free (rnd);
	
	return 0;
}

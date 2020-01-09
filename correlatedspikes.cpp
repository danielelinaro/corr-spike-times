// Generation of correlated spike trains
// R. Brette - Nov 2007
#include <cstdio>
#include "correlatedspikes.h"

gsl_rng *rnd;

int compareEvent (const void * a, const void * b)
{
  if ( (*(event*)a).t > (*(event*)b).t )
	  return 1;
  else
	  return -1;
}

// -------------
// Cox Processes
// -------------

void rectifiedGaussian(double mu,double sigma, double *mur,double *sigmar) {
	// Calculate the mean and standard deviation for a rectified
	// Gaussian distribution
	// mu, sigma: parameters of the original distribution
	// *mur, *sigmar: parameters of the rectified distribution
	//
	// CHECKED	
	double a=1.+gsl_sf_erf(mu/(sqrt(2)*sigma));
	
	*mur=(sigma/sqrt(2.*M_PI))*exp(-0.5*gsl_pow_2(mu/sigma))+.5*mu*a;
	*sigmar=sqrt((mu-*mur)*(*mur)+.5*gsl_pow_2(sigma)*a);
}

void invRectifiedGaussian(double mur,double sigmar, double *mu,double *sigma) {
	// Inverse of function rectifiedGaussian
	//
	// CHECKED
	double m,s;
	double x0,x1,fx0,fx1,x2;
	int n;

	// Use the secant method to find the root
	// Initialization
	x0=mur/sigmar;
	x1=1.1*(mur/sigmar);
	rectifiedGaussian(x0,1.,&m,&s);
	fx1=m/s-mur/sigmar;
	// Iterations
	n=0;
	while ((fabs(fx1)>1e-4) && (n<1000)) {
		rectifiedGaussian(x1,1.,&m,&s);
		fx0=fx1;
		fx1=m/s-mur/sigmar;
		x2=x1;
		x1=x1-((x1-x0)/(fx1-fx0))*fx1;
		x0=x2;
		n++;
	}
	
	*sigma=mur/(exp(-0.5*x1*x1)/(sqrt(2.*M_PI))+.5*x1*(1.+gsl_sf_erf(x1/sqrt(2.))));
	*mu=x1*(*sigma);
}

// --------------------------------
// Homogeneous pool - Cox processes
// --------------------------------

void homogeneousPoolCoxCD1(double r, double c, double tauc, double dt, int ndt, int n, int *spike) {
	// Clock-driven generation of homogeneously correlated spike trains
	// (Cox processes)
	// r = rate (Hz)
	// c = total correlation strength (in [0,1])
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// ndt = number of timesteps
	// n = number of neurons
	// spike = array of 0/1 integers, each row is a timestep, each column is a neuron
	double sigmar,sigma,s,lambda,x,mu;
	int i,j;
	
	// Correction of mu and sigma
	sigmar=sqrt(c*r/(2.*tauc*0.001));
	invRectifiedGaussian(r,sigmar,&mu,&sigma);
	
	// Simulation
	x=0.;
	lambda=exp(-dt/tauc);
	s=sigma*sqrt(1-exp(-2.*dt/tauc));
	mu=mu*dt*0.001;
	s=s*dt*0.001;
	for(j=0;j<ndt;j++) {
		x=x*lambda+gsl_ran_gaussian(rnd,s);
		for(i=0;i<n;i++)
			if (gsl_rng_uniform(rnd)<x+mu) {
				// Spike
				spike[j*n+i]=1;
			}
	}
}

void homogeneousPoolCoxCD2(double r, double c, double tauc, double dt, int ndt, int n, int *spike) {
	// Clock-driven generation of homogeneously correlated spike trains
	// (Cox processes)
	// Second method (faster if n is large)
	//
	// r = rate (Hz)
	// c = total correlation strength (in [0,1])
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// ndt = number of timesteps
	// n = number of neurons
	// spike = array of 0/1 integers, each row is a timestep, each column is a neuron
	double sigmar,sigma,s,lambda,x,mu;
	int i,j,nsp;
	int *buffer;
	int *numbers;
	
	// Correction of mu and sigma
	sigmar=sqrt(c*r/(2.*tauc*0.001));
	invRectifiedGaussian(r,sigmar,&mu,&sigma);
	
	// Initialization
	buffer=new int[n];
	numbers=new int[n];
	for(i=0;i<n;i++) {
		buffer[i]=0;
		numbers[i]=i;
	}
	
	// Simulation
	x=0.;
	lambda=exp(-dt/tauc);
	s=sigma*sqrt(1-exp(-2.*dt/tauc));
	mu=mu*dt*0.001;
	s=s*dt*0.001;
	for(j=0;j<ndt;j++) {
		x=x*lambda+gsl_ran_gaussian(rnd,s);
		if (x+mu>0) {
			nsp=gsl_ran_binomial(rnd,x+mu,n); // number of spikes in the timestep
			gsl_ran_choose(rnd,buffer,nsp,numbers,n, sizeof(int));
			for(i=0;i<nsp;i++)
				spike[n*j+buffer[i]]=1;
		}
	}
	
	delete[] buffer;
	delete[] numbers;
}

int homogeneousPoolCoxED(double r, double c, double tauc, double T, int n, event *spike, double *rate) {
	// Event-driven generation of homogeneously correlated spike trains
	// (Cox processes)
	//
	// r = rate (Hz)
	// c = total correlation strength (in [0,1])
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// T = duration (ms)
	// n = number of neurons
	// spike = array of events (= couples (neuron,time))
	// rate (optional) = value of the rate at every event time (set NULL if unused) 
	//
	// result = number of events
	//
	// CHECKED
	double sigmar,sigma,x,mu,t,delta;
	int k;
	double M;
	
	// Correction of mu and sigma
	sigmar=sqrt(c*r/(2.*tauc*0.001));
	invRectifiedGaussian(r,sigmar,&mu,&sigma);

	mu=mu*0.001*n;
	sigma=sigma*0.001*n;
	M=mu+10.*sigma;
	
	x=0.;
	// Iterations
	k=0;
	for(t=(delta=gsl_ran_exponential(rnd,1./M));t<T;t=t+(delta=gsl_ran_exponential(rnd,1./M))) {
		x=x*exp(-delta/tauc)+gsl_ran_gaussian(rnd,sigma*sqrt(1.-exp(-2.*delta/tauc)));
		if (gsl_rng_uniform(rnd)*M<x+mu) {
			spike[k].n=gsl_rng_uniform_int(rnd,n);
			spike[k].t=t;
			if (rate)
				rate[k]=(x+mu)*1000./n;
			k++;
		}
	}
	//printf("k=%d\n",k);
	
	return(k);
}

// --------------------------------------
// Arbitrary correlations - Cox processes
// --------------------------------------
void decomposeCox(gsl_matrix *C,gsl_matrix *L) {
	// Complete the diagonal of C and
	// find L such that C=LL^T
	// C is matrix of correlation coefficients with unspecified diagonal
	// c_ii contains r_i (rate)
	// C must be symmetric
	//
	// CHECKED
	gsl_eigen_nonsymm_workspace *evwork;
	gsl_vector_complex *eval;
	double alpha;
	double m;
	
	// 0) Replace r_i by r_i^2 (var_i(x) is should have magnitude r_i^2)
	for(size_t i=0;i<C->size1;i++)
		gsl_matrix_set(C,i,i,gsl_matrix_get(C,i,i)*gsl_matrix_get(C,i,i));
	
	// Completion
	// 1) Calculate D^{-1}C
	for(size_t i=0;i<C->size1;i++) {
		for(size_t j=0;j<C->size2;j++)
			gsl_matrix_set(L,i,j,gsl_matrix_get(C,i,j)/gsl_matrix_get(C,i,i));
		gsl_matrix_set(L,i,i,0.);
	}
	// 2) Find the smallest eigenvalue
	evwork=gsl_eigen_nonsymm_alloc(C->size1);
	eval=gsl_vector_complex_alloc(C->size1);
	gsl_eigen_nonsymm(L,eval,evwork);
	m=1e6;
	for(int i=0;i<(int)eval->size;i++) {
		gsl_complex c=gsl_vector_complex_get(eval,i);
		if ((fabs(c.dat[1])<fabs(c.dat[0])*1e-4) && // real number
				(c.dat[0]<m))
			m=c.dat[0];
	}
	alpha=-m;
	gsl_eigen_nonsymm_free(evwork);
	gsl_vector_complex_free(eval);
	// 3) Complete the diagonal with alpha*ri^2
	gsl_matrix_memcpy(L,C);
	//alpha=alpha+.01; // avoids bad conditioning problems (uncomment if Cholesky fails)
	for(size_t i=0;i<C->size1;i++)
		gsl_matrix_set(L,i,i,alpha*gsl_matrix_get(C,i,i));
	
	if (gsl_linalg_cholesky_decomp(L)==GSL_EDOM) // failure
		printf("Oops - Cholesky failure!\n");

	// 4) Writes zeros above the diagonal
	for(size_t i=0;i<L->size1;i++)
		for(size_t j=i+1;j<L->size2;j++)
			gsl_matrix_set(L,i,j,0.);
}

void generalCoxCD(gsl_matrix *C, double tauc, double dt, int ndt, int n, int *spike, double *Xstore, double *Ystore) {
	// Clock-driven generation of correlated spike trains
	// (Cox processes)
	// C = correlation matrix (symetric); the diagonal entries are rates
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// ndt = number of timesteps
	// n = number of neurons (must equal the size of C)
	// spike = array of 0/1 integers, each row is a timestep, each column is a neuron
	// Xstore = values of X at each time step (set NULL if unused)
	// Ystore = values of Y (set NULL if unused)
	//
	// CHECKED
	double lambda,s;

	// Initialization
	gsl_vector *Y=gsl_vector_alloc(n);
	gsl_vector *X=gsl_vector_alloc(n);
	gsl_vector *R=gsl_vector_alloc(n); // Rates
	for(int i=0;i<n;i++) {
		gsl_vector_set(R,i,gsl_matrix_get(C,i,i));
		gsl_vector_set(Y,i,gsl_ran_gaussian(rnd,1.));
	}

	// Calculate L (Cholesky)
	gsl_matrix *L=gsl_matrix_alloc(n,n);
	decomposeCox(C,L);
	
	// Simulation
	lambda=exp(-dt/tauc);
	s=sqrt(1-exp(-2.*dt/tauc));
	for(int j=0;j<ndt;j++) {
		for(int i=0;i<n;i++)
			gsl_vector_set(Y,i,gsl_vector_get(Y,i)*lambda+gsl_ran_gaussian(rnd,s));
		
		// X=R+LY
		gsl_vector_memcpy(X,R);
		gsl_blas_dgemv(CblasNoTrans, 1., L, Y, 1., X);
		
		for(int i=0;i<n;i++) {
			if (Xstore)
				Xstore[j*n+i]=gsl_vector_get(X,i);
			if (Ystore)
				Ystore[j*n+i]=gsl_vector_get(Y,i);
			if (gsl_rng_uniform(rnd)<gsl_vector_get(X,i)*dt*0.001)
				// Spike
				spike[j*n+i]=1;
		}
	}
	
	// Free objects
	gsl_matrix_free(L);
	gsl_vector_free(X);
	gsl_vector_free(Y);
	gsl_vector_free(R);
}

int generalCoxED(gsl_matrix *C, double tauc, double T, int n, event *spike) {
	// Event-driven generation of correlated spike trains
	// (Cox processes)
	//
	// C = correlation matrix (symetric); the diagonal entries are rates
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// T = duration (ms)
	// n = number of neurons
	// spike = array of events (= couples (neuron,time))
	//
	// result = number of events
	int k;
	double M,sStar;

	// Rates
	double *R=new double[n];
	double *B=new double[n+1];
	M=0.;
	B[0]=0.;
	for(int i=0;i<n;i++) {
		R[i]=gsl_matrix_get(C,i,i);
		M+=R[i];
		B[i+1]=B[i]+R[i];
	}
	sStar=M;
	
	// Calculate L (Cholesky)
	gsl_matrix *L=gsl_matrix_alloc(n,n);
	double *A=new double[(n+1)*n];
	decomposeCox(C,L);
	for(int i=0;i<n;i++)
		A[i]=0.;
	for(int i=1;i<n+1;i++)
		for(int j=0;j<n;j++)
			A[i*n+j]=A[(i-1)*n+j]+gsl_matrix_get(L,i-1,j);
	
	// Calculate the upper bound
	double msum=0.;
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			for(int p=0;p<n;p++)
				msum+=gsl_matrix_get(L,i,p)*gsl_matrix_get(L,j,p);
	M+=10.*sqrt(msum);
	
	// Initialize variables
	double *Y=new double[n];
	double s;
	for(int i=0;i<n;i++)
		Y[i]=gsl_ran_gaussian(rnd,1.);
	
	// Iterations
	double delta;
	k=0;
	for(double t=(delta=1000.*gsl_ran_exponential(rnd,1./M));t<T;t=t+(delta=1000.*gsl_ran_exponential(rnd,1./M))) {
		s=sStar;
		for(int i=0;i<n;i++) {
			Y[i]=Y[i]*exp(-delta/tauc)+gsl_ran_gaussian(rnd,sqrt(1.-exp(-2.*delta/tauc)));
			s+=Y[i]*A[n*n+i];
		}
		
		if (gsl_rng_uniform(rnd)*M<s) {
			double v=gsl_rng_uniform(rnd)*s;
			// Binary search 
			int a=0;
			int b=n;
			int m;
			double xm;
			
			while (b-a>1) {
				m=(int)((a+b)/2);
				xm=B[m];
				for(int i=0;i<n;i++)
					xm+=Y[i]*A[m*n+i];
				if (v>xm)
					a=m;
				else
					b=m;
			}
			m=a;			
			spike[k].n=m;
			spike[k].t=t;
			k++;
		}
	}
	
	gsl_matrix_free(L);
	delete[] R;
	delete[] B;
	delete[] A;
	delete[] Y;
	
	return(k);
}

// -----------------
// Mixture processes
// -----------------
// Global synchronization
int globalSyncMixture(double *R,int n,double c,double *nu,gsl_matrix *P) {
	double sr=0.,sr2=0.;
	int err=0;
	
	for(int i=0;i<n;i++) {
		sr+=R[i];
		sr2+=R[i]*R[i];
	}
	*nu=sr2/(sr*c*n);
	
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			double x=c*R[i]*R[j]*n/sr2;
			gsl_matrix_set(P,i,j,x);
			if (x>1.) // not valid
				err=1;
		}
	}
	
	return err;
}

// Topographic synchronization
void topographicSyncMixture(double r,int n,double cmax, double alpha,double *nu,gsl_matrix *P) {
	// alpha = log a is the spatial constant
	double a=exp(-alpha);
	double lambda=2.*cmax*(1.+a*a)/(1.+a);
	
	*nu=r*(1.-a)/(2.*lambda);
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			gsl_matrix_set(P,i,j,lambda*exp(-fabs(i-j)*alpha));
		}
	}
}

// General matrix, fixed column sums
int generalFixedColumnSumsMixture(double *r,int n,gsl_matrix *C,double *nu,gsl_matrix *P) {
	// vector nu has 2*n entries
	// matrix P has n rows and 2*n columns
	//
	// TESTED
	
	// Complete C
	for(int i=0;i<n;i++) {
		double x=0;
		for(int j=0;j<n;j++)
			x+=gsl_matrix_get(C,i,j);
		gsl_matrix_set(C,i,i,-x+gsl_matrix_get(C,i,i));
	}

	// Find c0 (smallest eigenvalue)
	gsl_matrix *A=gsl_matrix_alloc(n,n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			gsl_matrix_set(A,i,j,gsl_matrix_get(C,i,j)/r[i]);
	
	gsl_eigen_nonsymm_workspace *evwork;
	gsl_vector_complex *eval;	
	evwork=gsl_eigen_nonsymm_alloc(n);
	eval=gsl_vector_complex_alloc(n);
	gsl_eigen_nonsymm(A,eval,evwork);
	double m=1e6;
	for(int i=0;i<(int)eval->size;i++) {
		gsl_complex c=gsl_vector_complex_get(eval,i);
		if ((fabs(c.dat[1])<fabs(c.dat[0])*1e-4) && // real number
				(c.dat[0]<m))
			m=c.dat[0];
	}
	gsl_eigen_nonsymm_free(evwork);
	gsl_vector_complex_free(eval);
	double c0=-m;
	
	// Completion
	gsl_matrix_memcpy(A,C);
	for(int i=0;i<n;i++)
		gsl_matrix_set(A,i,i,gsl_matrix_get(C,i,i)+c0*r[i]);
	
	// Decomposition
	if (gsl_linalg_cholesky_decomp(A)==GSL_EDOM) // failure
		printf("Oops - Cholesky failure!\n");

	// Write zeros above the diagonal
	for(size_t i=0;i<n;i++)
		for(size_t j=i+1;j<n;j++)
			gsl_matrix_set(A,i,j,0.);
	
	// Calculate source rates
	for(int i=0;i<n;i++) {
		m=0.;
		for(int j=0;j<n;j++) {
			double x=gsl_pow_2(gsl_matrix_get(A,j,i));
			if (x>m)
				m=x;
		}
		nu[i]=m;
	}
	
	// Calculate the mixture matrix
	for(int i=0;i<n;i++)
		for(int k=0;k<n;k++)
			gsl_matrix_set(P,i,k,gsl_matrix_get(A,i,k)/sqrt(nu[k]));
	
	// Completion of rates
	for(int i=0;i<n;i++) {
		double x=0.;
		for(int j=0;j<n;j++)
			x+=gsl_matrix_get(P,i,j)*nu[j];
		nu[i+n]=r[i]-x;
	}
	
	// Completion of the mixture matrix
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			if (i==j)
				gsl_matrix_set(P,i,j+n,1.);
			else
				gsl_matrix_set(P,i,j+n,0.);
	
	// Verification
	// 1) Positive rates
	int err=0;
	for(int i=0;i<2*n;i++)
		if (nu[i]<0.)
			err=1;
	// 2) Mixture matrix in [0,1]
	for(int i=0;i<n;i++)
		for(int j=0;j<2*n;j++)
			if ((gsl_matrix_get(P,i,j)<0.) || (gsl_matrix_get(P,i,j)>1.))
				err=1;
	
	gsl_matrix_free(A);
	
	return err;
}

// General matrix, optimization
// Uses a simple gradient descent (not the most efficient strategy)
int generalOptimizationMixture(double *r,int n,gsl_matrix *C,double *nu,gsl_matrix *P) {
	// vector nu has 2*n entries
	// matrix P has n rows and 2*n columns
	//
	// TESTED
	
	gsl_matrix *A=gsl_matrix_calloc(n,n);
	
	// Initialization
	for(int i=0;i<n;i++) {
		nu[i]=r[i];
		gsl_matrix_set(P,i,i,1.);
	}
	
	// Steps
	double b=0.01/n;
	double a=(1./n)*b;
	
	// Iterations
	double U[n];
	for(int ns=0;ns<20000;ns++) {
		// Calculate A
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
				if (i!=j) {
					double x=0.;
					for(int k=0;k<n;k++)
						x+=gsl_matrix_get(P,i,k)*gsl_matrix_get(P,j,k)*nu[k];
					gsl_matrix_set(A,i,j,x-gsl_matrix_get(C,i,j));
				} else
					gsl_matrix_set(A,i,j,0.);
		
		// H(P*nu-R)
		for(int i=0;i<n;i++) {
			double x=0.;
			for(int j=0;j<n;j++)
				x+=nu[j]*gsl_matrix_get(P,i,j);
			if (x>=r[i])
				U[i]=1.;
			else
				U[i]=0.;
		}
		
		// dPE, dPF
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++) {
				double x=0.;
				for(int k=0;k<n;k++)
					x+=gsl_matrix_get(A,i,k)*gsl_matrix_get(P,k,j)*nu[j];
				double y=gsl_matrix_get(P,i,j)-a*4*x-b*U[i]*nu[j];
				if (y<0.) // Clipping
					y=0.;
				if (y>1.)
					y=1.;
				gsl_matrix_set(P,i,j,y);
			}
		
		// dNuE, dNuF
		for(int i=0;i<n;i++) {
			double x=0.;
			for(int k=0;k<n;k++)
				for(int l=0;l<n;l++)
					x+=gsl_matrix_get(P,k,i)*gsl_matrix_get(P,l,i)*gsl_matrix_get(A,k,l);
			double y=0.;
			for(int j=0;j<n;j++)
				y+=U[j]*gsl_matrix_get(P,j,i);
			nu[i]-=a*x+b*y;
			if (nu[i]<0.) // Clipping
				nu[i]=0.;
		}		
	}
	
	// Completion of rates
	for(int i=0;i<n;i++) {
		double x=0.;
		for(int j=0;j<n;j++)
			x+=gsl_matrix_get(P,i,j)*nu[j];
		nu[i+n]=r[i]-x;
	}
	
	// Completion of the mixture matrix
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			if (i==j)
				gsl_matrix_set(P,i,j+n,1.);
			else
				gsl_matrix_set(P,i,j+n,0.);
	
	// Verification
	// 1) Positive rates
	int err=0;
	for(int i=0;i<2*n;i++) {
		if (nu[i]<0.) {
			err=1;
                        fprintf(stderr, "nu[%d] = %f < 0.\n", i, nu[i]);
                }
        }
	// 2) Mixture matrix in [0,1]
	for(int i=0;i<n;i++) {
		for(int j=0;j<2*n;j++) {
			if ((gsl_matrix_get(P,i,j)<0.) || (gsl_matrix_get(P,i,j)>1.)) {
				err=1;
                                fprintf(stderr, "P[%d,%d] = %f.\n", i, j, gsl_matrix_get(P,i,j));
                        }
                }
        }

	gsl_matrix_free(A);
	
	return err;
}

// Generate random latencies (spike time shifts), giving rise to
// exponential correlations with unit time constant
double randomLatency() {
	return (gsl_ran_exponential(rnd,1.));
}

// Offline simulation of a mixture process
int offlineMixture(gsl_matrix *P,double *nu,int M, int N, double tauc, double T, event *spike) {
	// Offline simulation of a general mixture process
	// Returns the number of events
	//
	// P = mixture matrix
	// nu = rate vector
	// M = number of source spike trains
	// n = number of target spike trains
	// tauc = correlation time constant
	// T = duration (ms)
	// spike = array of events (= couples (neuron,time))
	double rMean=0.;
	double wSize;
	int nspikes=0;

	// Calculate average target rate
	for(int i=0;i<N;i++)
		for(int k=0;k<M;k++)
			rMean+=gsl_matrix_get(P,i,k)*nu[k];
	rMean=rMean/N;
	
	// Optimal window size
	wSize=M*1./rMean;
	// Windowing is not implemented here
	wSize=T*0.001;
	
	// Generate M Poisson spike trains
	int nsource[M];
	double *sourceTrain[M];
	int ntarget;
	double *targetTrain = new double[1000000];
	for(int k=0;k<M;k++) {
		nsource[k]=gsl_ran_poisson (rnd, nu[k]*wSize); // number of spikes in train k
		sourceTrain[k]=new double[nsource[k]];
		for(int i=0;i<nsource[k];i++) {
			sourceTrain[k][i]=gsl_ran_flat(rnd,0.,wSize)*1000.; // in ms
		}
		// Select and shift spikes
		for(int i=0;i<N;i++) {
			ntarget=gsl_ran_binomial(rnd, gsl_matrix_get(P,i,k), nsource[k]);
			gsl_ran_choose(rnd,targetTrain, ntarget, sourceTrain[k], nsource[k], sizeof(double));
			for(int j=0;j<ntarget;j++) {
				spike[nspikes].n=i;
				spike[nspikes].t=targetTrain[j]+randomLatency()*tauc;
				nspikes++;
			}
		}
	}
	
	// Sort the events
	qsort(spike,nspikes,sizeof(event),compareEvent);
	
        delete [] targetTrain;

	for(int k=0;k<M;k++)
		delete [] sourceTrain[k];
		
	return nspikes;
}

int homogeneousPoolMixtureA(double r, double c, double tauc, double T, int n, event *spike) {
	// Event-driven generation of homogeneously correlated spike trains
	// (Mixture processes - A)
	//
	// r = rate (Hz)
	// c = total correlation strength (in [0,1])
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// T = duration (ms)
	// n = number of neurons
	// spike = array of events (= couples (neuron,time))
	// rate (optional) = value of the rate at every event time (set NULL if unused) 
	//
	// result = number of events
	//
	gsl_matrix *P=gsl_matrix_alloc(n,1);
	double nu=r/c;
	for(int i=0;i<n;i++)
		gsl_matrix_set(P,i,0,c);
	
	int nevents=offlineMixture(P,&nu,1,n,tauc,T,spike);
	
	gsl_matrix_free(P);
	
	return nevents;
}

int homogeneousPoolMixtureB(double r, double c, double tauc, double T, int n, event *spike) {
	// Event-driven generation of homogeneously correlated spike trains
	// (Mixture processes - B)
	//
	// r = rate (Hz)
	// c = total correlation strength (in [0,1])
	// tauc = correlation time constant (ms)
	// dt = timestep (ms)
	// T = duration (ms)
	// n = number of neurons
	// spike = array of events (= couples (neuron,time))
	// rate (optional) = value of the rate at every event time (set NULL if unused) 
	//
	// result = number of events
	//
	gsl_matrix *P=gsl_matrix_calloc(n,n+1);
	double nu[n+1];
	for(int i=0;i<n;i++) {
		nu[i]=(1-c)*r;
		gsl_matrix_set(P,i,i,1.);
		gsl_matrix_set(P,i,n,1.);
	}
	nu[n]=c*r;
	
	int nevents=offlineMixture(P,nu,n+1,n,tauc,T,spike);
	
	gsl_matrix_free(P);
	
	return nevents;
}

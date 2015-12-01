#include <math.h>
#include "matrix.h"
#include "mex.h"   //--This one is required
#include <float.h>
#include <stdio.h>
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))


double termCD(double**** u, const int* N, int i0, int i1, int i2, int i3, double h, double va, double vd) {
// Computes upwind term
double tCD, abs_Dux, abs_Duy;
double Dux0, Dux1, Duy0, Duy1;
int N0, N1, N2, N3;

N0 = N[0];
N1 = N[1];
N2 = N[2];
N3 = N[3];

if (i0>0 && i0<N0-1) { 		Dux0 = (u[i0+1][i1][i2][i3] - u[i0-1][i1][i2][i3])/2.0/h;
} else if (i0 == N0-1) { 	Dux0 = (u[i0][i1][i2][i3] - u[i0-1][i1][i2][i3])/h; 
} else if (i0 == 0) { 		Dux0 = (u[i0+1][i1][i2][i3] - u[i0][i1][i2][i3])/h; }

if (i1>0 && i1<N1-1) { 		Dux1 = (u[i0][i1+1][i2][i3] - u[i0][i1-1][i2][i3])/2.0/h;
} else if (i1 == N1-1) { 	Dux1 = (u[i0][i1][i2][i3] - u[i0][i1-1][i2][i3])/h; 
} else if (i1 == 0) {		Dux1 = (u[i0][i1+1][i2][i3] - u[i0][i1][i2][i3])/h; }

if (i2>0 && i2<N2-1) { 		Duy0 = (u[i0][i1][i2+1][i3] - u[i0][i1][i2-1][i3])/2.0/h;
} else if (i2 == N2-1) { 	Duy0 = (u[i0][i1][i2][i3] - u[i0][i1][i2-1][i3])/h;
} else if (i2 == 0) { 		Duy0 = (u[i0][i1][i2+1][i3] - u[i0][i1][i2][i3])/h; }

if (i3>0 && i3<N3-1) { 		Duy1 = (u[i0][i1][i2][i3+1] - u[i0][i1][i2][i3-1])/2.0/h;
} else if (i3 == N3-1) { 	Duy1 = (u[i0][i1][i2][i3] - u[i0][i1][i2][i3-1])/h;
} else if (i3 == 0) { 		Duy1 = (u[i0][i1][i2][i3+1] - u[i0][i1][i2][i3])/h; }

abs_Dux = sqrt( Dux0*Dux0 + Dux1*Dux1 );
abs_Duy = sqrt( Duy0*Duy0 + Duy1*Duy1 );

tCD = -(-va*abs_Dux + vd*abs_Duy);

return tCD;
}



double termDiss(double**** u, const int* N, int i0, int i1, int i2, int i3, double h, double va, double vd) {
// 1-sided derivatives
double Dux0p, Dux0m, Dux1p, Dux1m;
double Duy0p, Duy0m, Duy1p, Duy1m;

// Dissipation scaling
double a0, a1, a2, a3;

// Central and absolute derivatives
double Dux0, Dux1, Duy0, Duy1;
double abs_Dux, abs_Duy;

// Dissipation term
double tDiss;

int N0, N1, N2, N3;

N0 = N[0];
N1 = N[1];
N2 = N[2];
N3 = N[3];

if (i0 < N0-1) {	Dux0p = (u[i0+1][i1][i2][i3] - u[i0][i1][i2][i3])/h;
} else { 			Dux0p = (u[i0][i1][i2][i3] - u[i0-1][i1][i2][i3])/h; }

if (i0 > 0) {		Dux0m = (u[i0][i1][i2][i3] - u[i0-1][i1][i2][i3])/h;
} else { 			Dux0m = (u[i0+1][i1][i2][i3] - u[i0][i1][i2][i3])/h; }

if (i1 < N1-1) {	Dux1p = (u[i0][i1+1][i2][i3] - u[i0][i1][i2][i3])/h;
} else { 			Dux1p = (u[i0][i1][i2][i3] - u[i0][i1-1][i2][i3])/h; }

if (i1 > 0) {		Dux1m = (u[i0][i1][i2][i3] - u[i0][i1-1][i2][i3])/h;
} else {			Dux1m = (u[i0][i1+1][i2][i3] - u[i0][i1][i2][i3])/h; }

if (i2 < N2-1) {	Duy0p = (u[i0][i1][i2+1][i3] - u[i0][i1][i2][i3])/h;
} else { 			Duy0p = (u[i0][i1][i2][i3] - u[i0][i1][i2-1][i3])/h; } 

if (i2 > 0) {		Duy0m = (u[i0][i1][i2][i3] - u[i0][i1][i2-1][i3])/h;
} else { 			Duy0m = (u[i0][i1][i2+1][i3] - u[i0][i1][i2][i3])/h; }

if (i3 < N3-1) {	Duy1p = (u[i0][i1][i2][i3+1] - u[i0][i1][i2][i3])/h;
} else {			Duy1p = (u[i0][i1][i2][i3] - u[i0][i1][i2][i3-1])/h; }

if (i3 > 0) {		Duy1m = (u[i0][i1][i2][i3] - u[i0][i1][i2][i3-1])/h;
} else { 			Duy1m = (u[i0][i1][i2][i3+1] - u[i0][i1][i2][i3])/h; }

// Dissipation scaling
Dux0 = 0.5*(Dux0p + Dux0m);
Dux1 = 0.5*(Dux1p + Dux1m);
Duy0 = 0.5*(Duy0p + Duy0m);
Duy1 = 0.5*(Duy1p + Duy1m);
abs_Dux = sqrt( Dux0*Dux0 + Dux1*Dux1 );
abs_Duy = sqrt( Duy0*Duy0 + Duy1*Duy1 );

if (abs_Dux>0.000001) { 	
	a0 = va*abs(Dux0)/abs_Dux; 
	a1 = va*abs(Dux1)/abs_Dux;
} else { 			
	a0 = 0.0; 
	a1 = 0.0;
}

if (abs_Duy>0.000001) { 	
	a2 = vd*abs(Duy0)/abs_Duy;
	a3 = vd*abs(Duy1)/abs_Duy;
} else {
	a2 = 0.0;
	a3 = 0.0;
}

// Dissipation term
tDiss = a0*(Dux0p - Dux0m)/2.0 + a1*(Dux1p - Dux1m)/2.0 + a2*(Duy0p - Duy0m)/2.0 + a3*(Duy1p - Duy1m)/2.0;

return tDiss;
}

void levelSet(double**** u, double**** u0, double**** u1, double**** a, const int* N, double h, double k, int T, double va, double vd) {
int i0, i1, i2, i3, N0, N1, N2, N3;
double tDiss, tCD;
int Tstep, Tsubstep;

N0 = N[0];
N1 = N[1];
N2 = N[2];
N3 = N[3];

for (Tstep=1; Tstep<=T; Tstep++) {
	printf("Step %d of %d\n", Tstep, T);
	
	// ----- FIRST SUBSTEP -----
	for (i0=0; i0<N0; i0++) {
		for (i1=0; i1<N1; i1++) {
			for (i2=0; i2<N2; i2++) {
				for (i3=0; i3<N3; i3++) {
					
					// Upwind and diffusion terms
					tCD = termCD(u, N, i0, i1, i2, i3, h, va, vd);
					tDiss = termDiss(u, N, i0, i1, i2, i3, h, va, vd);

					// Update u0 from u
					u0[i0][i1][i2][i3] = u[i0][i1][i2][i3] + min(0.0, - k*(tCD - tDiss));
					
					// Obstacles
					u0[i0][i1][i2][i3] = max(u0[i0][i1][i2][i3],-a[i0][i1][i2][i3]);
				}
			}
		}
	}
	
	// ----- SECOND SUBSTEP -----
	for (i0=0; i0<N0; i0++) {
		for (i1=0; i1<N1; i1++) {
			for (i2=0; i2<N2; i2++) {
				for (i3=0; i3<N3; i3++) {
					
					// Upwind and diffusion terms
					tCD = termCD(u0, N, i0, i1, i2, i3, h, va, vd);
					tDiss = termDiss(u0, N, i0, i1, i2, i3, h, va, vd);

					// Update u1 from u0
					u1[i0][i1][i2][i3] = u0[i0][i1][i2][i3] + min(0.0, - k*(tCD - tDiss));
					
					// Obstacles
					u1[i0][i1][i2][i3] = max(u1[i0][i1][i2][i3],-a[i0][i1][i2][i3]);
				}
			}
		}
	}
	
	// Update u by averaging u and u1
	for (i0=0; i0 < N0; i0++) {
		for (i1=0; i1 < N1; i1++) {
			for (i2=0; i2 < N2; i2++) {
				for (i3=0; i3<N3; i3++) {
					u[i0][i1][i2][i3] = 0.5*(u[i0][i1][i2][i3]+u1[i0][i1][i2][i3]);
				}
			}
		}
	}	
	

}

} // End of function


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //---Inside mexFunction---
    //Declarations
    double *uValues, *aValues;
    double ****u, ****u0, ****u1, ****a;
    const int *N; 
	int i0, i1, i2, i3, N0, N1, N2, N3; // 
    double h,k; // grid spacing in x and t
	double va, vd; // speeds
	int T;
	
	
    //Get the input
    uValues = mxGetPr(prhs[0]);
	aValues = mxGetPr(prhs[1]);
	h       = mxGetScalar(prhs[2]);
	k 		= mxGetScalar(prhs[3]);
	T 		= mxGetScalar(prhs[4]);
	va		= mxGetScalar(prhs[5]);
	vd		= mxGetScalar(prhs[6]);

	N    = mxGetDimensions(prhs[0]);
	N0 = N[0];
	N1 = N[1];
	N2 = N[2];
	N3 = N[3];

    // memory allocation for u
    u    = (double ****) malloc ( N0 * sizeof(double***));
	u0    = (double ****) malloc ( N0 * sizeof(double***));
	u1    = (double ****) malloc ( N0 * sizeof(double***));
	a    = (double ****) malloc ( N0 * sizeof(double***));
	for (i0=0; i0<N0; i0++) {
		u[i0]    = (double ***) malloc (N1 * sizeof(double**));
		u0[i0]    = (double ***) malloc (N1 * sizeof(double**));
		u1[i0]    = (double ***) malloc (N1 * sizeof(double**));
		a[i0]    = (double ***) malloc (N1 * sizeof(double**));
		for (i1=0; i1<N1; i1++) {
			u[i0][i1]    = (double **) malloc (N2 * sizeof(double*));
			u0[i0][i1]    = (double **) malloc (N2 * sizeof(double*));
			u1[i0][i1]    = (double **) malloc (N2 * sizeof(double*));
			a[i0][i1]    = (double **) malloc (N2 * sizeof(double*));
			for (i2=0; i2<N2; i2++) {
				u[i0][i1][i2]    = (double *) malloc (N3 * sizeof(double));
				u0[i0][i1][i2]    = (double *) malloc (N3 * sizeof(double));
				u1[i0][i1][i2]    = (double *) malloc (N3 * sizeof(double));
				a[i0][i1][i2]    = (double *) malloc (N3 * sizeof(double));			
			}
        }
	}
	
	// printf("%d %d %d\n", N0, N1, T);
    // assignment u
    for (i0=0; i0 < N0; i0++) {
		for (i1=0; i1 < N1; i1++) {
			for (i2=0; i2 < N2; i2++) {
				for (i3=0; i3<N3; i3++) {
						u[i0][i1][i2][i3]    =    uValues[((i3*N2+i2)*N1 + i1)*N0 + i0];
						u0[i0][i1][i2][i3]    =    uValues[((i3*N2+i2)*N1 + i1)*N0 + i0];
						u1[i0][i1][i2][i3]    =    uValues[((i3*N2+i2)*N1 + i1)*N0 + i0];
						a[i0][i1][i2][i3]    =    aValues[((i3*N2+i2)*N1 + i1)*N0 + i0];
				}
            }
		}
	}
    
	// Compute u
	levelSet(u, u0, u1, a, N, h, k, T, va, vd);
	
    for (i0=0; i0 < N0; i0++) {
		for (i1=0; i1 < N1; i1++) {
			for (i2=0; i2 < N2; i2++) {
				for (i3=0; i3<N3; i3++) {
					uValues[((i3*N2+i2)*N1 + i1)*N0 + i0]  = u[i0][i1][i2][i3];
					// u0Values[((i3*N2+i2)*N1 + i1)*N0 + i0]  = u0[i0][i1][i2][i3];
				}
            }
        }
	}
	
    // delete u
    for (i0=0; i0 < N0; i0++) {
		for (i1=0; i1 < N1; i1++) {
			for (i2=0; i2 < N2; i2++) {
				free(u[i0][i1][i2]);
				free(u0[i0][i1][i2]);
				free(u1[i0][i1][i2]);
				free(a[i0][i1][i2]);			
			}
			free(u[i0][i1]);
			free(u0[i0][i1]);
			free(u1[i0][i1]);
			free(a[i0][i1]);
        }
		free(u[i0]);
		free(u0[i0]);
		free(u1[i0]);
		free(a[i0]);
	}
	free(u);
	free(u0);
	free(u1);
	free(a);
}

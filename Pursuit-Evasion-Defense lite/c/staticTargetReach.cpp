#include <math.h>
#include "matrix.h"
#include "mex.h"   //--This one is required
#include <float.h>
#include <stdio.h>
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

double termEO(double uxm, double uxy, double uxp, double uym, double uyp, double h) {
// Computes upwind term
double tEO;
double Duxm, Duxp, Duym, Duyp;
double Dux, abs_Dux, abs_Duy, abs_Du;

Duxm = (uxy-uxm)/h;
Duxp = (uxp-uxy)/h;
Duym = (uxy-uym)/h;
Duyp = (uyp-uxy)/h;

Dux = max(Duxm,0.0)*max(Duxm,0.0) + min(Duxp,0.0)*min(Duxp,0.0);
abs_Dux = sqrt( max(Duxm,0.0)*max(Duxm,0.0) + min(Duxp,0.0)*min(Duxp,0.0) );
abs_Duy = sqrt( max(Duym,0.0)*max(Duym,0.0) + min(Duyp,0.0)*min(Duyp,0.0) );

abs_Du = sqrt(max(Duxm,0.0)*max(Duxm,0.0) + min(Duxp,0.0)*min(Duxp,0.0) + max(Duym,0.0)*max(Duym,0.0) + min(Duyp,0.0)*min(Duyp,0.0));

//if (1) {
if ( (Dux<=0) && (abs_Duy<=abs_Dux) ) {
	tEO = abs_Du;
} else {
	tEO = -1.0/sqrt(2.0) * Dux + 1.0/sqrt(2.0) * abs_Duy;
}

return tEO;
}

double termCD(double uxm, double uxy, double uxp, double uym, double uyp, double h) {
// Computes upwind term
double tCD;
double Dux, Duy, abs_Dux, abs_Duy, abs_Du;

Dux = (uxp-uxm)/2.0/h;
Duy = (uyp-uym)/2.0/h;

abs_Dux = sqrt( Dux*Dux );
abs_Duy = sqrt( Duy*Duy );

abs_Du = sqrt( Dux*Dux + Duy*Duy );

if (1) {
// if ( (Dux<=0) && (abs_Duy<=abs_Dux) ) {
	tCD = abs_Du;
} else {
	tCD = -1.0/sqrt(2.0) * Dux + 1.0/sqrt(2.0) * abs_Duy;
}

return tCD;
}


double termDiss(double uxm, double uxy, double uxp, double uym, double uyp, double upp, double umm, double upm, double ump, double h) {
double Duxp, Duxm, Duyp, Duym;
double a1, a2, Dux, Duy, abs_Du;
double tDiss;

/* STENCIL
ump uyp upp
uxm uxy uxp
umm uym upm
*/

Duxp = (uxp-uxy)/h;
Duxm = (uxy-uxm)/h;
Duyp = (uyp-uxy)/h;
Duym = (uxy-uym)/h;

// Dissipation scaling
Dux = 0.5*(Duxp + Duxm);
Duy = 0.5*(Duyp + Duym);
abs_Du = sqrt( Dux*Dux + Duy*Duy );

a1 = abs(Dux/abs_Du);
a2 = abs(Duy/abs_Du);

Duxp = (uxp-uxy)/h;
Duxm = (uxy-uxm)/h;
Duyp = (uyp-uxy)/h;
Duym = (uxy-uym)/h;

tDiss = (Duxp - Duxm)/2.0 + (Duyp - Duym)/2.0;

return tDiss;
}

void levelSet(double*** u, const int* N, double h, double k, double epsilon, double A) {
int i, j, n, Ni, Nj, T;
double uxm, uxy, uxp, uym, uyp, umm, ump, upm, upp;
double tHeat, tEO, tDiss, tCD;
int im, ip, jm, jp;

Ni = N[0];
Nj = N[1];
T = N[2];

for (n=0; n<T-1; n++){
	for (i=0; i<Ni; i++) {
		for (j=0; j<Nj; j++) {
			/* STENCIL
			ump uyp upp
			uxm uxy uxp
			umm uym upm
			*/
			
			if (i==0) { im = 1;
			} else { 	im = i-1; }
			
			if (i==Ni-1) { 	ip = Ni-2;
			} else { 		ip = i+1; }	
			
			if (j==0) { jm = 1;
			} else { 	jm = j-1; }
			
			if (j==Nj-1) { 	jp = Nj-2;
			} else { 	   	jp = j+1; }	
			
			
			uxm = u[im][j][n];
			uxy = u[i][j][n];
			uxp = u[ip][j][n];
			
			uym = u[i][jm][n];
			uyp = u[i][jp][n];
			
			umm = u[im][jm][n];
			upm = u[ip][jm][n];
			ump = u[im][jp][n];
			upp = u[ip][jp][n];
			
			// Upwind and diffusion terms
			// tEO = termEO(uxm, uxy, uxp, uym, uyp, h);
			// tHeat = termHeat(uxm, uxy, uxp, uym, uyp, upp, umm, upm, ump, h);
			tCD = termCD(uxm, uxy, uxp, uym, uyp, h);
			tDiss = termDiss(uxm, uxy, uxp, uym, uyp, upp, umm, upm, ump, h);

			// Update
			//u[i][j][n+1] = u[i][j][n] - k*A*tEO + k*epsilon*tHeat;
			// u[i][j][n+1] = u[i][j][n] - k*A*tEO;
			u[i][j][n+1] = u[i][j][n] - k*A*(tCD - tDiss);
			// u[i][j][n+1] = min(u[i][j][n+1], u[i][j][n]);

		}
	}
}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //---Inside mexFunction---
    //Declarations
    double *uValues, *lValues;
    double ***u;
    const int *N; 
	int i, j, n, Ni, Nj, T; // jth grid in x, nth grid in t
    double h,k; // grid spacing in x and t
	double A, epsilon;
	
    //Get the input
    uValues = mxGetPr(prhs[0]);
	h       = mxGetScalar(prhs[1]);
	k 		= mxGetScalar(prhs[2]);
	A		= mxGetScalar(prhs[3]);
	epsilon = mxGetScalar(prhs[4]);
	
	N    = mxGetDimensions(prhs[0]);
	Ni = N[0];
	Nj = N[1];
	T = N[2];

    // memory allocation for u
    u    = (double ***) malloc ( Ni * sizeof(double**));
	for (i=0;i<Ni;i++) {
		u[i]    = (double **) malloc (Nj * sizeof(double*));
		for (j=0;j<Nj;j++) {
			u[i][j]    = (double *) malloc (T * sizeof(double));
        }
	}
	
	// printf("%d %d %d\n", Ni, Nj, T);
    // assignment u
    for (i=0; i < Ni; i++) {
		for (j=0; j < Nj; j++) {
			for (n=0; n < T; n++) {
                u[i][j][n]    =    uValues[(n*Nj + j)*Ni + i];
            }
		}
	}
    printf("2\n");
	// Compute u
	levelSet(u, N, h, k, epsilon, A);
	
    for (i=0; i < Ni; i++) {
        for (j=0; j < Nj; j++) {
            for (n=0; n < T; n++) {
                uValues[(n*Nj+j)*Ni+i]  = u[i][j][n];
            }
        }
	}
	
    // delete u
	for(i=0; i< Ni; i++) {
		for(j=0; j< Nj; j++) {
			free(u[i][j]);
		}
		free(u[i]);
	}
	free(u);
}

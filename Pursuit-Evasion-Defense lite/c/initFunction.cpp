#include <math.h>
#include "matrix.h"
#include "mex.h"   //--This one is required
#include <float.h>
#include <stdio.h>
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

double inCurve(double *x, double *y, int length, double x0, double y0) {
// Computes whether point is inside or outside of a closed curve
// Indicator variable; -1 means inside, 1 means outside

double integral, dx, dy;
int i;

integral = 0.0;

for (i=0; i<length; i++) {
	// Point (x0, y0) is exactly on the curve
	if ((x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0) == 0.0) {
		return -1.0;
	}
	
	// Point (x0, y0) is inside or outside the curve
	if (i>0) {
		dx = x[i]-x[i-1];
		dy = y[i]-y[i-1];
	} else if (i==0) {
		dx = x[0]-x[length-1];
		dy = y[0]-y[length-1];
	}
	
	integral += (-(y[i]-y0)*dx + (x[i]-x0)*dy) / ((x[i]-x0)*(x[i]-x0) + (y[i]-y0)*(y[i]-y0));
	
}

if (abs(integral) > 5.0) {
	return -1.0;
} else {
	return 1.0;
}

}

double dist2Segment(double x0, double x1, double y0, double y1, double x, double y) {
// Computes the distance from a point (x,y) to the line segment (x0, y0)-(x1, x1)
double dist;

// P0 = (x0, y0), P1 = (x1, y1), P = (x, y)

if ((x1-x0)*(x-x0) + (y1-y0)*(y-y0) < 0.0) { 			// Angle between P0->P1 and P0->P > 90 degrees
	// Distance = P0->P
	dist = sqrt((x0-x)*(x0-x) + (y0-y)*(y0-y));
} else if ((x0-x1)*(x-x1) + (y0-y1)*(y-y1) < 0.0) {   // Angle between P1->P0 and P1->P > 90 degrees
	// Distance = P1->P
	dist = sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));
} else {
	// Distance = P->(P0-P1)
	dist = abs( (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0) ) / sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
}

return dist;
}

double minDist(double *x, double *y, int length, double x0, double y0) {
// Computes minimum distance from (x0, y0) to closed curve
int i;
double mDist, disti;
mDist = 1000000.0; // Initial very large value

for (i=0; i<length; i++) {
	if (i>0) {
		disti = dist2Segment(x[i], x[i-1], y[i], y[i-1], x0, y0);
	} else {
		disti = dist2Segment(x[0], x[length-1], y[0], y[length-1], x0, y0);
	}
	mDist = min(disti, mDist);
}

return mDist;
}

void createIC(double *x, double *y, int length, const int *N, double** xs, double** ys, double** u) {
int i, j, Ni, Nj;
double x0, y0;

Ni = N[0];
Nj = N[1];

for (i=0; i<Ni; i++){
	for (j=0; j<Nj; j++) {
		x0 = xs[i][j];
		y0 = ys[i][j];
		u[i][j] = inCurve(x, y, length, x0, y0)*minDist(x, y, length, x0, y0);
	}
}
	
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //---Inside mexFunction---
    //Declarations
    double *xValues, *yValues, *xsValues, *ysValues, *uValues;
	double *x, *y;
    double **xs, **ys, **u;
	double *s;
    const int *dims, *N; 
    int i, j, Ni, Nj, length; // jth grid in x, nth grid in t
	
    //Get the input
    xValues 	= mxGetPr(prhs[0]);
	yValues 	= mxGetPr(prhs[1]);
	xsValues 	= mxGetPr(prhs[2]);
	ysValues 	= mxGetPr(prhs[3]);
	uValues  	= mxGetPr(prhs[4]);
	
    dims    = mxGetDimensions(prhs[0]);
	length 	= dims[0];
	N       = mxGetDimensions(prhs[4]);
	Ni 		= N[0];
	Nj 		= N[1];
	
    // memory allocation for u
	xs = (double **) malloc ( Ni * sizeof(double*));
	ys = (double **) malloc ( Ni * sizeof(double*));
	u = (double **) malloc ( Ni * sizeof(double*));
	
	for (i=0; i<Ni; i++){
		xs[i] = (double *) malloc (Nj * sizeof(double));
		ys[i] = (double *) malloc (Nj * sizeof(double));
		u[i] = (double *) malloc (Nj * sizeof(double));
    }

    // Initialize u
    for (i=0; i < Ni; i++) {
		for (j=0; j < Nj; j++) {
			xs[i][j] = xsValues[j*Ni+i];
			ys[i][j] = ysValues[j*Ni+i];
			u[i][j] =  uValues[j*Ni+i];
        }
	}
    
	// memory allocation for x and y list
	x = (double *) malloc ( length * sizeof(double));
	y = (double *) malloc ( length * sizeof(double));

    // Initialize x and y list
    for (i=0; i < length; i++) {
		x[i] = xValues[i];
		y[i] = yValues[i];
	}	
	
	// Compute u
	createIC(x, y, length, N, xs, ys, u);

    // send the processed u to the output  
    for (i=0; i<Ni; i++) {
		for (j=0; j<Nj; j++) {
            uValues[(j*Ni)+i] = u[i][j];
		}
	}
    
    // delete u;
	for(i=0; i< Ni; i++){
		free(xs[i]);
		free(ys[i]);
		free(u[i]);
  	}
	free(xs);
	free(ys);
	free(u);
	
	free(x);
	free(y);
   
}

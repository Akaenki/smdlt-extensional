//
//  Header.h
//  KRBC
//
//  Created by Linling Miao on 5/1/16.
//  Copyright © 2017 Linling Miao. All rights reserved.
//

#ifndef Header_h
#define Header_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cblas.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define EP 0.9624
#define THETA0 31.7

#ifndef M_PI
#define M_PI 3.1415926535
#endif

/////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    double x,y,z;
} Vector3D_t;

typedef struct {
    double rx,ry,rz;
    double fx,fy,fz;
} Chain_t;

/* Mobility Matrix */
typedef struct {
    float *DiffMatrix; // binx*biny*binz*i*j*9
    float *selfMatrix; // i*j*9
    uint32_t *gridCounting; // binx*biny*binz
    float *DiffMatrixAvg; // preaveraged mobility matrix
    float *selfMatrixAvg; //preaveraged self matrix
    uint32_t *counter; //counter for preaveraged matrix
} Matrix_t;

/* Geyer-Winter Parameters */
typedef struct {
    double beta_ij;
    float *C_i;
} gwParm_t;

/* bin tag */
typedef struct {
    uint64_t x,y,z;
} bin_t;

/////////////////////////////////////////////////////////////////////////////
// Gloable Variables/Arrays
/////////////////////////////////////////////////////////////////////////////
long *idum; //SEED

int HImode;
int numThreads;

double c_normal, flowRate;
char *c_normalc, *flowRatec, *directory;

uint16_t iteration;
uint32_t N,NP;
uint64_t binxmax,binymax,binzmax; //max # of bins in all directions

Chain_t* Chain;
Matrix_t* matrix;
gwParm_t* gwParm;

double L; //inital box size

double point[3][2], point0[3][2];
double L1[2], L2[2]; //vector of the box sides

double xmax,ymax,zmax; //max point in all directions
int tp;
uint32_t step;

Vector3D_t **distCOM;

int trajStep;

char *outputName, *trajName;
bool ReadInitFromFile;

#endif /* Header_h */

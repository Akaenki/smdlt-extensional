//
//  main.h
//  KRBC
//  
//  Created by Linling Miao on 5/1/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef main_h
#define main_h

#include "Parameters.h"

/* update the dimension of lattece every time step */
void updateLattice(int t);

/* Forces */
void initForce();
void forceBond();
void forceLJ();

/* Move beads */
void updateChain();
void locateD_ij(Vector3D_t COMd, float *D_ij);
void locateSelfD_ij(float *D_ij);

/* Center of Mass calculations */
Vector3D_t CenterOfMass(int i);//return CoM of Chain i

/* Move beads into the box */
void applyPBD();

/* apply PBC to Center of Masses */
Vector3D_t comPBC(Vector3D_t COM);

/* Calculate NID between two beads */
Vector3D_t getNID(int i, int j);

/* Calculate NID between two CoMs */
Vector3D_t CenterOfMassNID(Vector3D_t i, Vector3D_t j);

/* calculate and store the distance between CoMs */
void distCenterOfMass();

/* locate the bin of give distance of CoMs */
bin_t binDistCenterOfMass(Vector3D_t COMd);

float gasdev(long *idum);
void printTrajectory(int t);

/* Wall-clock timmer unit of nano seconds */
long long timer();
#endif /* main_h */

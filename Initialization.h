//
//  Initialization.h
//  KRBC
//
//  Created by Linling Miao on 5/1/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Initialization_h
#define Initialization_h

void Initialization(int argc, char * argv[]);

/* Function used to read intial traj from a file */
void readChains();

void generateOutput();

/* Read the binary iteration file */
void readIteration();

/* Read pre-averaged matrix file */
void readMatrixAvg();

/* Read gw parameters from file 
 * To use this move the gw calculation function to the end */
void readGeyerWinterParameters();

long initRan();
float ran1(long *idum);

/* Initialize the parallel environment */
void initParaEvir();

void initChains();

/* Initialize the Lattice */
void initLattice();

/* Allocate memories for the pointers */
void initMatrix();

void printInitial();

/* Will generate a binary file to store the iteration */
void printIteration();
                       
#endif /* Initialization_h */

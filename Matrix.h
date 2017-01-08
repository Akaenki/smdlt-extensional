//
//  Matrix.h
//  KRBC
//
//  Created by Linling Miao on 5/13/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h

#include "Parameters.h"

void updateMatrix(uint32_t step);

/* Calculate mobility matrices */
void regularRPY(int i,int j,bin_t bin);
void EwaldSumRPY(int i,int j,bin_t bin);

/* will create a .txt file to record values of beta_ij*/
void GeyerWinterCalculation();

void printMatrix();
/* Header:
    * Ntotal: total number of beads (uint32)*1
    * N: number of beads per chain (uint32)*1
    * binmax: max number of bins in one dimension (uint64)*1
    * dr: binsize (double)*1
 
 * Body:
    * DiffMatrix: full averaged matrix (float)*9*N*N*binmax**3
    * selfMatrix: self interaction matrix (float)*9*N*N
    * gridCounting: counter (uint32)*binmax**3
 */

#endif /* Matrix_h */

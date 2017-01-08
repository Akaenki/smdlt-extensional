//
//  Matrix.c
//  KRBC
//
//  Created by Linling Miao on 5/13/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "Matrix.h"
#include "main.h"

void updateMatrix(uint32_t step){
    bin_t bin;
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<NP; ++j){
            if(i==j){
                bin.x = 0; bin.y = 0; bin.z = 0;
                regularRPY(i,j,bin);
            } else{
                Vector3D_t COMd = distCOM[i][j];
                
                bin = binDistCenterOfMass(COMd);
                regularRPY(i,j,bin);
            }
        }
    }
    printf("%d\n",step);
}

void regularRPY(int i,int j,bin_t bin){
    int nn = 3*N;
    double mm1,mm2;
    
    if(i==j){
        for(int ii = 0; ii<N; ++ii){
            for(int jj = 0; jj<N; ++jj){
                int start0 = ii*9*N+jj*9;
                if(ii==jj){
                  matrix->selfMatrix[start0] += 1.0;
                  matrix->selfMatrix[start0+4] += 1.0;
                  matrix->selfMatrix[start0+8] += 1.0;
                } else{
                    Vector3D_t NID = getNID(i*N+ii, j*N+jj);
                    double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                    double r = sqrt(rr);
                    if(r>=2.0){
                        mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                        mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                    } else{
                        mm1 = 1.0-9.0*r/32.0;
                        mm2 = 3.0/(32.0*r);
                    }
                    matrix->selfMatrix[start0] += mm1+mm2*NID.x*NID.x;
                    matrix->selfMatrix[start0+1] += mm2*NID.x*NID.y;
                    matrix->selfMatrix[start0+2] += mm2*NID.x*NID.z;
                    matrix->selfMatrix[start0+3] += mm2*NID.y*NID.x;
                    matrix->selfMatrix[start0+4] += mm1+mm2*NID.y*NID.y;
                    matrix->selfMatrix[start0+5] += mm2*NID.y*NID.z;
                    matrix->selfMatrix[start0+6] += mm2*NID.z*NID.x;
                    matrix->selfMatrix[start0+7] += mm2*NID.z*NID.y;
                    matrix->selfMatrix[start0+8] += mm1+mm2*NID.z*NID.z;
                }
            }
        }
    } else{
        for(int ii = 0; ii<N; ++ii){
            for(int jj = 0; jj<N; ++jj){
                uint64_t start = bin.x*binymax*binzmax*nn*nn+bin.y*binzmax*nn*nn+bin.z*nn*nn+ii*9*N+jj*9;
                Vector3D_t NID = getNID(i*N+ii, j*N+jj);
                double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                double r = sqrt(rr);
                if(r>=2.0){
                    mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                    mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                } else{
                    mm1 = 1.0-9.0*r/32.0;
                    mm2 = 3.0/(32.0*r);
                }
                matrix->DiffMatrix[start] += mm1+mm2*NID.x*NID.x;
                matrix->DiffMatrix[start+1] += mm2*NID.x*NID.y;
                matrix->DiffMatrix[start+2] += mm2*NID.x*NID.z;
                matrix->DiffMatrix[start+3] += mm2*NID.y*NID.x;
                matrix->DiffMatrix[start+4] += mm1+mm2*NID.y*NID.y;
                matrix->DiffMatrix[start+5] += mm2*NID.y*NID.z;
                matrix->DiffMatrix[start+6] += mm2*NID.z*NID.x;
                matrix->DiffMatrix[start+7] += mm2*NID.z*NID.y;
                matrix->DiffMatrix[start+8] += mm1+mm2*NID.z*NID.z;
            }
        }
        matrix->gridCounting[bin.x*binymax*binzmax+bin.y*binzmax+bin.z] += 1;
    }
}


void GeyerWinterCalculation(){
    uint64_t bintotal = binxmax*binymax*binzmax;
    int nn = 3*N;
    double eps = 0.0;
    
    int i,j,k,l;
    for(i = 0; i<bintotal; ++i){
        //double sum = 0.0;
        for(j = 0; j<nn*nn; ++j){
            eps += matrix->DiffMatrixAvg[i*nn*nn+j];//sum += DiffMatrix[i*nn*nn+j];
        }
        //eps += sum*gridCounting[i]
        //printf("epsilon=%f\n", eps);
    }
    for(i = 0; i<nn*nn; ++i){
        eps += matrix->selfMatrixAvg[i];
    }
    eps -= nn;
    uint64_t nt = bintotal*nn;//row
    eps /= nt*nn - nn;
    printf("epsilon=%f\n", eps);
    
    gwParm->beta_ij = (1-sqrt(1-((nt-1)*eps*eps-(nt-2)*eps)))/((nt-1)*eps*eps-(nt-2)*eps);
    printf("beta=%f\n", gwParm->beta_ij);
    
    gwParm->C_i = calloc(nn, sizeof(float));
    for(i = 0; i<bintotal; ++i){
        for(j = 0; j<N; ++j){
            for(k = 0; k<N; ++k){
                for(l = 0; l<9; ++l){
                    gwParm->C_i[j*3+l/3] += matrix->DiffMatrixAvg[i*nn*nn+j*N*9+k*9+l]*matrix->DiffMatrixAvg[i*nn*nn+j*N*9+k*9+l];// gwParm->C_i[j] += DiffMatrix[i*nn*nn+j*nn+k]*DiffMatrix[i*nn*nn+j*nn+k]
                }
            }
        }
    }
    
    for(i = 0; i<N; ++i){
        for(j = 0; j<N; ++j){
            for(k = 0; k<9; ++k){
                gwParm->C_i[i*3+k/3] += matrix->selfMatrixAvg[i*9*N+j*9+k]*matrix->selfMatrixAvg[i*9*N+j*9+k];
            }
        }
    }
    
    for(i = 0; i<nn; ++i){
        gwParm->C_i[i] -= 1.0;
        gwParm->C_i[i] = gwParm->C_i[i]*gwParm->beta_ij*gwParm->beta_ij+1.0;
        gwParm->C_i[i] = 1.0/sqrt(gwParm->C_i[i]);
    }
}

void printMatrix(){
    uint64_t dtotal = 9*N*N*binxmax*binymax*binzmax;
    uint32_t stotal = 9*N*N;
    uint64_t btotal = binxmax*binymax*binzmax;
    
    FILE *Diffusion;
    char* str = malloc(100*sizeof(char));
    sprintf(str, "%s/Matrix_%s_%u.bin",directory,c_normalc,iteration);
    
    Diffusion = fopen(str,"wb");
    
    printf(">Generating Binary Averaged Matrix File...\n");
    
    uint32_t Ntotal = N*NP;
    double dr = BIN_SIZE;
    
    fwrite(&Ntotal,sizeof(uint32_t),1,Diffusion);
    fwrite(&N,sizeof(uint32_t),1,Diffusion);
    fwrite(&binxmax,sizeof(uint64_t),1,Diffusion);
    fwrite(&binymax,sizeof(uint64_t),1,Diffusion);
    fwrite(&binzmax,sizeof(uint64_t),1,Diffusion);
    fwrite(&dr,sizeof(double),1,Diffusion);
    
    float *diffAvg = calloc(dtotal,sizeof(float));
    float *selfAvg = calloc(stotal,sizeof(float));
    
    for(int i = 0; i<btotal; ++i){
        if(matrix->gridCounting[i] != 0){
            for(int j = 0; j<9*N*N; ++j){
                diffAvg[i*9*N*N+j] = matrix->DiffMatrix[i*9*N*N+j]/matrix->gridCounting[i];
            }
        }
    }
    
    for(int i = 0; i<stotal; ++i){
        selfAvg[i] = matrix->selfMatrix[i]/step/NP;
    }
    fwrite(diffAvg,sizeof(float),dtotal,Diffusion);
    fwrite(selfAvg,sizeof(float),stotal,Diffusion);
    fwrite(matrix->gridCounting,sizeof(uint32_t),btotal,Diffusion);
    fwrite(&step,sizeof(int),1,Diffusion);
    fclose(Diffusion);
    free(diffAvg);
    free(selfAvg);
    
    free(str);
}





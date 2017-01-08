//
//  main.c
//  KRBC
//
//  Created by Linling Miao on 4/27/16.
//  Copyright © 2017 Linling Miao. All rights reserved.
//

#include "main.h"
#include "Initialization.h"
#include "Matrix.h"

int main(int argc, char * argv[]) {
    
    Initialization(argc, argv);
    
    for(int t = 0; t<TMAX; t++){
        updateLattice(t);
        
        initForce();
        forceBond();
        forceLJ();
        
        if(t%100==0) updateMatrix(++step);
        
        if(t%trajStep==0) printTrajectory(t);
    }
    
    printMatrix();
    printIteration();
    
    if(iteration < MAXITER){
        char* command = malloc(200*sizeof(char));
        sprintf(command,"bash run.sh %s %lf %s 1",c_normalc,L,flowRatec);
        system(command);
        
        free(command);
    }
    
    return 0;
}

void updateLattice(int t){
    int treal = t%tp;
    //if(treal==0) printf("t=%d\n",t);
    
    for(int i = 0; i<3; ++i){
        point[i][0] = point0[i][0]*exp(flowRate*treal*DT); point[i][1] = point0[i][1]*exp(-flowRate*treal*DT);
    }
    for(int i = 0; i<2; ++i){
        L1[i] = point[2][i] - point[1][i]; L2[i] = point[0][i] - point[1][i];
    }
    //printf("x = %f y = %f\n", point[1][0], point[1][1]);
}

void initForce(){
    for(int i = 0; i<N*NP; ++i){
        Chain[i].fx = 0.0;
        Chain[i].fy = 0.0;
        Chain[i].fz = 0.0;
    }
}

void forceBond(){
    for(int i = 0; i<NP; ++i){
        for(int j = 1; j<N; ++j){
            Vector3D_t NID = getNID(i*N+j,i*N+j-1);
            double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
            double r = sqrt(rr);
            double Fs = -KAPPA*(r-2.0);
            double fx = Fs*NID.x/r; double fy = Fs*NID.y/r; double fz = Fs*NID.z/r;
            Chain[i*N+j].fx += fx;
            Chain[i*N+j].fy += fy;
            Chain[i*N+j].fz += fz;
            Chain[i*N+j-1].fx -= fx;
            Chain[i*N+j-1].fy -= fy;
            Chain[i*N+j-1].fz -= fz;
        }
    }
}

void forceLJ(){
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<N; ++j){
            for(int k = i*N+j+1; k<N*NP; ++k){
                Vector3D_t NID = getNID(i*N+j,k);
                double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
                if(rr<25){
                    double ratio = 4.00/rr;
                    double r6 = ratio*ratio*ratio;
                    if(r6>3) r6 = 3;
                    double coeff = (12*EPSILON/rr)*(r6*r6-r6);
                    Chain[i*N+j].fx += coeff*NID.x;
                    Chain[i*N+j].fy += coeff*NID.y;
                    Chain[i*N+j].fz += coeff*NID.z;
                    Chain[k].fx -= coeff*NID.x;
                    Chain[k].fy -= coeff*NID.y;
                    Chain[k].fz -= coeff*NID.z;
                }
            }
        }
    }
}

Vector3D_t CenterOfMass(int i){
    Vector3D_t COM;
    COM.x = Chain[i*N].rx;
    COM.y = Chain[i*N].ry;
    COM.z = Chain[i*N].rz;
    for(int j = 1; j<N; ++j){
        Vector3D_t NID = getNID(i*N+j,i*N+j-1);
        COM.x += (N-j)*NID.x/N;
        COM.y += (N-j)*NID.y/N;
        COM.z += (N-j)*NID.z/N;
    }
    COM = comPBC(COM);
    return COM;
}

Vector3D_t getNID(int i, int j){
    Vector3D_t NID;
    double theta = atan(L1[1]/L1[0]);
    double rx1 = Chain[i].rx; double rx2 = Chain[j].rx;
    double ry1 = Chain[i].ry; double ry2 = Chain[j].ry;
    double rz1 = Chain[i].rz; double rz2 = Chain[j].rz;
    //Transform
    double rx1p = rx1*cos(theta) + ry1*sin(theta); double ry1p = -rx1*sin(theta) + ry1*cos(theta);
    double rx2p = rx2*cos(theta) + ry2*sin(theta); double ry2p = -rx2*sin(theta) + ry2*cos(theta);
    double L1xp = L1[0]*cos(theta) + L1[1]*sin(theta);
    double L2xp = L2[0]*cos(theta) + L2[1]*sin(theta); double L2yp = -L2[0]*sin(theta) + L2[1]*cos(theta);
    double dx = rx1p - rx2p; double dy = ry1p - ry2p; double dz = rz1 - rz2;
    dx -= L2xp*round(dy/L2yp);
    dy -= L2yp*round(dy/L2yp);
    dx -= L1xp*round(dx/L1xp);
    dz -= L*round(dz/L);
    //Reverse Transform
    NID.x = dx*cos(theta) - dy*sin(theta);
    NID.y = dx*sin(theta) + dy*cos(theta);
    NID.z = dz;
    return NID;
}

Vector3D_t CenterOfMassNID(Vector3D_t i,Vector3D_t j){
    Vector3D_t COMd;
    double theta = atan(L1[1]/L1[0]);
    //Transform
    double rx1p = i.x*cos(theta) + i.y*sin(theta); double ry1p = -i.x*sin(theta) + i.y*cos(theta);
    double rx2p = j.x*cos(theta) + j.y*sin(theta); double ry2p = -j.x*sin(theta) + j.y*cos(theta);
    double L1xp = L1[0]*cos(theta) + L1[1]*sin(theta);
    double L2xp = L2[0]*cos(theta) + L2[1]*sin(theta); double L2yp = -L2[0]*sin(theta) + L2[1]*cos(theta);
    double dx = rx1p - rx2p; double dy = ry1p - ry2p; double dz = i.z - j.z;
    dx -= L2xp*round(dy/L2yp);
    dy -= L2yp*round(dy/L2yp);
    dx -= L1xp*round(dx/L1xp);
    dz -= L*round(dz/L);
    //Reverse Tranform;
    COMd.x = dx*cos(theta) - dy*sin(theta);
    COMd.y = dx*sin(theta) + dy*cos(theta);
    COMd.z = dz;
    return COMd;
}

void distCenterOfMass(){
    if(distCOM==NULL){
        distCOM = calloc(NP,sizeof(distCOM[0]));
        for(int i = 0; i<NP; ++i){
            distCOM[i] = calloc(NP,sizeof(distCOM[0][0]));
        }
    }
    
    for(int i = 0; i<NP; ++i){
        for(int j = 0; j<NP; ++j){
            Vector3D_t COM_i = CenterOfMass(i); Vector3D_t COM_j = CenterOfMass(j);
            distCOM[i][j] = CenterOfMassNID(COM_j, COM_i);
        }
    }
}

bin_t binDistCenterOfMass(Vector3D_t COMd){
    bin_t bin;
    double dr = BIN_SIZE;
    if(COMd.x>=0.0) bin.x = (int)(COMd.x/dr)+binxmax/2; else bin.x = floor(COMd.x/dr)+binxmax/2;
    if(COMd.y>=0.0) bin.y = (int)(COMd.y/dr)+binymax/2; else bin.y = floor(COMd.y/dr)+binymax/2;
    if(COMd.z>=0.0) bin.z = (int)(COMd.z/dr)+binzmax/2; else bin.z = floor(COMd.z/dr)+binzmax/2;
    return bin;
}

void applyPBC(){
    double theta = atan(L1[1]/L1[0]);
    //printf("l1x = %f l1y = %f\n",L1[0],L1[1]);
    //printf("%f\n",theta);
    //Transform
    double L1xp = L1[0]*cos(theta) + L1[1]*sin(theta);
    double L2xp = L2[0]*cos(theta) + L2[1]*sin(theta); double L2yp = -L2[0]*sin(theta) + L2[1]*cos(theta);
    for(int i = 0; i<N*NP; ++i){
        double rx = Chain[i].rx; double ry = Chain[i].ry;
        Chain[i].rx = rx*cos(theta) + ry*sin(theta);
        Chain[i].ry = -rx*sin(theta) + ry*cos(theta);
    }
    //PBC
    for(int i = 0; i<N*NP; ++i){
        Chain[i].rx -= L2xp*round(Chain[i].ry/L2yp);
        Chain[i].ry -= L2yp*round(Chain[i].ry/L2yp);
        Chain[i].rx -= L1xp*round((Chain[i].rx-Chain[i].ry*L2xp/L2yp)/L1xp);
        Chain[i].ry -= L*round(Chain[i].rz/L);
    }
    //Inverse Transform
    for(int i = 0; i<N*NP; ++i){
        double rx = Chain[i].rx; double ry = Chain[i].ry;
        Chain[i].rx = rx*cos(theta) - ry*sin(theta);
        Chain[i].ry = rx*sin(theta) + ry*cos(theta);
    }
    
}

Vector3D_t comPBC(Vector3D_t COM){
    Vector3D_t COM2;
    double theta = atan(L1[1]/L1[0]);
    double L1xp = L1[0]*cos(theta) + L1[1]*sin(theta);
    double L2xp = L2[0]*cos(theta) + L2[1]*sin(theta); double L2yp = -L2[0]*sin(theta) + L2[1]*cos(theta);
    //Transform
    COM2.x = COM.x*cos(theta) + COM.y*sin(theta);
    COM2.y = -COM.x*sin(theta) + COM.y*cos(theta);
    COM2.z = COM.z;
    //PBC
    COM2.x -= L2xp*round(COM2.x/L2yp);
    COM2.y -= L2yp*round(COM2.y/L2yp);
    COM2.x -= L1xp*round((COM2.x-COM2.y*L2xp/L2yp)/L1xp);
    COM2.z -= L*round(COM2.z/L);
    //Inverse Transform
    COM.x = COM2.x*cos(theta) - COM2.y*sin(theta);
    COM.y = COM2.x*sin(theta) + COM2.y*cos(theta);
    COM.z = COM2.z;
    
    return COM;
}

void updateChain(){
    double *RR = calloc(3*N*NP, sizeof(double));
    double p = sqrt(2*DT);
    for(int i = 0; i<3*N*NP; ++i){
        RR[i] = gasdev(idum);
    }
    int nn = 3*N;
    
    if(HImode){
#pragma omp parallel for schedule(static)
        for(int i = 0; i<NP; i++){
            float *D_noise = calloc(nn,sizeof(float));
            float *D_force = calloc(nn,sizeof(float));
            
            float *D_ij = calloc(nn*nn,sizeof(float));
            float *force = calloc(nn, sizeof(float));
            float *RRR = calloc(nn, sizeof(float));
            for(int j = 0; j<NP; ++j){
                if(i != j){
                    locateD_ij(distCOM[i][j], D_ij);
                } else{
                    locateSelfD_ij(D_ij);
                }
                
                for(int k = 0; k<N; ++k){
                    force[3*k] = Chain[j*N+k].fx;
                    force[3*k+1] = Chain[j*N+k].fy;
                    force[3*k+2] = Chain[j*N+k].fz;
                }
                
                for(int k = 0; k<nn; ++k){
                    RRR[k] = RR[j*nn+k];
                }
                
                cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,D_ij,nn,force,1,1,D_force,1);
                
                if(i == j){
                    for(int ii = 0; ii<nn; ++ii){
                        D_ij[ii*nn+ii] /= gwParm->beta_ij;
                    }
                }
                
                cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,gwParm->beta_ij,D_ij,nn,RRR,1,1,D_noise,1);
            }
            free(force);
            free(D_ij);
            free(RRR);
            
            //update
            for(int j = 0; j<N; ++j){
                Chain[i*N+j].rx += DT*(D_force[j*3]+flowRate*Chain[i*N+j].rx)+p*gwParm->C_i[j*3]*D_noise[j*3];
                Chain[i*N+j].ry += DT*(D_force[j*3+1]-flowRate*Chain[i*N+j].ry)+p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                Chain[i*N+j].rz += DT*D_force[j*3+2]+p*gwParm->C_i[j*3+2]*D_noise[j*3+1];
            }
            
            free(D_force);
            free(D_noise);
        }
        free(RR);
    } else{
#pragma omp parallel for schedule(static)
        for(int i = 0; i<NP*N; ++i){
            Chain[i].rx += DT*(Chain[i].fx+flowRate*Chain[i].rx)+p*RR[3*i];
            Chain[i].ry += DT*(Chain[i].fy-flowRate*Chain[i].ry)+p*RR[3*i+1];
            Chain[i].rz += DT*Chain[i].fz+p*RR[3*i+2];
        }
        free(RR);
    }
}

void locateD_ij(Vector3D_t COMd, float *D_ij){
    int nn = 3*N;
    bin_t bin = binDistCenterOfMass(COMd);
    uint64_t start = bin.x*binymax*binzmax*nn*nn+bin.y*binzmax*nn*nn+bin.z*nn*nn;
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            int start2 = i*9*N+j*9;
            for(int k = 0; k<9; ++k){
                D_ij[(i*3+k/3)*nn+j*3+k%3] = matrix->DiffMatrixAvg[start+start2+k];
            }
        }
        
    }
}

void locateSelfD_ij(float *D_ij){
    int nn = 3*N;
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            int start = i*9*N+j*9;
            for(int k = 0; k<9; ++k){
                D_ij[(i*3+k/3)*nn+j*3+k%3] = matrix->selfMatrixAvg[start+k];
            }
        }
    }
}

float gasdev(long *idum){
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0){
        do{
            v1 = 2.0*ran1(idum)-1.0;
            v2 = 2.0*ran1(idum)-1.0;
            rsq = v1*v1+v2*v2;
        }
        while(rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else{
        iset = 0;
        return gset;
    }
}

void printTrajectory(int t){
    char* str = malloc(sizeof(char)*30);
    sprintf(str, "KRBC_%d_%.3f.xyz", NP, flowRate);
    FILE *Trajectory;
    int isColor = 0;
    Trajectory = fopen(str, "a");
    fprintf(Trajectory, "%d\n%d\n", N*NP, t);
    for(int i = 0; i<N*NP; ++i){
        if(isColor){
            if(i>=10*N && i<20*N) fprintf(Trajectory, "O %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
            else fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
        } else{
            fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
        }
    }
    fclose(Trajectory);
    free(str);
}

long long timer(){
    struct timespec ts;
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return 1000000000*ts.tv_sec + ts.tv_nsec;
}
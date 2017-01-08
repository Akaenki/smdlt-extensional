//
//  Initialization.c
//  KRBC
//
//  Created by Linling Miao on 5/1/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "Initialization.h"
#include "Parameters.h"
#include "Matrix.h"
#include "util.h"

void Initialization(int argc, char * argv[]){
    ParseInput(argc, argv);
    
    matrix = (Matrix_t*)malloc(sizeof(Matrix_t));
    gwParm = (gwParm_t*)malloc(sizeof(gwParm_t));
    
    trajStep = 1000;
        
    directory = malloc(50*sizeof(char));
    directory = "./output";
    
    if(HImode){
        readIteration();
        readMatrixAvg();
        GeyerWinterCalculation();

        iteration++;
        
        //printIteration();
    } else{
        iteration = 1;
        
        //printIteration();
    }
    
    outputName = malloc(100*sizeof(char));
    sprintf(outputName,"%s/KRBC_%u.txt",directory,NP);
    trajName = malloc(100*sizeof(char));
    sprintf(trajName,"%s/KRBC_%u_%u.xyz",directory,NP,iteration);
    
    ReadInitFromFile = false;

    idum = malloc(sizeof(long));
    *idum = -1;
    ran1(idum);
    *idum = -1*initRan();

    //generateOutput();
    
    initParaEvir();
    
    initChains();
    
    if(ReadInitFromFile) readChains();
    
    initLattice();
    initMatrix();
    
    step = 0;
}

void generateOutput(){
    FILE *outputfile;
    outputfile = fopen(outputName, "w");
    fprintf(outputfile, "epsilon = %lf\n", EPSILON);
    fprintf(outputfile, "kappa = %lf\n", KAPPA);
    fprintf(outputfile, "N = %d\n", N);
    fprintf(outputfile, "dt = %lf\n", DT);
    fprintf(outputfile, "tmax = %d\n", TMAX);
    fprintf(outputfile, "NP = %d\n", NP);
    fprintf(outputfile, "L = %lf\n", L);
    fprintf(outputfile, "KMAX = %d\n", KMAX);
    fprintf(outputfile, "BinSize = %f\n", BIN_SIZE);
    fprintf(outputfile, "c_star = %lf\n\n",C_STAR);
    fprintf(outputfile, "SEED %ld\n",*idum);
    fclose(outputfile);
}

void readChains(){
    Chain = calloc(N*NP,sizeof(Chain_t));
    
    printf(">Reading Conformations from File...\n");
    
    FILE *trajec;
    trajec = fopen("./output/Trajectory.xyz","r");
    if(trajec==NULL){
        printf("!!Trajectory File not Founded, please use .xyz file. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i<N*NP; ++i){
        fscanf(trajec, "A %lf %lf %lf\n", &Chain[i].rx, &Chain[i].ry, &Chain[i].rz);
    }
    fclose(trajec);
    
    printf(">Read Conformation Sucessfully.\n");
}

void readIteration(){
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/iteration_%s.bin",directory,c_normalc);
    FILE *iter = fopen(str,"rb");
    if(iter==NULL){
        printf("!!Iteration File not Founded. Exiting...\n");
        exit(EXIT_FAILURE);
    } else{
        printf(">Reading Iteration...\n");
    }
    fread(&iteration,sizeof(uint16_t),1,iter);
    fclose(iter);
    free(str);
}

void initLattice(){
    tp = (int)(EP/flowRate/DT);
    printf("%d\n",tp);
    double theta = 90.0-(45.0+THETA0);
    point[0][0] = -sqrt(2)/2*L*cos(theta*M_PI/180); point[0][1] = sqrt(2)/2*L*sin(theta*M_PI/180);
    point[1][0] = -sqrt(2)/2*L*sin(theta*M_PI/180); point[1][1] = -sqrt(2)/2*L*cos(theta*M_PI/180);
    point[2][0] = sqrt(2)/2*L*cos(theta*M_PI/180); point[2][1] = -sqrt(2)/2*L*sin(theta*M_PI/180);
    for(int i = 0; i<3; ++i){
        for(int j = 0; j<2; ++j){
            point0[i][j] = point[i][j];
        }
    }
    for(int i = 0; i<2; ++i){
        L1[i] = point[2][i] - point[1][i]; L2[i] = point[0][i] - point[1][i];
    }
}

void readMatrixAvg(){
    FILE *Matrix;
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/Matrix_%s_%u.bin",directory,c_normalc,iteration);
    
    Matrix = fopen(str, "rb");
    if(Matrix==NULL){
        printf("!!Matrix File not Founded. Exiting...\n");
        exit(EXIT_FAILURE);
    }
    
    uint32_t Ntotal,Nc;
    uint64_t binsX,binsY,binsZ;
    double binsize;
    //header
    fread(&Ntotal,sizeof(uint32_t),1,Matrix);
    fread(&Nc,sizeof(uint32_t),1,Matrix);
    fread(&binsX,sizeof(uint64_t),1,Matrix);
    fread(&binsY,sizeof(uint64_t),1,Matrix);
    fread(&binsZ,sizeof(uint64_t),1,Matrix);
    fread(&binsize,sizeof(double),1,Matrix);
    //matrix
    uint64_t bintotal = binsX*binsY*binsZ;
    uint64_t mtotal = bintotal*Nc*Nc*9;
    int stotal = 9*Nc*Nc;
    
    printf(">total N:%d Nc:%d\n>total Bin:%llu total Elements:%llu\n",Ntotal,Nc,bintotal,mtotal);
    
    matrix->DiffMatrixAvg = (float*)malloc(mtotal*sizeof(float));
    matrix->selfMatrixAvg = (float*)malloc(stotal*sizeof(float));
    matrix->counter = (uint32_t*)malloc(bintotal*sizeof(uint32_t));
    
    //Original Structure
    fread(matrix->DiffMatrixAvg,sizeof(float),mtotal,Matrix);
    fread(matrix->selfMatrixAvg,sizeof(float),stotal,Matrix);
    fread(matrix->counter, sizeof(uint32_t),bintotal,Matrix);
    
    fclose(Matrix);
    free(str);
    
    printf(">Read Matrix Sucessfully.\n");
}

void readGeyerWinterParameters(){
    FILE* gwParameter;
    char *str = malloc(100*sizeof(char));
    sprintf(str,"%s/GeyerWinterParam_%d.bin",directory,iteration);
    gwParameter = fopen(str, "rb");
    
    gwParm->C_i = (float*)malloc(3*N*sizeof(float));
    fread(&gwParm->beta_ij,sizeof(double),1,gwParameter);
    fread(gwParm->C_i,sizeof(float),3*N,gwParameter);
    fclose(gwParameter);
    free(str);
}

void initParaEvir(){
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
    //numThreads = omp_get_num_procs();
    printf(">OpenMP Detected, Runing with %d threads\n",numThreads);
#else
    printf(">OpenMP not detected\n");
#endif
}

void initChains(){
    Chain = calloc(N*NP,sizeof(Chain_t));
    
    for(int i = 0; i<NP; ++i){
        int test = 0;
        while(test==0){
            test = 1;
            Chain[i*N].rx = ran1(idum)*L - L/2.0;
            Chain[i*N].ry = ran1(idum)*L - L/2.0;
            Chain[i*N].rz = ran1(idum)*L - L/2.0;
            
            Chain[i*N].rx -= round(Chain[i*N].rx/L)*L;
            Chain[i*N].ry -= round(Chain[i*N].ry/L)*L;
            Chain[i*N].rz -= round(Chain[i*N].rz/L)*L;
            
            for(int j = 0; j<i; ++j){
                double dx = Chain[i*N].rx - Chain[j*N].rx;
                double dy = Chain[i*N].ry - Chain[j*N].ry;
                double dz = Chain[i*N].rz - Chain[j*N].rz;
                dx -= round(dx/L)*L;
                dy -= round(dy/L)*L;
                dz -= round(dz/L)*L;
                
                if(dx*dx+dy*dy+dz*dz<2.0){
                    test = 0;
                }
            }
        }
    }
    
    for(int i = 0; i<NP; ++i){
        for(int j = 1; j<N; ++j){
            int test = 0;
            while(test==0){
                test = 1;
                double theta = ran1(idum)*2.0*3.14158;
                double phi = acos(2.0*ran1(idum)-1.0);
                Chain[i*N+j].rx = Chain[i*N+j-1].rx + 2.05*cos(theta)*sin(phi);
                Chain[i*N+j].ry = Chain[i*N+j-1].ry + 2.05*sin(theta)*sin(phi);
                Chain[i*N+j].rz = Chain[i*N+j-1].rz + 2.05*cos(phi);
                Chain[i*N+j].rx -= round(Chain[i*N+j].rx/L)*L;
                Chain[i*N+j].ry -= round(Chain[i*N+j].ry/L)*L;
                Chain[i*N+j].rz -= round(Chain[i*N+j].rz/L)*L;
                
                /*if(CheckOverlap){
                    int k;
                    for(k = 0; k<i*N+j; ++k){
                        double dx = Chain[i*N+j].rx-Chain[k].rx;
                        double dy = Chain[i*N+j].ry-Chain[k].ry;
                        double dz = Chain[i*N+j].rz-Chain[k].rz;
                        dx -= round(dx/L)*L;
                        dy -= round(dy/L)*L;
                        dz -= round(dz/L)*L;
                        
                        if(dx*dx+dy*dy+dz*dz<4.0){
                            test = 0;
                        }
                    }
                }*/
            }
        }
    }
}

void initMatrix(){
    xmax = point0[2][0]*exp(EP);
    ymax = -point0[1][1];
    zmax = L/2;
    //printf("%f %f\n",xmax,ymax);
    binxmax = 2*ceil(xmax/BIN_SIZE);
    binymax = 2*ceil(ymax/BIN_SIZE);
    binzmax = ceil(L/BIN_SIZE);
    int nn = 3*N;
    matrix->DiffMatrix = calloc(nn*nn*binxmax*binymax*binzmax,sizeof(float));
    matrix->selfMatrix = calloc(nn*nn,sizeof(float));
    matrix->gridCounting = calloc(binxmax*binymax*binzmax,sizeof(uint32_t));
}

void printInitial(){
    FILE *Trajectory;
    Trajectory = fopen(trajName, "w");
    fprintf(Trajectory, "%d\n-1\n", N*NP);
    for(int i = 0; i<N*NP; ++i){
        fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
    }
    fclose(Trajectory);
}

void printIteration(){
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/iteration_%s.bin",directory,c_normalc);
    FILE *iterat = fopen(str,"wb");
    printf(">Generating Iteration File...\n");
    fwrite(&iteration,sizeof(uint16_t),1,iterat);
    fclose(iterat);
    free(str);
}

long initRan(){
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000+1000000000;
    
    /*time_t seconds;
    time(&seconds); //switched from time due to crashes... unsure why...
    return -1*(unsigned long)(seconds); //Dividing by 100 keeps within correct bounds? not sure why this works */
}

float ran1(long *idum){
    long j;
    long k;
    static long idum2 = 123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if(*idum <= 0){
        if(-(*idum)<1) *idum=1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for(j=NTAB+7;j>=0;--j)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum<0) *idum+=IM1;
            if(j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum<0) *idum += IM1;
    k=idum2/IQ2;
    if(*idum<0) idum2+= IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if(iy<1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}





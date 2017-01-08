//
//  util.c
//  SingleChain
//
//  Created by Linling Miao on 12/18/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* this file only used to parse input arguments */

#include <getopt.h>
#include "util.h"

void ParseInput(int argc, char * argv[]){
    /* Default Values */
    HImode = 0;
    numThreads = 1;
    c_normal = 1.0;
    flowRate = 0;
    L = 20;
    N = 20;
    
    int option = 0;
    
    if(argc < 9) printf("Will use default values if an input is not specified\n");
    while((option = getopt(argc, argv, "n:t:m:l:c:f:")) != -1){
        switch(option){
            case 'n':
                sscanf(optarg, "%u", &N);
                break;
            case 't':
                sscanf(optarg, "%d", &numThreads);
                break;
            case 'l':
                sscanf(optarg, "%lf", &L);
                break;
            case 'c':
                c_normalc = optarg;
                sscanf(optarg, "%lf", &c_normal);
                break;
            case 'm':
                sscanf(optarg, "%d", &HImode);
                break;
            case 'f':
                flowRatec = optarg;
                sscanf(optarg, "%lf", &flowRate);
            case '?':
                printf("Unknown option -%c.\n Execution abort...", optopt);
                exit(EXIT_FAILURE);
        }
    }
    
    /* Calculate NP based on L and c_normal input */
    NP = ceil(L*L*L*(C_STAR*c_normal)/N);
    
    printf("Chain length: %d\n", N);
    printf("Box size: %lf, Concentration: %lf\n", L, c_normal);
    printf("Number of chains: %d\n", NP);
    printf("Flow rate: %lf\n", flowRate);
    //printf("Running with number of %d OpenMP threads...\n", numThreads);
}

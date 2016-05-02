// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"

int 
main(int argc, char **argv)
{
    signal(SIGSEGV, handler);
 	if(argc!=5){
 		fprintf(stderr, "Usage : ./ROUTE.exe -d=<num> -n=<num> <input_benchmark_name> <output_file_name> \n");
 		return 1;
 	}

 	int status;
    int netDecomp = atoi(argv[1]);
    int netOrdering = atoi(argv[2]);
	char *inputFileName = argv[3];
 	char *outputFileName = argv[4];

 	/// create a new routing instance
 	routingInst *rst = new routingInst;
	
 	/// read benchmark
    status = readBenchmark(inputFileName, rst);
 	if(status==0){
 		fprintf(stderr, "ERROR: reading input file \n");
 		return 1;
 	}

 	/// run actual routing
 	status = solveRouting(rst);
 	if(status==0){
 		fprintf(stderr, "ERROR: running routing \n");
 		release(rst);
 		return 1;
 	}

    if (netOrdering == 1){
        status = RRR(rst);
        if(status==0){
            fprintf(stderr, "ERROR: running rip-up and re-route\n");
            release(rst);
            return 1;
        }
    }

 	/// write the result
 	status = writeOutput(outputFileName, rst);
 	if(status==0){
 		fprintf(stderr, "ERROR: writing the result \n");
 		release(rst);
 		return 1;
 	}

 	release(rst);
    printf("\nDONE!\n\n");
 	return 0;
}

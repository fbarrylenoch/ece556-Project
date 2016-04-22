// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"

int main(int argc, char **argv)
{
    time_t init_time, curr_time;
    double seconds;

    signal(SIGSEGV, handler);
 	if(argc!=3){
 		printf("Usage : ./ROUTE.exe <input_benchmark_name> <output_file_name> \n");
 		return 1;
 	}
    time(&init_time);

 	int status;
	char *inputFileName = argv[1];
 	char *outputFileName = argv[2];

 	/// create a new routing instance
 	routingInst *rst = new routingInst;
	
 	/// read benchmark
    status = readBenchmark(inputFileName, rst);
 	if(status==0){
 		printf("ERROR: reading input file \n");
 		return 1;
 	}

    //printf("we passed reading in the benchmark\n");

    /// print benchmart
    //status = printRoutingInstince(rst);
	
 	/// run actual routing
 	status = solveRouting(rst);
 	if(status==0){
 		printf("ERROR: running routing \n");
 		release(rst);
 		return 1;
 	}

    
//    status = calcEdgeWeights();
    status = calcEdgeWeights(0, NULL, rst);
 	/// write the result
 	status = writeOutput(outputFileName, rst);
 	if(status==0){
 		printf("ERROR: writing the result \n");
 		release(rst);
 		return 1;
 	}

 	release(rst);
    //usleep(50000);
    //time(&curr_time);
    //seconds = difftime(curr_time, init_time);
    //printf("That took %.3f seconds", seconds);
 	printf("\nDONE!\n\n");
 	return 0;
}

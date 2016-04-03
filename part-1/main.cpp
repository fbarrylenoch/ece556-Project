// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"

int main(int argc, char **argv)
{
    signal(SIGSEGV, handler);
 	if(argc!=3){
 		printf("Usage : ./ROUTE.exe <input_benchmark_name> <output_file_name> \n");
 		return 1;
 	}

 	int status;
	char *inputFileName = argv[1];
 	char *outputFileName = argv[2];

 	/// create a new routing instance
 	routingInst *rst = new routingInst;
	
 	/// read benchmark
    printf("about to enter read Benchmark\n");
    status = readBenchmark(inputFileName, rst);
    for(int i = 0; i < rst->numEdges; i++)
        cout << "edge weight of e" << i << " = " << rst->edgeCaps[i] << "\t";
    cout << endl;
 	if(status==0){
 		printf("ERROR: reading input file \n");
 		return 1;
 	}
	
 	/// run actual routing
 	status = solveRouting(rst);
 	if(status==0){
 		printf("ERROR: running routing \n");
 		release(rst);
 		return 1;
 	}
	
 	/// write the result
 	status = writeOutput(outputFileName, rst);
 	if(status==0){
 		printf("ERROR: writing the result \n");
 		release(rst);
 		return 1;
 	}

 	release(rst);
 	printf("\nDONE!\n\n");
 	return 0;
}

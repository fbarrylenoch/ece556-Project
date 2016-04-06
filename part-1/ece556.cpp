#include "ece556.h"

void handler(int sig) {
    void *array[10];
    size_t size;

    size = backtrace(array, 10);

    fprintf(stderr, "ERROR: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

int readBenchmark(const char *fileName, routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/  
    ifstream fin;
    fin.open(fileName);
    if(!fin.good()) return 0;
    int v = 0;
    //read each line of the file
    const char* token[20] = {}; // initialize to 0
    while(!fin.eof()){
        v++;
        // read an entire line into memory
        char buf[512];
        fin.getline(buf, 512);

        // parse the line into blank-delimited tokens
        int n = 0; // a for-loop index

        // array to store memory addresses of the tokens in buff

        //parse the line
        token[0] = strtok(buf, " "); // subsequent tokens
        if(token[0]){ // zerof if line is blank
            for (n = 1; n < 20; n++){
                token[n] = strtok(0, " "); //subsequent tokens
                if (!token[n]) break; // no more tokens
            }
        }else break;



        if(strncmp(token[0], "grid", 64) == 0){
            int x = atoi(token[1]);
            int y = atoi(token[2]);

            rst->gx = x;
            rst->gy = y;
            rst->numEdges = y*(x-1) + x*(y-1);
            rst->edgeCaps = (int *)malloc(rst->numEdges * sizeof(int));
        }
        else if(strncmp(token[0], "capacity",64) == 0){
            rst->cap = atoi(token[1]);
            for (int i = 0; i < rst->numEdges; i++) 
                rst->edgeCaps[i] = rst->cap;
        }
        else if(strncmp(token[0], "num", 64) == 0){
            rst->numNets = atoi(token[2]);
            rst->nets = (net *)malloc(rst->numNets*sizeof(net));
        }
        else if(strcmp(token[0], "n0") == 0){
            for(int i = 0; i < rst->numNets; i++){
                int num = atoi(token[1]);
                net *tempNet = new net;
                tempNet->id = i;
                tempNet->numPins = num;
                tempNet->pins = (point *)malloc(num*sizeof(point));
                for(int j = 0; j < num; j++){
                    fin.getline(buf, 512);

                    // parse the line into blank-delimited tokens
                    int m = 0; // a for-loop index

                    //parse the line
                    token[0] = strtok(buf, " "); // subsequent tokens
                    if(token[0]){ // zerof if line is blank
                        for (m = 1; m < 20; m++){
                            token[m] = strtok(0, " "); //subsequent tokens
                            if (!token[m]) break; // no more tokens
                        }
                    }
                    point *tempPoint = new point;
                    tempPoint->x = atoi(token[0]);
                    tempPoint->y = atoi(token[1]);
                    tempNet->pins[j] = *tempPoint;
                    delete tempPoint;
                }
                rst->nets[i] = *tempNet;
                delete tempNet;
                if(i < rst->numNets -1) {// get the new net
                    fin.getline(buf, 512);

                    // parse the line into blank-delimited tokens
                    int m = 0; // a for-loop index

                    //parse the line
                    token[0] = strtok(buf, " "); // subsequent tokens
                    if(token[0]){ // zerof if line is blank
                        for (m = 1; m < 20; m++){
                            token[m] = strtok(0, " "); //subsequent tokens
                            if (!token[m]) break; // no more tokens
                        }
                    }
                }
            }
        }
        else{
            int num = atoi(token[0]);
            for(int i = 0; i < num; i++){
                fin.getline(buf, 512);
                // parse the line into blank-delimited tokens
                int m = 0; // a for-loop index

                //parse the line
                token[0] = strtok(buf, " "); // subsequent tokens
                if(token[0]){ // zerof if line is blank
                    for (m = 1; m < 20; m++){
                        token[m] = strtok(0, " "); //subsequent tokens
                        if (!token[m]) break; // no more tokens
                    }
                }
                int x1 = atoi(token[0]);
                int y1 = atoi(token[1]);
                int x2 = atoi(token[2]);
                int y2 = atoi(token[3]);
                int updatedCap = atoi(token[4]);

                // check if horizontal line
                if (y1 == y2) {
                    // calculate index in rst->edgeCaps
                    int xCap = 0;
                    if (x1 < x2)
                        xCap = x1;
                    else  
                        xCap = x2;

                    int index = 2*y1 + xCap;
                    rst->edgeCaps[index] = updatedCap;
                }
                // else vertical line
                else {
                    // calculate index in rst->edgeCaps
                    int yCap = 0;
                    if (y1 < y2)
                        yCap = y1;
                    else
                        yCap = y2;

                    // number of horizontal indexes
                    int horizontalIndexes = (rst->gy)*(rst->gx - 1);

                    int index = horizontalIndexes + 3*yCap + x1;
                    rst->edgeCaps[index] = updatedCap;
                }
            }
        }
    }
    //cout << "the first point of the last net = " << rst->nets[rst->numNets-1].pins[0].x << " " << rst->nets[rst->numNets-1].pins[0].y  << endl;
    return 1;
}

int solveRouting(routingInst *rst){
    /*********** TO BE FILLED BY YOU **********/
    //get every net in rst

    for (int i = 0; i < rst->numNets; i++) {
        net *tempNet = &(rst->nets[i]);
        tempNet->nroute.numSegs = 0;
        tempNet->nroute.segments = (segment *)malloc(tempNet->nroute.numSegs*sizeof(segment));

        for (int j = 0; j < tempNet->numPins-1; j++) {

            point p1 = tempNet->pins[j];
            point p2 = tempNet->pins[j+1];
            if ((p1.x == p2.x || p1.y == p2.y)) {
                segment tempSeg = tempNet->nroute.segments[tempNet->nroute.numSegs];
                tempNet->nroute.numSegs++;
		tempSeg.p1 = p1;
                tempSeg.p2 = p2;
                tempSeg.numEdges = 1;
                tempSeg.edges = (int *)malloc(4*sizeof(int));
                tempSeg.edges[0] = p1.x;
		tempSeg.edges[1] = p1.y;
		tempSeg.edges[2] = p2.x;
		tempSeg.edges[3] = p2.y;
            }
            else {
                point topP;
                point bottomP;
                if (p1.y < p2.y) {
                    topP = p1;
                    bottomP = p2;
                }
                else {
                    topP = p2;
                    bottomP = p1;
                 }
                // add topP to midPoint
                segment tempSeg = tempNet->nroute.segments[tempNet->nroute.numSegs];
		tempNet->nroute.numSegs++;
                tempSeg.p1 = topP;
                point *midPoint = new point;
                midPoint->x = topP.x;
                midPoint->y = bottomP.y;
                tempSeg.p2 = *midPoint;
                tempSeg.numEdges = 1;
                tempSeg.edges = (int *)malloc(4*sizeof(int));
                tempSeg.edges[0] = topP.x;
		tempSeg.edges[1] = topP.y;
		tempSeg.edges[2] = midPoint->x;
		tempSeg.edges[3] = midPoint->y;

                //add midPoint to bottomP
                tempSeg = tempNet->nroute.segments[tempNet->nroute.numSegs];
		tempNet->nroute.numSegs++;
                tempSeg.p1 = *midPoint;
                tempSeg.p2 = bottomP;
                tempSeg.numEdges = 1;
                tempSeg.edges = (int *)malloc(4*sizeof(int));
                tempSeg.edges[0] = midPoint->x;
                tempSeg.edges[1] = midPoint->y; 
                tempSeg.edges[2] = bottomP.x;
                tempSeg.edges[3] = bottomP.y;
            }
        }
    }
  return 1;
}

int writeOutput(const char *outRouteFile, routingInst *rst){
    /*********** TO BE FILLED BY YOU **********/
	try {
                ofstream output;
                output.open(outRouteFile);

                for (int i = 0; i < rst->numNets; i++) {

                        net tmpNet = rst->nets[i];
                        route tmpRoute = tmpNet.nroute;
                        output << "n" + tmpNet.id + "\n";

                        for (int j = 0; j < tmpRoute.numSegs; j++) {
                                segment tmpSeg = tmpRoute.segments[j];
                                output << "("+tmpSeg.p1.x+","+tmpSeg.p1.y+")-("+
                                        tmpSeg.p2.x+","+tmpSeg.p2.y+")\n";
                        }
                        output << "!\n";
                }
                output.close();
        } catch(int e) {
                return 0;
        }
    return 1;
}


int release(routingInst *rst){
    /*********** TO BE FILLED BY YOU **********/
	try {
                for (int i = 0; i < rst->numNets; i++) {
                        net deleteNet = rst->nets[i];
                        delete deleteNet.nroute.segments;
                        delete deleteNet.pins;
                }

                delete rst->edgeCaps;
                delete rst->edgeUtils;

                delete rst;

        } catch(int e) {
                return 0;
        }
    return 1;
}



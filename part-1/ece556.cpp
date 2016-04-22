#include "ece556.h"

int* getIndex(point p1, point p2, routingInst *rst){
    int *index;
    int distance;
    index = (int *)malloc(sizeof(int));
    index[0] = -1;

    // check if horizontal line
    if (p1.y == p2.y) {
        // calculate index in rst->edgeCaps
        if ((p1.x-p2.x) == -1) {
            //printf("\t\t1) distance = %d\n", distance);
            delete index;
            index = (int *)malloc(sizeof(int));
            index[0] = (rst->gx-1)*p1.y + p1.x;
        }
        else if((p1.x-p2.x) < -1) {
            distance = abs(p1.x-p2.x);
            //printf("\t\t2) distance = %d\n", distance);
            delete index;
            index = (int *)malloc(distance*sizeof(int));
            for (int i=0; i < distance; i++) {
                index[i] = (rst->gx-1)*p1.y + p1.x + i;
            }
        }            
        else if((p1.x-p2.x) == 1) {
            delete index;
            //printf("\t\t3) distance = %d\n", distance);
            index = (int *)malloc(sizeof(int));
            index[0] = (rst->gx-1)*p1.y + p2.x;
        }
        else {
            distance = abs(p1.x-p2.x);
            delete index;
            //pprintf("\t\t4) distance = %d\n", distance);
            index = (int *)malloc(distance*sizeof(int));
            for (int i=0; i < distance; i++) {
                index[i] = (rst->gx-1)*p1.y + p2.x + i;
            }
        } // case where p1.x-p2.x > 1
    }
    // else vertical line
    else if(p1.x == p2.x){
        // calculate index in rst->edgeCaps
        if ((p1.y-p2.y) == -1) {
            delete index;
            //printf("\t\t5) distance = %d\n", distance);
            index = (int *)malloc(sizeof(int));
            index[0] = (rst->gy)*(rst->gx-1) + rst->gy*p1.y + p1.x;
        }
        else if ((p1.y-p2.y) < -1) {
            distance = abs(p1.y-p2.y);
            //printf("\t\t6) distance = %d\n", distance);
            delete index;
            index = (int *)malloc(distance*sizeof(int));
            for (int i=0; i < distance; i++) {
                index[i] = (rst->gy)*(rst->gx-1) + rst->gy*(p1.y + i) + p1.x;
            }
        }
        else if ((p1.y-p2.y) == 1) {
            delete index;
            //printf("\t\t7) distance = %d\n", distance);
            index = (int *)malloc(sizeof(int));
            index[0] = (rst->gy)*(rst->gx-1) + rst->gy*p2.y + p1.x;
        }
        else {
            distance = abs(p1.y-p2.y);
            //printf("\t\t8) distance = %d\n", distance);
            delete index;
            index = (int *)malloc(distance*sizeof(int));
            for (int i=0; i < distance; i++) {
                index[i] = (rst->gy)*(rst->gx-1) + rst->gy*(p2.y + i) + p1.x;
            }
        } // case where p1.y-p2.y > 1
    }
    return index;
}

//will be called to rip up and reroute a net
bool RRR(routingInst *rst, net *netRRR){
    int* indices;
    point p1;
    point p2;
    route routeRRR;
    segment *segRRR;
    // routes for each shape attempt - not sure if I will use these
    route routeL; /* used to store the best L route option */
    route routeZ; /* used to store the best Z route option */
    route routeRZ; /* used to store the best Rotated Z route option */
    route routeU; /* used to store the best U route option */
    route routeC; /* used to store the best C route option */

    routeRRR = netRRR->nroute;
    // rip up
    for (int i = 0; i < routeRRR.numSegs; i++){ /* get each segment in the net*/
        segRRR = &(routeRRR.segments[i]);

        // assign endpoints of the segment
        p1 = segRRR->p1;
        p2 = segRRR->p2;

        indices = getIndex(p1,p2,rst);  /* get the edge indices of the segment*/

        // update utilizations and weights for each edge
        for (int i = 0; i < segRRR->numEdges; i++){
            rst->edgeCaps[segRRR->edges[i]]++;
            rst->edgeUtils[segRRR->edges[i]]--;
        }
    }

    // clear the segments of the route
    routeRRR.numSegs = 0;
    delete routeRRR.segments;

    // reroute
    /* allocate enough memory for every connection from pin to pin to be the largest number of segments */
    routeRRR.segments = (segment *)malloc((netRRR->numPins-1)*3*sizeof(segment));
     
        // try L shape - maybe seperate method?

        // try Z shape - maybe seperate method?

        // try Rotated Z shape - maybe seperate method?

        // try U shape - maybe seperate method?

        // try Rotated U shape - maybe seperate method?

        // try C shape - maybe seperate method?

        // try Rotated C shape - maybe seperate method?

        // compare each attempt, take best option


    // if error for any reason, return false
    /* return false; */

    // else return successful
    return true;
}

// returns best L shape route - return straight line if in line
route shapeL(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeL; /* used to store the best L route option */
    point topP; /* used to calculate route */
    point bottomP; /* used to calculate route */
    point midPoint; /* used to calculate route */
    
    /* check if straight line, assign routeL if so */
    if (p1.x == p2.x || p1.y == p2.y) {
        routeL.segments = (segment *)malloc(sizeof(segment));
        routeL.numSegs = 0;

        segment *segL = &(routeL.segments[routeL.numSegs]);
        routeL.numSegs++;
        segL->p1 = p1;
        segL->p2 = p2;
        if (p1.x == p2.x) {
            segL->numEdges = abs(p1.y - p2.y);
        }
        else {
            segL->numEdges = abs(p1.x - p2.x);
        }

        int* edgesL = getIndex(p1,p2,rst);
        segL->edges = (int *)malloc(segL->numEdges*sizeof(int));
        segL->edges = edgesL;
    }
    else{

        /* assign L shape route */
        if (p1.y < p2.y) {
            topP = p2;
            bottomP = p1;
        }
        else {
            topP = p1;
            bottomP = p2;
        }
        // allocate routeL
        routeL.segments = (segment *)malloc(2*sizeof(segment));
        routeL.numSegs = 0;
        // assign top and bottom points
        segment *segL = &(routeL.segments[routeL.numSegs]);
        routeL.numSegs++;
        segL->p1 = topP;
        // create a midPoint
        midPoint.x = topP.x;
        midPoint.y = bottomP.y;
        segL->p2 = midPoint;
        segL->numEdges = abs(topP.y - midPoint.y);
        // save edges to the segment
        int* edgesL = getIndex(topP,midPoint,rst);
        segL->edges = (int *)malloc(segL->numEdges*sizeof(int));
        segL->edges = edgesL;
        // add midPoint to bottomP
        segL = &(routeL.segments[routeL.numSegs]);
        routeL.numSegs++;
        segL->p1 = midPoint;
        segL->p2 = bottomP;
        segL->numEdges = abs(midPoint.x - bottomP.x);
        // save edges to the segment
        edgesL = getIndex(midPoint,bottomP,rst);
        segL->edges = (int *)malloc(segL->numEdges*sizeof(int));
        segL->edges = edgesL;

        /* assign L in other direction */
        if (p1.y < p2.y) {
            topP = p2;
            bottomP = p1;
        }
        else {
            topP = p1;
            bottomP = p2;
        }
        // allocate routeL
        routeComp.segments = (segment *)malloc(2*sizeof(segment));
        routeComp.numSegs = 0;
        // assign top and bottom points
        segment *segComp = &(routeComp.segments[routeComp.numSegs]);
        routeComp.numSegs++;
        segComp->p1 = topP;
        // create a midPoint
        midPoint.y = topP.y;
        midPoint.x = bottomP.x;
        segComp->p2 = midPoint;
        segComp->numEdges = abs(topP.x - midPoint.x);
        // save edges to the segment
        int* edgesComp = getIndex(topP,midPoint,rst);
        segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
        segComp->edges = edgesComp;
        // add midPoint to bottomP
        segComp = &(routeComp.segments[routeComp.numSegs]);
        routeComp.numSegs++;
        segComp->p1 = midPoint;
        segComp->p2 = bottomP;
        segComp->numEdges = abs(midPoint.y - bottomP.y);
        // save edges to the segment
        edgesComp = getIndex(midPoint,bottomP,rst);
        segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
        segComp->edges = edgesComp;

        /* compare and take route with lowest weight */
        // need edge weight to be completed before this can be completed
    }
    // return the route
    return routeL;
}

// returns best Z shape route - call shapeL if too close to make a Z
route shapeZ(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeZ; /* used to store the best Z route option */
    point topPoint; /* used to calculate route */
    point bottomPoint; /* used to calculate route */
    point topMidPoint; /* used to calculate route */
    point bottomMidPoint; /* used to calculate route */
    int multiplier; /* used to get proper dirction for route generation */

    // call shapeL if too close
    if (abs(p1.x - p2.x) < 2) {
        routeZ = shapeL(p1,p2,rst);
    }

    // else assign Z shape route
    else {
        // allocate space for routeZ and routeComp
        routeZ.segments = (segment *)malloc(3*sizeof(segment));
        routeZ.numSegs = 0;
        routeComp.segments = (segment *)malloc(3*sizeof(segment));
        routeComp.numSegs = 0;
        /* assign outside most Z route */
        // assign top point
        if (p1.y < p2.y) {
            topPoint = p2;
            bottomPoint = p1;
        }
        else {
            topPoint = p1;
            bottomPoint = p2;
        }
        // assign multiplier
        if (topPoint.x < bottomPoint.x) {
            multiplier = 1;
        }
        else {
            multiplier = -1;
        }
        // create midpoints
        topMidPoint.x = (topPoint.x + (1 * multiplier));
        topMidPoint.y = topPoint.y;
        bottomMidPoint.x = (topPoint.x + (1 * multiplier));
        bottomMidPoint.y = bottomPoint.y;
        // create top segment
        segment *segZ = &(routeZ.segments[routeZ.numSegs]);
        routeZ.numSegs++;
        segZ->p1 = topPoint;
        segZ->p2 = topMidPoint;
        segZ->numEdges = abs(topPoint.x - topMidPoint.x);
        // get top edges
        int* edgesZ = getIndex(topPoint,topMidPoint,rst);
        segZ->edges = (int *)malloc(segZ->numEdges*sizeof(int));
        segZ->edges = edgesZ;
        // create mid segment
        segZ = &(routeZ.segments[routeZ.numSegs]);
        routeZ.numSegs++;
        segZ->p1 = topMidPoint;
        segZ->p2 = bottomMidPoint;
        segZ->numEdges = abs(topMidPoint.y - bottomMidPoint.y);
        // get mid edges
        edgesZ = getIndex(topMidPoint,bottomMidPoint,rst);
        segZ->edges = (int *)malloc(segZ->numEdges*sizeof(int));
        segZ->edges = edgesZ;
        // create bottom segment
        segZ = &(routeZ.segments[routeZ.numSegs]);
        routeZ.numSegs++;
        segZ->p1 = bottomMidPoint;
        segZ->p2 = bottomPoint;
        segZ->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
        // get bottom edges
        edgesZ = getIndex(bottomMidPoint,bottomPoint,rst);
        segZ->edges = (int *)malloc(segZ->numEdges*sizeof(int));
        segZ->edges = edgesZ;

        //check every other possible Z route and compare weights
        for (int i = 1; i < (abs(p1.x - p2.x)-1); i++){
            routeComp.numSegs = 0;
            // create midpoints
            topMidPoint.x = (topPoint.x + ((1 + i)*multiplier));
            topMidPoint.y = topPoint.y;
            bottomMidPoint.x = (topPoint.x + ((1 + i)*multiplier));
            bottomMidPoint.y = bottomPoint.y;
            // create top segment
            segment *segComp = &(routeComp.segments[routeComp.numSegs]);
            routeComp.numSegs++;
            segComp->p1 = topPoint;
            segComp->p2 = topMidPoint;
            segComp->numEdges = abs(topPoint.x - topMidPoint.x);
            // get top edges
            int* edgesComp = getIndex(topPoint,topMidPoint,rst);
            segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
            segComp->edges = edgesComp;
            // create mid segment
            segComp = &(routeComp.segments[routeComp.numSegs]);
            routeComp.numSegs++;
            segComp->p1 = topMidPoint;
            segComp->p2 = bottomMidPoint;
            segComp->numEdges = abs(topMidPoint.y - bottomMidPoint.y);
            // get mid edges
            edgesComp = getIndex(topMidPoint,bottomMidPoint,rst);
            segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
            segComp->edges = edgesComp;
            // create bottom segment
            segComp = &(routeComp.segments[routeComp.numSegs]);
            routeComp.numSegs++;
            segComp->p1 = bottomMidPoint;
            segComp->p2 = bottomPoint;
            segComp->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
            // get bottom edges
            edgesComp = getIndex(bottomMidPoint,bottomPoint,rst);
            segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
            segComp->edges = edgesComp;

            // check if weight is lower than routeZ, assign to routeZ if true
            /* 
            if (routeComp_weight < routeZ_weight){
                routeZ = routeComp;
            }
            */
        }
    }
    // return the route
    return routeZ;
}

// returns best Rotated Z shape route - call shapeL if too close to make a RZ
route shapeRZ(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeRZ; /* used to store the best Z route option */
    point leftMidPoint; /* used to calculate route */
    point rightMidPoint; /* used to calculate route */
    point leftPoint; /* used to calculate route */
    point rightPoint; /* used to calculate route */
    int multiplier; /* used to get proper direction for route generation */


    // call shapeL if too close
    if (abs(p1.y - p2.y) < 2){
        routeRZ = shapeL(p1,p2,rst);
    }

    // else loop for all possible RZ routes
    else {
        // allocate space for routeRZ and routeComp
        routeRZ.segments = (segment *)malloc(3*sizeof(segment));
        routeRZ.numSegs = 0;
        routeComp.segments = (segment *)malloc(3*sizeof(segment));
        routeComp.numSegs = 0;
        // assign farthest left oriented Z route
        if (p1.x < p2.x) {
            rightPoint = p2;
            leftPoint = p1;
        }
        else {
            rightPoint = p1;
            leftPoint = p2;
        }
        // assign multiplier
        if (rightPoint.y < leftPoint.y) {
            multiplier = 1;
        }
        else {
            multiplier = -1;
        }
        // create midpoints
        leftMidPoint.x = leftPoint.x;
        leftMidPoint.y = leftPoint.y - (1 * multiplier);
        rightMidPoint.x = rightPoint.x;
        rightMidPoint.y = leftPoint.y - (1 * multiplier);
        // create left segment
        segment *segRZ = &(routeRZ.segments[routeRZ.numSegs]);
        routeRZ.numSegs++;
        segRZ->p1 = leftPoint;
        segRZ->p2 = leftMidPoint;
        segRZ->numEdges = abs(leftPoint.y - leftMidPoint.y);
        // get left edges
        int* edgesRZ = getIndex(leftPoint,leftMidPoint,rst);
        segRZ->edges = (int *)malloc(segRZ->numEdges*sizeof(int));
        segRZ->edges = edgesRZ;
        // create mid segment
        segRZ = &(routeRZ.segments[routeRZ.numSegs]);
        routeRZ.numSegs++;
        segRZ->p1 = leftMidPoint;
        segRZ->p2 = rightMidPoint;
        segRZ->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
        // get mid edges
        edgesRZ = getIndex(leftMidPoint,rightMidPoint,rst);
        segRZ->edges = (int *)malloc(segRZ->numEdges*sizeof(int));
        segRZ->edges = edgesRZ;
        // create right segment
        segRZ = &(routeRZ.segments[routeRZ.numSegs]);
        routeRZ.numSegs++;
        segRZ->p1 = rightMidPoint;
        segRZ->p2 = rightPoint;
        segRZ->numEdges = abs(rightMidPoint.y - rightPoint.y);
        // get right edges
        edgesRZ = getIndex(rightMidPoint,rightPoint,rst);
        segRZ->edges = (int *)malloc(segRZ->numEdges*sizeof(int));
        segRZ->edges = edgesRZ;

        //check every other possible RZ route and compare weights
        for (int i = 1; i < (abs(p1.y - p2.y)-1); i++){
            routeComp.numSegs = 0;
            // create midpoints
            leftMidPoint.x = leftPoint.x;
            leftMidPoint.y = leftPoint.y - ((1 + i) * multiplier);
            rightMidPoint.x = rightPoint.x;
            rightMidPoint.y = leftPoint.y - ((1 + i) * multiplier);
            // create left segment
            segment *segComp = &(routeComp.segments[routeComp.numSegs]);
            routeComp.numSegs++;
            segComp->p1 = leftPoint;
            segComp->p2 = leftMidPoint;
            segComp->numEdges = abs(leftPoint.y - leftMidPoint.y);
            // get left edges
            int* edgesComp = getIndex(leftPoint,leftMidPoint,rst);
            segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
            segComp->edges = edgesComp;
            // create mid segment
            segComp = &(routeComp.segments[routeComp.numSegs]);
            routeComp.numSegs++;
            segComp->p1 = leftMidPoint;
            segComp->p2 = rightMidPoint;
            segComp->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
            // get mid edges
            edgesComp = getIndex(leftMidPoint,rightMidPoint,rst);
            segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
            segComp->edges = edgesComp;
            // create right segment
            segComp = &(routeComp.segments[routeComp.numSegs]);
            routeComp.numSegs++;
            segComp->p1 = rightMidPoint;
            segComp->p2 = rightPoint;
            segComp->numEdges = abs(rightMidPoint.y - rightPoint.y);
            // get right edges
            edgesComp = getIndex(rightMidPoint,rightPoint,rst);
            segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
            segComp->edges = edgesComp;

            // check if weight is lower than routeRZ, assign to routeRZ if true
            /* 
            if (routeComp_weight < routeRZ_weight){
                routeRZ = routeComp;
            }
            */
        }
    }
    // return the route
    return routeRZ;
}

// returns best U shape route
route shapeU(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeU; /* used to store the best U route option */
    
    // return the route
    return routeU;
}

// returns the best Rotated U shape route
route shapeRU(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeRU; /* used to store the best U route option */
    
    // return the route
    return routeRU;
}

// returns best C shape route
route shapeC(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeC; /* used to store the best C route option */
    
    // return the route
    return routeC;
}

// returns best Rotated C shape route
route shapeRC(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeRC; /* used to store the best C route option */
    
    // return the route
    return routeRC;
}

void handler(int sig) {
    void *array[10];
    size_t size;

    size = backtrace(array, 10);

    fprintf(stderr, "ERROR: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

int readBenchmark(const char *fileName, routingInst *rst){
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
        token[0] = strtok(buf, "\t "); // subsequent tokens
        if(token[0]){ // zerof if line is blank
            for (n = 1; n < 20; n++){
                token[n] = strtok(0, "\t "); //subsequent tokens
                if (!token[n]) break; // no more tokens
            }
        }else break;

        // check for grid size
        if(strncmp(token[0], "grid", 64) == 0){
            //printf("check for grid size\n");
            int x = atoi(token[1]);
            int y = atoi(token[2]);

            rst->gx = x;
            rst->gy = y;
            rst->numEdges = y*(x-1) + x*(y-1);
            rst->edgeCaps = (int *)malloc(rst->numEdges * sizeof(int));
            rst->edgeUtils = (int *)malloc(rst->numEdges * sizeof(int));
            rst->edgeWeight = (int *) malloc(rst->numEdges * sizeof(int));
            rst->edgeHis = (int *) malloc(rst->numEdges * sizeof(int));
            // going to set edgeUtils to 0 as a default and 1 if the edge is utilized
            for (int i = 0; i < rst->numEdges; i++) {
                rst->edgeUtils[i] = 0;
            }
        }
        // check for edge caps
        else if(strncmp(token[0], "capacity",64) == 0){
            //printf("check for edge caps\n");
            rst->cap = atoi(token[1]);
            for (int i = 0; i < rst->numEdges; i++) 
                rst->edgeCaps[i] = rst->cap;
        }
        // check for number of nets
        else if(strncmp(token[0], "num", 64) == 0){
            //printf("check for number of nets\n");
            rst->numNets = atoi(token[2]);
            rst->nets = (net *)malloc(rst->numNets*sizeof(net));
        }
        // check for nets
        else if(strcmp(token[0], "n0") == 0){
            //printf("we are reading in nets\n");
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
                    token[0] = strtok(buf, "\t "); // subsequent tokens
                    if(token[0]){ // zerof if line is blank
                        for (m = 1; m < 20; m++){
                            token[m] = strtok(0, "\t "); //subsequent tokens
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
                    token[0] = strtok(buf, "\t "); // subsequent tokens
                    if(token[0]){ // zerof if line is blank
                        for (m = 1; m < 20; m++){
                            token[m] = strtok(0, "\t "); //subsequent tokens
                            if (!token[m]) break; // no more tokens
                        }
                    }
                }
            }
        }
        else{
            int num = atoi(token[0]);
            //printf("starting blockage constraints, there are %d blockages\n", num);
            for(int i = 0; i < num; i++){
                fin.getline(buf, 512);
                // parse the line into blank-delimited tokens
                token[0] = strtok(buf, "\t "); // subsequent tokens
                if(token[0]){ // zerof if line is blank
                    for (int m = 1; m < 20; m++){
                        token[m] = strtok(0, "\t "); //subsequent tokens
                        if (!token[m]) break; // no more tokens
                    }
                }
                int *index;
                point* p1 = new point;
                point* p2 = new point;
                p1->x = atoi(token[0]);
                p1->y = atoi(token[1]);
                p2->x = atoi(token[2]);
                p2->y = atoi(token[3]);
                int updatedCap = atoi(token[4]);
                //printf("\tblockage: %d\n", i);
                //printf("\tblockage (%d/%d) between (%d,%d)->(%d,%d) with cap %d\n",
                //        i, num, p1->x, p1->y, p2->x, p2->y, updatedCap);
                index = getIndex(*p1, *p2, rst);
                rst->edgeCaps[index[0]] = updatedCap;
                delete p1;
                delete p2;
            }
        }
    }
    return 1;
}

int solveRouting(routingInst *rst){
    //get every net in rst
    int TOF = 0; // temporarily being used to calculate total overflow
    int TWL = 0; // temporarily being used to calculate total wire length
    for (int i = 0; i < rst->numNets; i++) {
        net *tempNet = &(rst->nets[i]);
        tempNet->nroute.numSegs = 0;

        tempNet->nroute.segments = (segment *)malloc((tempNet->numPins-1)*2*sizeof(segment));

        for (int j = 0; j < tempNet->numPins-1; j++) {

            point p1 = tempNet->pins[j];
            point p2 = tempNet->pins[j+1];
            if (p1.x == p2.x || p1.y == p2.y) {
                segment *tempSeg = &(tempNet->nroute.segments[tempNet->nroute.numSegs]);
                tempNet->nroute.numSegs++;
                tempSeg->p1 = p1;
                tempSeg->p2 = p2;
                if (p1.x == p2.x) {
                    tempSeg->numEdges = abs(p1.y - p2.y);
                }
                else {
                    tempSeg->numEdges = abs(p1.x - p2.x);
                }

                int* tempEdges = getIndex(p1,p2,rst);
                tempSeg->edges = (int *)malloc(tempSeg->numEdges*sizeof(int));
                tempSeg->edges = tempEdges;

                for (int i = 0; i < tempSeg->numEdges; i++){
                    rst->edgeCaps[tempSeg->edges[i]]--;
                    rst->edgeUtils[tempSeg->edges[i]]++;
                }

            }
            else {
                point topP;
                point bottomP;
                if (p1.y < p2.y) {
                    topP = p2;
                    bottomP = p1;
                }
                else {
                    topP = p1;
                    bottomP = p2;
                }

                // add topP to midPoint
                segment *tempSeg = &(tempNet->nroute.segments[tempNet->nroute.numSegs]);
                tempNet->nroute.numSegs++;
                tempSeg->p1 = topP;
                point midPoint;
                midPoint.x = topP.x;
                midPoint.y = bottomP.y;
                tempSeg->p2 = midPoint;
                tempSeg->numEdges = abs(topP.y - midPoint.y);

                // save edges to the segment
                int* tempEdges = getIndex(topP,midPoint,rst);
                tempSeg->edges = (int *)malloc(tempSeg->numEdges*sizeof(int));
                tempSeg->edges = tempEdges;

                for (int i = 0; i < tempSeg->numEdges; i++){
                    rst->edgeCaps[tempSeg->edges[i]]--;
                    rst->edgeUtils[tempSeg->edges[i]]++;
                }

                // uncomment for test purposes
                /*cout << "Start Segment:\n";
                for (int i = 0; i < tempSeg->numEdges; i++) {
                    cout << "edge: " << tempSeg->edges[i] << "\n";
                }
                cout << "\n";*/

                // add midPoint to bottomP
                tempSeg = &(tempNet->nroute.segments[tempNet->nroute.numSegs]);
                tempNet->nroute.numSegs++;
                tempSeg->p1 = midPoint;
                tempSeg->p2 = bottomP;
                tempSeg->numEdges = abs(midPoint.x - bottomP.x);

                // save edges to the segment
                tempEdges = getIndex(midPoint,bottomP,rst);
                tempSeg->edges = (int *)malloc(tempSeg->numEdges*sizeof(int));
                tempSeg->edges = tempEdges;

                for (int i = 0; i < tempSeg->numEdges; i++){
                    rst->edgeCaps[tempSeg->edges[i]]--;
                    rst->edgeUtils[tempSeg->edges[i]]++;
                }

                // uncomment for test purposes
                /*cout << "Start Segment:\n";
                for (int i = 0; i < tempSeg->numEdges; i++) {
                    cout << "edge: " << tempSeg->edges[i] << "\n";
                }
                cout << "\n";*/
            }
        }
    }

    // I am like 80% sure this works...
    for (int i = 0; i < rst->numEdges; i++){
        if (rst->edgeCaps[i] < 0)
            TOF = TOF + abs(rst->edgeCaps[i]);
        if (rst->edgeUtils[i] > 0)
            TWL = TWL + rst->edgeUtils[i];
    }
    cout << "Total Overflow: " << TOF << "\n";
    cout << "Total Wire Length: " << TWL << "\n";
    return 1;
}

int calcEdgeWeights(int mode, route *newRoute, routingInst *rst){

    int Overflow;
    // calculate all edge weights at the beginning
    if(mode == 0){
        for(int i = 0; i < rst->numEdges; ++i){
            Overflow = max((rst->edgeUtils[i] - rst->edgeCaps[i]), 0);
            rst->edgeHis[i] = 0;
            if(Overflow > 0)
                rst->edgeHis[i] += 1;
            rst->edgeWeight[i] = Overflow * rst->edgeHis[i];
            //printf("%d, ", rst->edgeWeight[i]);
        }

        cout << endl;
        return 1;
    }

    else if(mode == 1){

        return 1;
    }
    else if(mode == 2){

        return 1;
    }
    return 0;
}

int writeOutput(const char *outRouteFile, routingInst *rst){
	try {
        ofstream output;
        output.open(outRouteFile);

        for (int i = 0; i < rst->numNets; i++) {

            net tmpNet = rst->nets[i];
            route tmpRoute = tmpNet.nroute;
            output << "n" << tmpNet.id << "\n";

            for (int j = 0; j < tmpRoute.numSegs; j++) {
                segment tmpSeg = tmpRoute.segments[j];
                output << "(" << tmpSeg.p1.x << "," << tmpSeg.p1.y << ")-(" << 
                    tmpSeg.p2.x << "," << tmpSeg.p2.y << ")\n";
            }
            output << "!\n";
        }
        output.close();
    } catch(int e) {
        return 0;
    }
    return 1;
}

int printRoutingInstince (routingInst *rst){
    for(int i = 0; i < rst->numNets; ++i){
        printf("net n%d, has %d pins\n", i, rst->nets[i].numPins);
       /* printf("the pin connections are:\n");
        for(int j = 0; j < rst->nets[i].numPins; ++j){
            printf("\t(%d,%d)\n", rst->nets[i].pins[j].x,rst->nets[i].pins[j].y);
        }*/
    }
    return 1;
}


int release(routingInst *rst){
	try {
        for (int i = 0; i < rst->numNets; i++) {
            net deleteNet = rst->nets[i];
            delete deleteNet.nroute.segments;
            delete deleteNet.pins;
        }
        delete rst->edgeCaps;
        delete rst;

    } catch(int e) {
        return 0;
    }
    return 1;
}

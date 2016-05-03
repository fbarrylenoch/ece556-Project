#include "ece556.h"

void handler(int sig) {
    void *array[10];
    size_t size;

    size = backtrace(array, 10);

    fprintf(stderr, "ERROR: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

bool compare (const net &n1, const net &n2){
    return n1.cost > n2.cost;
}

int readBenchmark(const char *fileName, routingInst *rst){
    vector<int> indices;
    int status;
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
            int x = atoi(token[1]);
            int y = atoi(token[2]);

            rst->gx = x;
            rst->gy = y;
            rst->numEdges = y*(x-1) + x*(y-1);
            rst->edgeCaps = (int *)malloc(rst->numEdges * sizeof(int));
            rst->edgeUtils = (int *)malloc(rst->numEdges * sizeof(int));
            rst->edgeWeight = (int *) malloc(rst->numEdges * sizeof(int));
            rst->edgeOver = (int *) malloc(rst->numEdges * sizeof(int));
            rst->edgeHis = (int *) malloc(rst->numEdges * sizeof(int));
            // going to set edgeUtils to 0 as a default and 1 if the edge is utilized
            for (int i = 0; i < rst->numEdges; i++) {
                rst->edgeUtils[i] = 0;
            }
        }
        // check for edge caps
        else if(strncmp(token[0], "capacity",64) == 0){
            rst->cap = atoi(token[1]);
            for (int i = 0; i < rst->numEdges; i++) 
                rst->edgeCaps[i] = rst->cap;
        }
        // check for number of nets
        else if(strncmp(token[0], "num", 64) == 0){
            rst->numNets = atoi(token[2]);
        }
        // check for nets
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
                rst->nets.push_back(*tempNet);
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
                point* p1 = new point;
                point* p2 = new point;
                p1->x = atoi(token[0]);
                p1->y = atoi(token[1]);
                p2->x = atoi(token[2]);
                p2->y = atoi(token[3]);
                int updatedCap = atoi(token[4]);

                status = getIndex(*p1, *p2, rst, &indices);
                if(status==0)
                    return 0;
                rst->edgeCaps[indices.front()] = updatedCap;
                indices.clear();
                delete p1;
                delete p2;
            }
        }
    }
    return 1;
}

int getIndex(point p1, point p2, routingInst *rst, std::vector<int>* indicies){
    int distance;
    indicies->push_back(-1);
    // check if horizontal line
    if (p1.y == p2.y) {
        // calculate index in rst->edgeCaps
        if ((p1.x-p2.x) == -1) {
            indicies->front() = (rst->gx-1)*p1.y + p1.x;
        }
        else if((p1.x-p2.x) < -1) {
            distance = abs(p1.x-p2.x);
            indicies->front() = (rst->gx-1)*p1.y + p1.x;
            for (int i=1; i < distance; ++i) {
                indicies->push_back((rst->gx-1)*p1.y + p1.x + i);
            }
        }            
        else if((p1.x-p2.x) == 1) {
            indicies->front() = (rst->gx-1)*p1.y + p2.x;
        }
        else {
            distance = abs(p1.x-p2.x);
            indicies->front() = (rst->gx-1)*p1.y + p2.x;
            for (int i=1; i < distance; ++i) {
                indicies->push_back((rst->gx-1)*p1.y + p2.x + i);
            }
        } // case where p1.x-p2.x > 1
    }
    // else vertical line
    else if(p1.x == p2.x){
        // calculate index in rst->edgeCaps
        if ((p1.y-p2.y) == -1) {
            indicies->front() = (rst->gy)*(rst->gx-1) + rst->gx*p1.y + p1.x;
        }
        else if ((p1.y-p2.y) < -1) {
            distance = abs(p1.y-p2.y);
            indicies->front() = (rst->gy)*(rst->gx-1) + rst->gx*p1.y + p1.x;
            for (int i=0; i < distance; i++)
                indicies->push_back((rst->gy)*(rst->gx-1) + rst->gx*(p1.y + i) + p1.x);
        }
        else if ((p1.y-p2.y) == 1) {
            indicies->front() = (rst->gy)*(rst->gx-1) + rst->gx*p2.y + p1.x;
        }
        else {
            distance = abs(p1.y-p2.y);
            indicies->front() = (rst->gy)*(rst->gx-1) + rst->gx*p2.y + p1.x;
            for (int i=0; i < distance; i++) {
                indicies->push_back((rst->gy)*(rst->gx-1) + rst->gx*(p2.y + i) + p1.x);
            }
        } // case where p1.y-p2.y > 1
    }
    if(indicies->front() == -1)
        return 0;
    return 1;
}

int printRoutingInstince (routingInst *rst){
    for(int i = 0; i < rst->numNets; ++i){
        printf("net n%d, has %d pins\n", i, rst->nets.at(i).numPins);
       /* printf("the pin connections are:\n");
        for(int j = 0; j < rst->nets.at(i).numPins; ++j){
            printf("\t(%d,%d)\n", rst->nets.at(i).pins[j].x,rst->nets.at(i).pins[j].y);
        }*/
    }
    return 1;
}

int solveRouting(routingInst *rst){
    vector<int> indices;
    int status;
    //get every net in rst
    for (int i = 0; i < rst->numNets; ++i) {
        net *tempNet = &(rst->nets[i]);
        tempNet->nroute.numSegs = 0;
        for (int j = 0; j < tempNet->numPins-1; j++) {

            point p1 = tempNet->pins[j];
            point p2 = tempNet->pins[j+1];
            if (p1.x == p2.x || p1.y == p2.y) {
                segment *tempSeg = new segment;
                tempNet->nroute.numSegs++;
                tempSeg->p1 = p1;
                tempSeg->p2 = p2;
                if (p1.x == p2.x) {
                    tempSeg->numEdges = abs(p1.y - p2.y);
                }
                else {
                    tempSeg->numEdges = abs(p1.x - p2.x);
                }

                status = getIndex(tempSeg->p1, tempSeg->p2, rst, &indices);
                if(status==0)
                    return 0;
                for(int k = 0; k < (int)indices.size(); ++k){
                    tempSeg->edges.push_back(indices[k]);
                    rst->edgeUtils[indices[k]]++;
                }
                indices.clear();
                tempNet->nroute.segments.push_back(*tempSeg);
                delete tempSeg;
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
                segment *tempSeg = new segment; 
                tempNet->nroute.numSegs++;
                tempSeg->p1 = topP;
                point midPoint;
                midPoint.x = topP.x;
                midPoint.y = bottomP.y;
                tempSeg->p2 = midPoint;
                tempSeg->numEdges = abs(topP.y - midPoint.y);

                // save edges to the segment
                status = getIndex(tempSeg->p1, tempSeg->p2, rst, &indices);
                if(status==0)
                    return 0;
                for(int k = 0; k < tempSeg->numEdges; ++k){
                    tempSeg->edges.push_back(indices[k]);
                    rst->edgeUtils[indices[k]]++;
                }
                indices.clear();
                tempNet->nroute.segments.push_back(*tempSeg);
                delete tempSeg;
                // add midPoint to bottomP
                tempSeg = new segment;
                tempNet->nroute.numSegs++;
                tempSeg->p1 = midPoint;
                tempSeg->p2 = bottomP;
                tempSeg->numEdges = abs(midPoint.x - bottomP.x);
                // save edges to the segment
                
                status = getIndex(tempSeg->p1, tempSeg->p2, rst, &indices);
                if(status==0)
                    return 0;
                for(int k = 0; k < tempSeg->numEdges; ++k){
                    tempSeg->edges.push_back(indices[k]);
                    rst->edgeUtils[tempSeg->edges[k]]++;
                }
                indices.clear();
                tempNet->nroute.segments.push_back(*tempSeg);
                delete tempSeg;
            }
        }
    }
    // calculate all edge weights at the beginning
    for(int i = 0; i < rst->numEdges; ++i){
        rst->edgeHis[i] = 0;
        rst->edgeOver[i] = max((rst->edgeUtils[i] - rst->edgeCaps[i]), 0);
        if(rst->edgeOver[i] > 0)
            rst->edgeHis[i] += 1;
        rst->edgeWeight[i] = rst->edgeOver[i] * rst->edgeHis[i];
    }
    return 1;
}

int RRR(routingInst *rst){
    time_t init_time, curr_time;
    double seconds;
    bool done = false;
    int status;

    time(&init_time);
    while(!done){
        calcNetCost(rst);
        printf("sending net n%d", rst->nets[0].id);
        status = RRRnet(rst, &rst->nets[0]);
		if(status==0)
			return 0;
        time(&curr_time);
        seconds = difftime(curr_time, init_time);
        printf("\tcurrent time %.2f\n", seconds);
        if(seconds >= 120)
            done = true;
    }
    return 1;
}

int calcEdgeWeights(route *newRoute, routingInst *rst){
    vector<int> indices;
    int status;
    for(int i = 0; i < newRoute->numSegs; ++i){
        status = getIndex(newRoute->segments[i].p1, newRoute->segments[i].p2, rst, &indices);
		if(status==0)
		    return 0;
        for(int j = 0; j < newRoute->segments[i].numEdges; ++j){
            if(rst->edgeOver[indices[j]] > 0)
                rst->edgeHis[indices[j]] += 1;
            rst->edgeOver[indices[j]] = max((rst->edgeUtils[indices[j]] 
                        - rst->edgeCaps[indices[j]]), 0);
            rst->edgeWeight[indices[j]] = rst->edgeOver[indices[j]] * 
                                           rst->edgeHis[indices[j]];
        }
        indices.clear();
    }
    return 1;
}

int calcRouteCost(route *newRoute, routingInst *rst){
    vector<int> indices;
    int cost, status;
    cost = 0;
    // over all segments in the route
    for(int i = 0; i < newRoute->numSegs; ++i){
        status = getIndex(newRoute->segments[i].p1, newRoute->segments[i].p2, rst, &indices);
        if(status==0)
            return 0;
        // over all edges in the segment
        for(int j = 0; j < newRoute->segments[i].numEdges; ++j){
            cost += rst->edgeWeight[indices[j]];
        }
        indices.clear();
    }
    return cost;
}

int calcNetCost(routingInst *rst){
    vector<int> indices;
    int status;
    // calc cost
    // over all nets
    for(int i = 0; i < rst->numNets; ++i){
        rst->nets.at(i).cost = 0;
        // over all segments in the route
        for(int j = 0; j < rst->nets.at(i).nroute.numSegs; ++j){
            status = getIndex(rst->nets.at(i).nroute.segments[j].p1, 
                    rst->nets.at(i).nroute.segments[j].p2, rst, &indices);
			if(status==0)
			    return 0;
            // over all edges in the segment
            for(int k = 0; k < rst->nets.at(i).nroute.segments[j].numEdges; ++k){
                rst->nets.at(i).cost += rst->edgeWeight[indices[k]];
            }
            indices.clear();
        }
    }
    //sort nets
    sort(rst->nets.begin(), rst->nets.end(), compare);
    return 1;
}

//will be called to rip up and reroute a net
int RRRnet(routingInst *rst, net *netRRR){
    int status;
    vector<int> indices;
    route routeRRR;
    //printf("\twe're in RRR for net n%d\n", netRRR->id);
    
    // routes for each shape attempt
    route routeBest;
    route routeL; // used to store the best L route option  
    route routeZ; // used to store the best Z route option  
    route routeRZ; // used to store the best Rotated Z route option  
    route routeU; // used to store the best U route option  
    route routeRU; // used to store the best U route option  
    int routeBest_weight,
        routeL_weight,
        routeZ_weight,
        routeRZ_weight,
        routeU_weight,
        routeRU_weight;
        

    routeRRR = netRRR->nroute;
    // rip up
    // over all segments in the route
    for(int i = 0; i < routeRRR.numSegs; ++i){
        status = getIndex(routeRRR.segments[i].p1, routeRRR.segments[i].p2, 
                            rst, &indices);
        if(status==0)
            return 0;
        // over all edges in the segment
        for(int j = 0; j < routeRRR.segments[i].numEdges; ++j)
            rst->edgeUtils[indices.at(j)]--;
        indices.clear();
    }
    //printf("\twe have ripped up the current route\n");

    // call edgeWeight
    calcEdgeWeights(&routeRRR,rst);
    // clear the segments of the route
    for(int i = 0; i < routeRRR.numSegs; ++i)
        routeRRR.segments[i].edges.clear();
    routeRRR.segments.clear();
    routeRRR.numSegs = 0;
    

    // reroute
    // loop for all pins
    for (int i = 0; i < (netRRR->numPins-1); i++) {
        //printf("\tchecking best route for pin%d to %d\n", i, (i+1));
        // try L shape
        routeL = shapeL(netRRR->pins[i],netRRR->pins[i+1],rst);
        // try Z shape 
        routeZ = shapeZ(netRRR->pins[i],netRRR->pins[i+1],rst);
        // try Rotated Z shape
        routeRZ = shapeRZ(netRRR->pins[i],netRRR->pins[i+1],rst);
        // try U shape
        routeU = shapeU(netRRR->pins[i],netRRR->pins[i+1],rst);
        // try Rotated U shape
        routeRU = shapeRU(netRRR->pins[i],netRRR->pins[i+1],rst);
        //printf("\twe have all of the candidate routes\n");
        
        
        // get weights of each route
        routeL_weight = calcRouteCost(&routeL, rst);
        routeZ_weight = calcRouteCost(&routeZ, rst);
        routeRZ_weight = calcRouteCost(&routeRZ, rst);
        routeU_weight = calcRouteCost(&routeU, rst);
        routeRU_weight = calcRouteCost(&routeRU, rst);
        //printf("\twe have calculated all weights for each route\n");
        
        
        //compare weights and take lowest cost route
        routeBest = routeL;
        routeBest_weight = routeL_weight;
        //printf("\tfor pins %d - %d, L is the best\n",i,i+1);
        if (routeZ_weight < routeBest_weight) {
            //printf("\tfor pins %d - %d, Z is the best\n",i,i+1);
            routeBest = routeZ;
            routeBest_weight = routeZ_weight;
        }
        if (routeRZ_weight < routeBest_weight) {
            //printf("\tfor pins %d - %d, RZ is the best\n",i,i+1);
            routeBest = routeRZ;
            routeBest_weight = routeRZ_weight;
        }
        if (routeU_weight < routeBest_weight) {
            //printf("\tfor pins %d - %d, U is the best\n",i,i+1);
            routeBest = routeU;
            routeBest_weight = routeU_weight;
        }
        if (routeRU_weight < routeBest_weight) {
            //printf("\tfor pins %d - %d, RU is the best\n",i,i+1);
            routeBest = routeRU;
            routeBest_weight = routeRU_weight;
        }
        // add to routeRRR_segments
        //printf("\twe have chosen the best route of them all\n");

        for (int j = 0; j < routeBest.numSegs; ++j)
            routeRRR.segments.push_back(routeBest.segments[j]);
        routeRRR.numSegs += routeBest.numSegs;
        //printf("\twe have added the best route to the route\n");
        
        // delete each route
        for(int j = 0; j < routeL.numSegs; ++j)
            routeL.segments[j].edges.clear();
        routeL.segments.clear();
        for(int j = 0; j < routeZ.numSegs; ++j)
            routeZ.segments[j].edges.clear();
        routeZ.segments.clear();
        for(int j = 0; j < routeRZ.numSegs; ++j)
            routeRZ.segments[j].edges.clear();
        routeRZ.segments.clear();
        for(int j = 0; j < routeU.numSegs; ++j)
            routeU.segments[j].edges.clear();
        routeU.segments.clear();
        for(int j = 0; j < routeRU.numSegs; ++j)
            routeRU.segments[j].edges.clear();
        routeRU.segments.clear();
    }
    // print routeRRR
    netRRR->nroute = routeRRR;
    // recalculate edgeWeight after Net has been rerouted
    calcEdgeWeights(&routeRRR,rst);
    return 1;
}

// returns best L shape route - return straight line if in line
route shapeL(point p1, point p2, routingInst *rst){
    vector<int> edgeL;
    vector<int> edgeComp;
    route routeComp; // used to compare different variations of each route shape  
    route routeL; // used to store the best L route option  
    point topP; // used to calculate route  
    point bottomP; // used to calculate route  
    point midPoint; // used to calculate route  
    int routeComp_weight, routeL_weight, status;
    
    // check if straight line, assign routeL if so  
    if (p1.x == p2.x || p1.y == p2.y) {
        routeL.numSegs = 0;

        segment *segL = new segment;
        routeL.numSegs++;
        segL->p1 = p1;
        segL->p2 = p2;
        if (p1.x == p2.x) {
            segL->numEdges = abs(p1.y - p2.y);
        }
        else {
            segL->numEdges = abs(p1.x - p2.x);
        }
        status = getIndex(segL->p1, segL->p2, rst, &edgeL);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeL.size(); ++i)
            segL->edges.push_back(edgeL[i]);
        edgeL.clear();
        routeL.segments.push_back(*segL);
        delete segL;
    }
    else{
        //printf("\t\tThe two points are an L shape\n");
        // assign L shape route  
        if (p1.y < p2.y) {
            topP = p2;
            bottomP = p1;
        }
        else {
            topP = p1;
            bottomP = p2;
        }
        
        routeL.numSegs = 0;
        // create a midPoint
        midPoint.x = topP.x;
        midPoint.y = bottomP.y;
        
        // assign top and bottom points
        segment *segL = new segment;
        routeL.numSegs++;
        segL->p1 = topP;        
        segL->p2 = midPoint;
        segL->numEdges = abs(topP.y - midPoint.y);
        // save edges to the segment
        status = getIndex(segL->p1, segL->p2, rst, &edgeL);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeL.size(); ++i)
            segL->edges.push_back(edgeL[i]);
        edgeL.clear();
        routeL.segments.push_back(*segL);
        delete segL;
        
        // add midPoint to bottomP
        segL = new segment;
        routeL.numSegs++;
        segL->p1 = midPoint;
        segL->p2 = bottomP;
        segL->numEdges = abs(midPoint.x - bottomP.x);
        // save edges to the segment
        status = getIndex(segL->p1, segL->p2, rst, &edgeL);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeL.size(); ++i)
            segL->edges.push_back(edgeL[i]);
        edgeL.clear();
        routeL.segments.push_back(*segL);
        delete segL;
        
        // assign L in other direction  
        if (p1.y < p2.y) {
            topP = p2;
            bottomP = p1;
        }
        else {
            topP = p1;
            bottomP = p2;
        }
        
        routeComp.numSegs = 0;
        // create a midPoint
        midPoint.y = topP.y;
        midPoint.x = bottomP.x;
        
        // assign top and bottom points
        segment *segComp = new segment;
        routeComp.numSegs++;
        segComp->p1 = topP;
        segComp->p2 = midPoint;
        segComp->numEdges = abs(topP.x - midPoint.x);
        // save edges to the segment
        status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeComp.size(); ++i)
            segComp->edges.push_back(edgeComp[i]);
        edgeComp.clear();
        routeComp.segments.push_back(*segComp);
        delete segComp;
       
       // add midPoint to bottomP
        segComp = new segment;
        routeComp.numSegs++;
        segComp->p1 = midPoint;
        segComp->p2 = bottomP;
        segComp->numEdges = abs(midPoint.y - bottomP.y);
        // save edges to the segment
        status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeComp.size(); ++i)
            segComp->edges.push_back(edgeComp[i]);
        edgeComp.clear();
        routeComp.segments.push_back(*segComp);
        delete segComp;

        routeComp_weight = calcRouteCost(&routeComp, rst); 
        routeL_weight = calcRouteCost(&routeL, rst);
        if (routeComp_weight < routeL_weight)
            routeL = routeComp;
        routeComp.segments.clear();
        
        // compare and take route with lowest weight  
        // need edge weight to be completed before this can be completed
    }
    // return the route
    return routeL;
}

// returns best Z shape route - call shapeL if too close to make a Z
route shapeZ(point p1, point p2, routingInst *rst){
    vector<int> edgeZ;
    vector<int> edgeComp;
    route routeComp; // used to compare different variations of each route shape  
    route routeZ; // used to store the best Z route option  
    point topPoint; // used to calculate route  
    point bottomPoint; // used to calculate route  
    point topMidPoint; // used to calculate route  
    point bottomMidPoint; // used to calculate route  
    int multiplier; // used to get proper dirction for route generation  
    int routeComp_weight, routeZ_weight, status;

    // call shapeL if too close
    if (abs(p1.x - p2.x) < 2) {
        routeZ = shapeL(p1,p2,rst);
    }

    // else assign Z shape route
    else {
        routeZ.numSegs = 0;
        // assign outside most Z route  
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
        segment *segZ = new segment;
        routeZ.numSegs++;
        segZ->p1 = topPoint;
        segZ->p2 = topMidPoint;
        segZ->numEdges = abs(topPoint.x - topMidPoint.x);
        // get top edges
        status = getIndex(segZ->p1, segZ->p2, rst, &edgeZ);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeZ.size(); ++i)
            segZ->edges.push_back(edgeZ[i]);
        edgeZ.clear();
        routeZ.segments.push_back(*segZ);
        delete segZ;

        // create mid segment
        segZ = new segment;
        routeZ.numSegs++;
        segZ->p1 = topMidPoint;
        segZ->p2 = bottomMidPoint;
        segZ->numEdges = abs(topMidPoint.y - bottomMidPoint.y); 
        // get mid edges
        status = getIndex(segZ->p1, segZ->p2, rst, &edgeZ);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeZ.size(); ++i)
            segZ->edges.push_back(edgeZ[i]);
        edgeZ.clear();
        routeZ.segments.push_back(*segZ);
        delete segZ;
        
        // create bottom segment
        segZ = new segment;
        routeZ.numSegs++;
        segZ->p1 = bottomMidPoint;
        segZ->p2 = bottomPoint;
        segZ->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
        // get bottom edges
        status = getIndex(segZ->p1, segZ->p2, rst, &edgeZ);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeZ.size(); ++i)
            segZ->edges.push_back(edgeZ[i]);
        edgeZ.clear();
        routeZ.segments.push_back(*segZ);
        delete segZ;

        //check every other possible Z route and compare weights
        for (int i = 1; i < (abs(p1.x - p2.x)-1); i++){
            routeComp.numSegs = 0;
            // create midpoints
            topMidPoint.x = (topPoint.x + ((1 + i)*multiplier));
            topMidPoint.y = topPoint.y;
            bottomMidPoint.x = (topPoint.x + ((1 + i)*multiplier));
            bottomMidPoint.y = bottomPoint.y;

            // create top segment
            segment *segComp = new segment;
            routeComp.numSegs++;
            segComp->p1 = topPoint;
            segComp->p2 = topMidPoint;
            segComp->numEdges = abs(topPoint.x - topMidPoint.x);
            // get top edges
            status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
            if(status==0)
                fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
            for(int j = 0; j < (int)edgeComp.size(); ++j){
                segComp->edges.push_back(edgeComp[j]);
            }
            edgeComp.clear();
            routeComp.segments.push_back(*segComp);
            delete segComp;

            // create mid segment
            segComp = new segment;
            routeComp.numSegs++;
            segComp->p1 = topMidPoint;
            segComp->p2 = bottomMidPoint;
            segComp->numEdges = abs(topMidPoint.y - bottomMidPoint.y);
            // get mid edges
            status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
            if(status==0)
                fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
            for(int j = 0; j < (int)edgeComp.size(); ++j)
                segComp->edges.push_back(edgeComp[j]);
            edgeComp.clear();
            routeComp.segments.push_back(*segComp);
            delete segComp;

            // create bottom segment
            segComp = new segment;
            routeComp.numSegs++;
            segComp->p1 = bottomMidPoint;
            segComp->p2 = bottomPoint;
            segComp->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
            // get bottom edges
            status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
            if(status==0)
                fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
            for(int j = 0; j < (int)edgeComp.size(); ++j)
                segComp->edges.push_back(edgeComp[j]);
            edgeComp.clear();
            routeComp.segments.push_back(*segComp);
            delete segComp;

            // check if weight is lower than routeZ, assign to routeZ if true
            routeComp_weight = calcRouteCost(&routeComp, rst); 
            routeZ_weight = calcRouteCost(&routeZ, rst);
            if (routeComp_weight < routeZ_weight)
                routeZ = routeComp;
            
            for(int j = 0; j < routeComp.numSegs; ++j)
                routeComp.segments[j].edges.clear();
            routeComp.segments.clear();
            
        }
    }
    // return the route
    return routeZ;
}

// returns best Rotated Z shape route - call shapeL if too close to make a RZ
route shapeRZ(point p1, point p2, routingInst *rst){
    vector<int> edgeRZ;
    vector<int> edgeComp;
    route routeComp; // used to compare different variations of each route shape  
    route routeRZ; // used to store the best Z route option  
    point leftMidPoint; // used to calculate route  
    point rightMidPoint; // used to calculate route  
    point leftPoint; // used to calculate route  
    point rightPoint; // used to calculate route  
    int multiplier; // used to get proper direction for route generation  
    int routeComp_weight, routeRZ_weight, status;

    // call shapeL if too close
    if (abs(p1.y - p2.y) < 2){
        routeRZ = shapeL(p1,p2,rst);
    }

    // else loop for all possible RZ routes
    else {
        routeRZ.numSegs = 0;
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
        segment *segRZ = new segment;
        routeRZ.numSegs++;
        segRZ->p1 = leftPoint;
        segRZ->p2 = leftMidPoint;
        segRZ->numEdges = abs(leftPoint.y - leftMidPoint.y);
        // get left edges
        status = getIndex(segRZ->p1, segRZ->p2, rst, &edgeRZ);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeRZ.size(); ++i)
            segRZ->edges.push_back(edgeRZ[i]);
        edgeRZ.clear();
        routeRZ.segments.push_back(*segRZ);
        delete segRZ;

        // create mid segment
        segRZ = new segment;
        routeRZ.numSegs++;
        segRZ->p1 = leftMidPoint;
        segRZ->p2 = rightMidPoint;
        segRZ->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
        // get mid edges
        status = getIndex(segRZ->p1, segRZ->p2, rst, &edgeRZ);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeRZ.size(); ++i)
            segRZ->edges.push_back(edgeRZ[i]);
        edgeRZ.clear();
        routeRZ.segments.push_back(*segRZ);
        delete segRZ;

        // create right segment
        segRZ = new segment;
        routeRZ.numSegs++;
        segRZ->p1 = rightMidPoint;
        segRZ->p2 = rightPoint;
        segRZ->numEdges = abs(rightMidPoint.y - rightPoint.y);
        // get right edges
        status = getIndex(segRZ->p1, segRZ->p2, rst, &edgeRZ);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeRZ.size(); ++i)
            segRZ->edges.push_back(edgeRZ[i]);
        edgeRZ.clear();
        routeRZ.segments.push_back(*segRZ);
        delete segRZ;

        //check every other possible RZ route and compare weights
        for (int i = 1; i < (abs(p1.y - p2.y)-1); i++){
            routeComp.numSegs = 0;
            // create midpoints
            leftMidPoint.x = leftPoint.x;
            leftMidPoint.y = leftPoint.y - ((1 + i) * multiplier);
            rightMidPoint.x = rightPoint.x;
            rightMidPoint.y = leftPoint.y - ((1 + i) * multiplier);

            // create left segment
            segment *segComp = new segment;
            routeComp.numSegs++;
            segComp->p1 = leftPoint;
            segComp->p2 = leftMidPoint;
            segComp->numEdges = abs(leftPoint.y - leftMidPoint.y);
            // get left edges
            status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
            if(status==0)
                fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
            for(int j = 0; j < (int)edgeComp.size(); ++j)
                segComp->edges.push_back(edgeComp[j]);
            edgeComp.clear();
            routeComp.segments.push_back(*segComp);
            delete segComp;

            // create mid segment
            segComp = new segment;
            routeComp.numSegs++;
            segComp->p1 = leftMidPoint;
            segComp->p2 = rightMidPoint;
            segComp->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
            // get mid edges
            status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
            if(status==0)
                fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
            for(int j = 0; j < (int)edgeComp.size(); ++j)
                segComp->edges.push_back(edgeComp[j]);
            edgeComp.clear();
            routeComp.segments.push_back(*segComp);
            delete segComp;

            // create right segment
            segComp = new segment;
            routeComp.numSegs++;
            segComp->p1 = rightMidPoint;
            segComp->p2 = rightPoint;
            segComp->numEdges = abs(rightMidPoint.y - rightPoint.y);
            // get right edges
            status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
            if(status==0)
                fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
            for(int j = 0; j < (int)edgeComp.size(); ++j){;
                segComp->edges.push_back(edgeComp[j]);
            }
            edgeComp.clear();
            routeComp.segments.push_back(*segComp);
            delete segComp;

            // check if weight is lower than routeRZ, assign to routeRZ if true
            routeComp_weight = calcRouteCost(&routeComp, rst); 
            routeRZ_weight = calcRouteCost(&routeRZ, rst);
            if (routeComp_weight < routeRZ_weight)
                routeRZ = routeComp;

            for(int j = 0; j < routeComp.numSegs; ++j)
                routeComp.segments[j].edges.clear();
            routeComp.segments.clear();
        }
    }
    // return the route
    return routeRZ;
}


// returns best U shape route
route shapeU(point p1, point p2, routingInst *rst){
    vector<int> edgeU;
    vector<int> edgeComp;
    route routeComp; // used to compare different variations of each route shape  
    route routeU; // used to store the best U route option  
    point bottomPoint; // used to store the point with a lower y value  
    point leftPoint;
    point topPoint;
    point rightPoint;
    point leftMidPoint;
    point rightMidPoint;
    int routeComp_weight, routeU_weight, status;
    
     // return empty route if same x values
    if (p1.x == p2.x){
        return routeU;
    }

    // find bottom point and top point
    if (p1.y < p2.y) {
        bottomPoint = p1;
        topPoint = p2;
    }
    else {
        bottomPoint = p2;
        topPoint = p1;
    }

    // find left point and right point
    if (p1.x < p2.x) {
        leftPoint = p1;
        rightPoint = p2;
    }
    else {
        leftPoint = p2;
        rightPoint = p1;
    }

    // return empty route if bottomPoint.y == 0
    if (bottomPoint.y == 0) {
        return routeU;
    }
    else {
        // allocate size
        routeU.numSegs = 0;
        //find the best ShapeU route
        leftMidPoint.x = leftPoint.x;
        leftMidPoint.y = (bottomPoint.y - 1);
        rightMidPoint.x = rightPoint.x;
        rightMidPoint.y = (bottomPoint.y - 1);

        // create left segment
        segment *segU = new segment;
        routeU.numSegs++;
        segU->p1 = leftPoint;
        segU->p2 = leftMidPoint;
        segU->numEdges = abs(leftPoint.y - leftMidPoint.y);
        // get left edges
        status = getIndex(segU->p1, segU->p2, rst, &edgeU);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeU.size(); ++i)
            segU->edges.push_back(edgeU[i]);
        edgeU.clear();
        routeU.segments.push_back(*segU);
        delete segU;
        
        // create mid segment
        segU = new segment;
        routeU.numSegs++;
        segU->p1 = leftMidPoint;
        segU->p2 = rightMidPoint;
        segU->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
        // get mid edges
        status = getIndex(segU->p1, segU->p2, rst, &edgeU);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeU.size(); ++i)
            segU->edges.push_back(edgeU[i]);
        edgeU.clear();
        routeU.segments.push_back(*segU);
        delete segU;
        
        // create right segment
        segU = new segment;
        routeU.numSegs++;
        segU->p1 = rightMidPoint;
        segU->p2 = rightPoint;
        segU->numEdges = abs(rightMidPoint.y - rightPoint.y);
        // get right edges
        status = getIndex(segU->p1, segU->p2, rst, &edgeU);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeU.size(); ++i)
            segU->edges.push_back(edgeU[i]);
        edgeU.clear();
        routeU.segments.push_back(*segU);
        delete segU;
        
        //check every other possible U route and compare weights
        if (bottomPoint.y - 2 >= 0) { // may have to be -1  
            for (int i = 2; i <= bottomPoint.y; i++) {
                routeComp.numSegs = 0;
                // create midpoints
                leftMidPoint.x = leftPoint.x;
                leftMidPoint.y = (bottomPoint.y - i);
                rightMidPoint.x = rightPoint.x;
                rightMidPoint.y = (bottomPoint.y - i);
                
                // create left segment
                segment *segComp = new segment;
                routeComp.numSegs++;
                segComp->p1 = leftPoint;
                segComp->p2 = leftMidPoint;
                segComp->numEdges = abs(leftPoint.y - leftMidPoint.y);
                // get left edges
                status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
                if(status==0)
                    fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
                for(int j = 0; j < (int)edgeComp.size(); ++j)
                    segComp->edges.push_back(edgeComp[j]);
                edgeComp.clear();
                routeComp.segments.push_back(*segComp);
                delete segComp;
                
                // create mid segment
                segComp = new segment;
                routeComp.numSegs++;
                segComp->p1 = leftMidPoint;
                segComp->p2 = rightMidPoint;
                segComp->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
                // get mid edges
                status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
                if(status==0)
                    fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");                
                for(int j = 0; j < (int)edgeComp.size(); ++j)
                    segComp->edges.push_back(edgeComp[j]);
                edgeComp.clear();
                routeComp.segments.push_back(*segComp);
                delete segComp;

                
                // create right segment
                segComp = new segment;
                routeComp.numSegs++;
                segComp->p1 = rightMidPoint;
                segComp->p2 = rightPoint;
                segComp->numEdges = abs(rightMidPoint.y - rightPoint.y);
                // get right edges
                status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
                if(status==0)
                    fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");                
                for(int j = 0; j < (int)edgeComp.size(); ++j)
                    segComp->edges.push_back(edgeComp[j]);
                edgeComp.clear();
                routeComp.segments.push_back(*segComp);
                delete segComp;

                // check if weight is lower than routeRZ, assign to routeRZ if true
                routeComp_weight = calcRouteCost(&routeComp, rst); 
                routeU_weight = calcRouteCost(&routeU, rst);
                if (routeComp_weight < routeU_weight)
                    routeU = routeComp;
            
                for(int j = 0; j < routeComp.numSegs; ++j)
                    routeComp.segments[j].edges.clear();
                routeComp.segments.clear();
            }
        }
        // return the route
        return routeU;
    }
}

// returns the best Rotated U shape route
route shapeRU(point p1, point p2, routingInst *rst){
    printf("in RU with p1: (%d,%d) and p2: (%d,%d)\n",p1.x,p1.y,p2.x,p2.y);
    vector<int> edgeRU;
    vector<int> edgeComp;
    route routeComp; // used to compare different variations of each route shape  
    route routeRU; // used to store the best U route option  
    point bottomPoint; // used to store the point with a lower y value  
    point leftPoint;
    point topPoint;
    point rightPoint;
    point leftMidPoint;
    point rightMidPoint;
    int routeComp_weight, routeRU_weight, status;
    
     // return empty route if same x values
    if (p1.x == p2.x){
        printf("returning empty route\n");
        return routeRU;
    }

    // find bottom point and top point
    if (p1.y < p2.y) {
        bottomPoint = p1;
        topPoint = p2;
    }
    else {
        bottomPoint = p2;
        topPoint = p1;
    }

    // find left point and right point
    if (p1.x < p2.x) {
        leftPoint = p1;
        rightPoint = p2;
    }
    else {
        leftPoint = p2;
        rightPoint = p1;
    }

    // return empty route if topPoint.y == (rst->gy - 1)
    if (topPoint.y == (rst->gy - 1)) { // may have to be -1  
        return routeRU;
    }
    else {
        // allocate size
        routeRU.numSegs = 0;
        //find the best ShapeU route
        leftMidPoint.x = leftPoint.x;
        leftMidPoint.y = (topPoint.y + 1);
        rightMidPoint.x = rightPoint.x;
        rightMidPoint.y = (topPoint.y + 1);

        // create left segment
        segment *segRU = new segment;
        routeRU.numSegs++;
        segRU->p1 = leftPoint;
        segRU->p2 = leftMidPoint;
        segRU->numEdges = abs(leftPoint.y - leftMidPoint.y);
        // get left edges
        status = getIndex(segRU->p1, segRU->p2, rst, &edgeRU);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeRU.size(); ++i){
            //printf("\tedgeRU[i]= %d\n",edgeRU[i]);
            segRU->edges.push_back(edgeRU[i]);
        }
        edgeRU.clear();
        routeRU.segments.push_back(*segRU);
        delete segRU;
        
        // create mid segment
        segRU = new segment;
        routeRU.numSegs++;
        segRU->p1 = leftMidPoint;
        segRU->p2 = rightMidPoint;
        segRU->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
        // get mid edges
        status = getIndex(segRU->p1, segRU->p2, rst, &edgeRU);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeRU.size(); ++i){
            //printf("\tedgeRU[i]= %d\n",edgeRU[i]);
            segRU->edges.push_back(edgeRU[i]);
        }
        edgeRU.clear();
        routeRU.segments.push_back(*segRU);
        delete segRU;
        
        // create right segment
        segRU = new segment;
        routeRU.numSegs++;
        segRU->p1 = rightMidPoint;
        segRU->p2 = rightPoint;
        segRU->numEdges = abs(rightMidPoint.y - rightPoint.y);
        // get right edges
        status = getIndex(segRU->p1, segRU->p2, rst, &edgeRU);
        if(status==0)
            fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
        for(int i = 0; i < (int)edgeRU.size(); ++i){
            //printf("\tedgeRU[i]= %d\n",edgeRU[i]);
            segRU->edges.push_back(edgeRU[i]);
        }
        edgeRU.clear();
        routeRU.segments.push_back(*segRU);
        delete segRU;
        
        //check every other possible RU route and compare weights
        if (topPoint.y <= (rst->gy - 2)) {
            for (int i = 2; i <= (rst->gx - topPoint.y); i++) {
                routeComp.numSegs = 0;
                // create midpoints
                leftMidPoint.x = leftPoint.x;
                leftMidPoint.y = (topPoint.y + i);
                rightMidPoint.x = rightPoint.x;
                rightMidPoint.y = (topPoint.y + i);
               
                // create left segment
                segment *segComp = new segment;
                routeComp.numSegs++;
                segComp->p1 = leftPoint;
                segComp->p2 = leftMidPoint;
                segComp->numEdges = abs(leftPoint.y - leftMidPoint.y);
                // get left edges
                status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
                if(status==0)
                    fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
                for(int j = 0; j < (int)edgeComp.size(); ++j){
                    //printf("\tedgeComp[j]= %d\n",edgeComp[j]);
                    segComp->edges.push_back(edgeComp[j]);
                }
                edgeComp.clear();
                routeComp.segments.push_back(*segComp);
                delete segComp;
                
                // create mid segment
                segComp = new segment;
                routeComp.numSegs++;
                segComp->p1 = leftMidPoint;
                segComp->p2 = rightMidPoint;
                segComp->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
                // get mid edges
                status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
                if(status==0)
                    fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
                for(int j = 0; j < (int)edgeComp.size(); ++j){
                    //printf("\tedgeComp[j]= %d\n",edgeComp[j]);
                    segComp->edges.push_back(edgeComp[j]);
                }
                edgeComp.clear();
                routeComp.segments.push_back(*segComp);
                delete segComp;
                
                // create right segment
                segComp = new segment;
                routeComp.numSegs++;
                segComp->p1 = rightMidPoint;
                segComp->p2 = rightPoint;
                segComp->numEdges = abs(rightMidPoint.y - rightPoint.y);
                // get right edges
                status = getIndex(segComp->p1, segComp->p2, rst, &edgeComp);
                if(status==0)
                    fprintf(stderr, "\n\t==>\twe got an error from getIndex\n");
                for(int j = 0; j < (int)edgeComp.size(); ++j){
                    //printf("\tedgeComp[j]= %d\n",edgeComp[j]);
                    segComp->edges.push_back(edgeComp[j]);
                }
                edgeComp.clear();
                routeComp.segments.push_back(*segComp);
                delete segComp;

                // check if weight is lower than routeRZ, assign to routeRZ if true
                routeComp_weight = calcRouteCost(&routeComp, rst); 
                routeRU_weight = calcRouteCost(&routeRU, rst);
                if (routeComp_weight < routeRU_weight)
                    routeRU = routeComp;
                for(int j = 0; j < routeComp.numSegs; ++j)
                    routeComp.segments[j].edges.clear();
                routeComp.segments.clear();
            }
        }
        // return the route
        printf("leaving RU\n");
        return routeRU;
    }
}

int writeOutput(const char *outRouteFile, routingInst *rst){
	try {
        ofstream output;
        output.open(outRouteFile);

        for (int i = 0; i < rst->numNets; i++) {

            output << "n" << rst->nets.at(i).id << "\n";

            for (int j = 0; j < rst->nets.at(i).nroute.numSegs; j++) {
                segment tmpSeg = rst->nets.at(i).nroute.segments[j];
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

int printWireLen(routingInst *rst){
    int TWL = 0, TOF = 0;
    for(int i = 0; i < rst->numEdges; ++i)
        if((rst->edgeUtils[i] - rst->edgeCaps[i]) > 0)
            TOF += (rst->edgeUtils[i] - rst->edgeCaps[i]);
    for(int i = 0; i < (int)rst->nets.size(); ++i)
        for(int j = 0; j < rst->nets[i].nroute.numSegs; ++j)
            TWL += rst->nets[i].nroute.segments[j].numEdges;
    
    printf("total wirelength is %d\ntotal overflow is %d\n", TWL, TOF);
    return 1;
}

int release(routingInst *rst){
	try {
        for (int i = 0; i < rst->numNets; i++) {
            net deleteNet = rst->nets.at(i);
            deleteNet.nroute.segments.clear();
            free(deleteNet.pins);
        }
        free(rst->edgeCaps);
        free(rst->edgeUtils);
        free(rst->edgeWeight);
        free(rst->edgeOver);
        free(rst->edgeHis);
        
        delete rst;

    } catch(int e) {
        return 0;
    }
    return 1;
}

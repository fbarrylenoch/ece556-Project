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
            index[0] = (rst->gy)*(rst->gx-1) + rst->gx*p1.y + p1.x;
        }
        else if ((p1.y-p2.y) < -1) {
            distance = abs(p1.y-p2.y);
            //printf("\t\t6) distance = %d\n", distance);
            delete index;
            index = (int *)malloc(distance*sizeof(int));
            for (int i=0; i < distance; i++) {
                index[i] = (rst->gy)*(rst->gx-1) + rst->gx*(p1.y + i) + p1.x;
            }
        }
        else if ((p1.y-p2.y) == 1) {
            delete index;
            //printf("\t\t7) distance = %d\n", distance);
            index = (int *)malloc(sizeof(int));
            index[0] = (rst->gy)*(rst->gx-1) + rst->gx*p2.y + p1.x;
        }
        else {
            distance = abs(p1.y-p2.y);
            //printf("\t\t8) distance = %d\n", distance);
            delete index;
            index = (int *)malloc(distance*sizeof(int));
            for (int i=0; i < distance; i++) {
                index[i] = (rst->gy)*(rst->gx-1) + rst->gx*(p2.y + i) + p1.x;
            }
        } // case where p1.y-p2.y > 1
    }
    return index;
}

//will be called to rip up and reroute a net
void RRR(routingInst *rst, net *netRRR){
    int* indices;
    point p1;
    point p2;
    route routeRRR;
    segment *segRRR;
    
    // routes for each shape attempt
    route routeL; /* used to store the best L route option */
    route routeZ; /* used to store the best Z route option */
    route routeRZ; /* used to store the best Rotated Z route option */
    route routeU; /* used to store the best U route option */
    route routeRU; /* used to store the best U route option */
    route routeC; /* used to store the best U route option */
    route routeRC; /* used to store the best C route option */
    

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
            //rst->edgeCaps[segRRR->edges[i]]++;
            rst->edgeUtils[segRRR->edges[i]]--;
        }
    }

    // call edgeWeight

    // clear the segments of the route
    routeRRR.numSegs = 0;
    delete routeRRR.segments;

    // reroute
    /* allocate enough memory for every connection from pin to pin to be the largest number of segments */
    routeRRR.segments = (segment *)malloc((netRRR->numPins-1)*3*sizeof(segment));
    
    // loop for all pins

        // try L shape
        routeL = shapeL(p1,p2,rst);
        // try Z shape 
        routeZ = shapeZ(p1,p2,rst);
        // try Rotated Z shape
        routeRZ = shapeRZ(p1,p2,rst);
        // try U shape
        routeU = shapeU(p1,p2,rst);
        // try Rotated U shape
        routeRU = shapeRU(p1,p2,rst);
        // try C shape
        routeC = shapeC(p1,p2,rst);
        // try Rotated C shape
        routeRC = shapeRC(p1,p2,rst);
        // compare each attempt, take best option

    // recalculate edgeWeight after Net has been rerouted

    // get new net ordering

    // recall RRR function


    return;
}

// returns best L shape route - return straight line if in line
route shapeL(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeL; /* used to store the best L route option */
    point topP; /* used to calculate route */
    point bottomP; /* used to calculate route */
    point midPoint; /* used to calculate route */
    int routeComp_weight, routeL_weight;
    
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

        routeComp_weight = calcRouteCost(&routeComp, rst); 
        routeL_weight = calcRouteCost(&routeL, rst);
        if (routeComp_weight < routeL_weight)
            routeL = routeComp;
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
    int routeComp_weight, routeZ_weight; 

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
            routeComp_weight = calcRouteCost(&routeComp, rst); 
            routeZ_weight = calcRouteCost(&routeZ, rst);
            if (routeComp_weight < routeZ_weight)
                routeZ = routeComp;
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
    int routeComp_weight, routeRZ_weight;


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
            routeComp_weight = calcRouteCost(&routeComp, rst); 
            routeRZ_weight = calcRouteCost(&routeRZ, rst);
            if (routeComp_weight < routeRZ_weight)
                routeRZ = routeComp;
        }
    }
    // return the route
    return routeRZ;
}

// returns best U shape route
route shapeU(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeU; /* used to store the best U route option */
    point bottomPoint; /* used to store the point with a lower y value */
    point leftPoint;
    point topPoint;
    point rightPoint;
    point leftMidPoint;
    point rightMidPoint;
    int routeComp_weight, routeU_weight;
    
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
        routeU.segments = (segment *)malloc(3*sizeof(segment));
        routeU.numSegs = 0;
        routeComp.segments = (segment *)malloc(3*sizeof(segment));
        routeComp.numSegs = 0;
        //find the best ShapeU route
        leftMidPoint.x = leftPoint.x;
        leftMidPoint.y = (bottomPoint.y - 1);
        rightMidPoint.x = rightPoint.x;
        rightMidPoint.y = (bottomPoint.y - 1);

        // create left segment
        segment *segU = &(routeU.segments[routeU.numSegs]);
        routeU.numSegs++;
        segU->p1 = leftPoint;
        segU->p2 = leftMidPoint;
        segU->numEdges = abs(leftPoint.y - leftMidPoint.y);
        // get left edges
        int* edgesU = getIndex(leftPoint,leftMidPoint,rst);
        segU->edges = (int *)malloc(segU->numEdges*sizeof(int));
        segU->edges = edgesU;
        // create mid segment
        segU = &(routeU.segments[routeU.numSegs]);
        routeU.numSegs++;
        segU->p1 = leftMidPoint;
        segU->p2 = rightMidPoint;
        segU->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
        // get mid edges
        edgesU = getIndex(leftMidPoint,rightMidPoint,rst);
        segU->edges = (int *)malloc(segU->numEdges*sizeof(int));
        segU->edges = edgesU;
        // create right segment
        segU = &(routeU.segments[routeU.numSegs]);
        routeU.numSegs++;
        segU->p1 = rightMidPoint;
        segU->p2 = rightPoint;
        segU->numEdges = abs(rightMidPoint.y - rightPoint.y);
        // get right edges
        edgesU = getIndex(rightMidPoint,rightPoint,rst);
        segU->edges = (int *)malloc(segU->numEdges*sizeof(int));
        segU->edges = edgesU;
        //check every other possible U route and compare weights
        if (bottomPoint.y - 2 >= 0) { /* may have to be -1 */
            for (int i = 2; i <= bottomPoint.y; i++) {
                routeComp.numSegs = 0;
                // create midpoints
                leftMidPoint.x = leftPoint.x;
                leftMidPoint.y = (bottomPoint.y - i);
                rightMidPoint.x = rightPoint.x;
                rightMidPoint.y = (bottomPoint.y - i);
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
                routeComp_weight = calcRouteCost(&routeComp, rst); 
                routeU_weight = calcRouteCost(&routeU, rst);
                if (routeComp_weight < routeU_weight)
                    routeU = routeComp;
            }
        }
        // return the route
        return routeU;
    }
}

// returns the best Rotated U shape route
route shapeRU(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeRU; /* used to store the best U route option */
    point bottomPoint; /* used to store the point with a lower y value */
    point leftPoint;
    point topPoint;
    point rightPoint;
    point leftMidPoint;
    point rightMidPoint;
    int routeComp_weight, routeRU_weight;
    
     // return empty route if same x values
    if (p1.x == p2.x){
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
    if (topPoint.y == (rst->gy - 1)) { /* may have to be -1 */
        return routeRU;
    }
    else {
        // allocate size
        routeRU.segments = (segment *)malloc(3*sizeof(segment));
        routeRU.numSegs = 0;
        routeComp.segments = (segment *)malloc(3*sizeof(segment));
        routeComp.numSegs = 0;
        //find the best ShapeU route
        leftMidPoint.x = leftPoint.x;
        leftMidPoint.y = (topPoint.y + 1);
        rightMidPoint.x = rightPoint.x;
        rightMidPoint.y = (topPoint.y + 1);

        // create left segment
        segment *segRU = &(routeRU.segments[routeRU.numSegs]);
        routeRU.numSegs++;
        segRU->p1 = leftPoint;
        segRU->p2 = leftMidPoint;
        segRU->numEdges = abs(leftPoint.y - leftMidPoint.y);
        // get left edges
        int* edgesRU = getIndex(leftPoint,leftMidPoint,rst);
        segRU->edges = (int *)malloc(segRU->numEdges*sizeof(int));
        segRU->edges = edgesRU;
        // create mid segment
        segRU = &(routeRU.segments[routeRU.numSegs]);
        routeRU.numSegs++;
        segRU->p1 = leftMidPoint;
        segRU->p2 = rightMidPoint;
        segRU->numEdges = abs(leftMidPoint.x - rightMidPoint.x);
        // get mid edges
        edgesRU = getIndex(leftMidPoint,rightMidPoint,rst);
        segRU->edges = (int *)malloc(segRU->numEdges*sizeof(int));
        segRU->edges = edgesRU;
        // create right segment
        segRU = &(routeRU.segments[routeRU.numSegs]);
        routeRU.numSegs++;
        segRU->p1 = rightMidPoint;
        segRU->p2 = rightPoint;
        segRU->numEdges = abs(rightMidPoint.y - rightPoint.y);
        // get right edges
        edgesRU = getIndex(rightMidPoint,rightPoint,rst);
        segRU->edges = (int *)malloc(segRU->numEdges*sizeof(int));
        segRU->edges = edgesRU;
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
                routeComp_weight = calcRouteCost(&routeComp, rst); 
                routeRU_weight = calcRouteCost(&routeRU, rst);
                if (routeComp_weight < routeRU_weight)
                    routeRU = routeComp;
            }
        }
        // return the route
        return routeRU;
    }
}

// returns best C shape route
route shapeC(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeC; /* used to store the best U route option */
    point bottomPoint; /* used to store the point with a lower y value */
    point leftPoint;
    point topPoint;
    point rightPoint;
    point topMidPoint;
    point bottomMidPoint;
    int routeComp_weight, routeC_weight;
    
     // return empty route if same y values
    if (p1.y == p2.y){
        return routeC;
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

    // return empty route if bottomPoint.x == 0
    if (leftPoint.x == 0) {
        return routeC;
    }
    else {
        // allocate size
        routeC.segments = (segment *)malloc(3*sizeof(segment));
        routeC.numSegs = 0;
        routeComp.segments = (segment *)malloc(3*sizeof(segment));
        routeComp.numSegs = 0;
        //find the best ShapeC route
        topMidPoint.x = (leftPoint.x - 1);
        topMidPoint.y = topPoint.y;
        bottomMidPoint.x = (leftPoint.x - 1);
        bottomMidPoint.y = bottomPoint.y;

        // create top segment
        segment *segC = &(routeC.segments[routeC.numSegs]);
        routeC.numSegs++;
        segC->p1 = topPoint;
        segC->p2 = topMidPoint;
        segC->numEdges = abs(topMidPoint.x - topPoint.x);
        // get top edges
        int* edgesC = getIndex(topPoint,topMidPoint,rst);
        segC->edges = (int *)malloc(segC->numEdges*sizeof(int));
        segC->edges = edgesC;
        // create mid segment
        segC = &(routeC.segments[routeC.numSegs]);
        routeC.numSegs++;
        segC->p1 = topMidPoint;
        segC->p2 = bottomMidPoint;
        segC->numEdges = abs(topMidPoint.y - bottomMidPoint.y);
        // get mid edges
        edgesC = getIndex(topMidPoint,bottomMidPoint,rst);
        segC->edges = (int *)malloc(segC->numEdges*sizeof(int));
        segC->edges = edgesC;
        // create bottom segment
        segC = &(routeC.segments[routeC.numSegs]);
        routeC.numSegs++;
        segC->p1 = bottomMidPoint;
        segC->p2 = bottomPoint;
        segC->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
        // get right edges
        edgesC = getIndex(bottomMidPoint,bottomPoint,rst);
        segC->edges = (int *)malloc(segC->numEdges*sizeof(int));
        segC->edges = edgesC;
        //check every other possible C route and compare weights
        if ((leftPoint.x - 2) >= 0) { /* may have to be -1 */
            for (int i = 2; i <= bottomPoint.y; i++) {
                routeComp.numSegs = 0;
                // create midpoints
                topMidPoint.x = (leftPoint.x - i);
                topMidPoint.y = topPoint.y;
                bottomMidPoint.x = (leftPoint.x - i);
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
                // create right segment
                segComp = &(routeComp.segments[routeComp.numSegs]);
                routeComp.numSegs++;
                segComp->p1 = bottomMidPoint;
                segComp->p2 = bottomPoint;
                segComp->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
                // get right edges
                edgesComp = getIndex(bottomMidPoint,bottomPoint,rst);
                segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
                segComp->edges = edgesComp;

                // check if weight is lower than routeC, assign to routeC if true
                routeComp_weight = calcRouteCost(&routeComp, rst); 
                routeC_weight = calcRouteCost(&routeC, rst);
                if (routeComp_weight < routeC_weight)
                    routeC = routeComp;
            }
        }
        // return the route
        return routeC;
    }
}

// returns best Rotated C shape route
route shapeRC(point p1, point p2, routingInst *rst){
    route routeComp; /* used to compare different variations of each route shape */
    route routeRC; /* used to store the best U route option */
    point bottomPoint; /* used to store the point with a lower y value */
    point leftPoint;
    point topPoint;
    point rightPoint;
    point topMidPoint;
    point bottomMidPoint;
    int routeComp_weight, routeRC_weight;
    
     // return empty route if same y values
    if (p1.y == p2.y){
        return routeRC;
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

    // return empty route if rightPoint.x == (rst->gx - 1)
    if (rightPoint.x == (rst->gx - 1)) {
        return routeRC;
    }
    else {
        // allocate size
        routeRC.segments = (segment *)malloc(3*sizeof(segment));
        routeRC.numSegs = 0;
        routeComp.segments = (segment *)malloc(3*sizeof(segment));
        routeComp.numSegs = 0;
        //find the best ShapeC route
        topMidPoint.x = (rightPoint.x + 1);
        topMidPoint.y = topPoint.y;
        bottomMidPoint.x = (rightPoint.x + 1);
        bottomMidPoint.y = bottomPoint.y;

        // create top segment
        segment *segRC = &(routeRC.segments[routeRC.numSegs]);
        routeRC.numSegs++;
        segRC->p1 = topPoint;
        segRC->p2 = topMidPoint;
        segRC->numEdges = abs(topMidPoint.x - topPoint.x);
        // get top edges
        int* edgesRC = getIndex(topPoint,topMidPoint,rst);
        segRC->edges = (int *)malloc(segRC->numEdges*sizeof(int));
        segRC->edges = edgesRC;
        // create mid segment
        segRC = &(routeRC.segments[routeRC.numSegs]);
        routeRC.numSegs++;
        segRC->p1 = topMidPoint;
        segRC->p2 = bottomMidPoint;
        segRC->numEdges = abs(topMidPoint.y - bottomMidPoint.y);
        // get mid edges
        edgesRC = getIndex(topMidPoint,bottomMidPoint,rst);
        segRC->edges = (int *)malloc(segRC->numEdges*sizeof(int));
        segRC->edges = edgesRC;
        // create bottom segment
        segRC = &(routeRC.segments[routeRC.numSegs]);
        routeRC.numSegs++;
        segRC->p1 = bottomMidPoint;
        segRC->p2 = bottomPoint;
        segRC->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
        // get right edges
        edgesRC = getIndex(bottomMidPoint,bottomPoint,rst);
        segRC->edges = (int *)malloc(segRC->numEdges*sizeof(int));
        segRC->edges = edgesRC;
        //check every other possible RC route and compare weights
        if ((rightPoint.x + 2) <= (rst->gx - 1)) { /* may have to be -1 */
            for (int i = 2; i <= (rst->gx - (1 + rightPoint.x)); i++) {
                routeComp.numSegs = 0;
                // create midpoints
                topMidPoint.x = (rightPoint.x + i);
                topMidPoint.y = topPoint.y;
                bottomMidPoint.x = (rightPoint.x + i);
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
                // create right segment
                segComp = &(routeComp.segments[routeComp.numSegs]);
                routeComp.numSegs++;
                segComp->p1 = bottomMidPoint;
                segComp->p2 = bottomPoint;
                segComp->numEdges = abs(bottomMidPoint.x - bottomPoint.x);
                // get right edges
                edgesComp = getIndex(bottomMidPoint,bottomPoint,rst);
                segComp->edges = (int *)malloc(segComp->numEdges*sizeof(int));
                segComp->edges = edgesComp;

                // check if weight is lower than routeRC, assign to routeRC if true
                routeComp_weight = calcRouteCost(&routeComp, rst); 
                routeRC_weight = calcRouteCost(&routeRC, rst);
                if (routeComp_weight < routeRC_weight)
                    routeRC = routeComp;
            }
        }
        // return the route
        return routeRC;
    }
}

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
            rst->edgeOver = (int *) malloc(rst->numEdges * sizeof(int));
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
            //rst->nets = (net *)malloc(rst->numNets*sizeof(net));
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
                //printf("updating edgeCaps at: %d...\n", index[0]);
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
        net *tempNet = &(rst->nets.at(i));
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
                    //rst->edgeCaps[tempSeg->edges[i]]--;
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
                    //rst->edgeCaps[tempSeg->edges[i]]--;
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
                    //rst->edgeCaps[tempSeg->edges[i]]--;
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
    //cout << "Total Overflow: " << TOF << endl;
    //cout << "Total Wire Length: " << TWL << endl;
    return 1;
}

int calcEdgeWeights(int mode, route *newRoute, routingInst *rst){

    int* indecies;
    // calculate all edge weights at the beginning
    if(mode == 0){
        //printf("initializing edge weights\n");
        for(int i = 0; i < rst->numEdges; ++i){
            rst->edgeHis[i] = 0;
            rst->edgeOver[i] = max((rst->edgeUtils[i] - rst->edgeCaps[i]), 0);
            if(rst->edgeOver[i] > 0)
                rst->edgeHis[i] += 1;
            rst->edgeWeight[i] = rst->edgeOver[i] * rst->edgeHis[i];
            //printf("%d, ", rst->edgeWeight[i]);
        }

        //cout << endl << endl;
        return 1;
    }

    else if(mode == 1){
        //segments
        for(int i = 0; i < newRoute->numSegs; ++i){
            indecies = getIndex(newRoute->segments[i].p1, newRoute->segments[i].p2, rst);
            for(int j = 0; j < newRoute->segments[i].numEdges; ++j){
                if(rst->edgeOver[indecies[j]] > 0)
                    rst->edgeHis[indecies[j]] += 1;
                rst->edgeOver[indecies[j]] = max((rst->edgeUtils[indecies[j]] 
                            - rst->edgeCaps[indecies[j]]), 0);
                rst->edgeWeight[indecies[j]] = rst->edgeOver[indecies[j]] * 
                                               rst->edgeHis[indecies[j]];
            }
        }

        return 1;
    }
    return 0;
}

int calcNetCost(routingInst *rst){
    int* indecies;
    // calc cost
    // over all nets
    for(int i = 0; i < rst->numNets; ++i){
        rst->nets.at(i).cost = 0;
        // over all segments in the route
        for(int j = 0; j < rst->nets.at(i).nroute.numSegs; ++j){
            indecies = getIndex(rst->nets.at(i).nroute.segments[j].p1, rst->nets.at(i).nroute.segments[j].p2, rst);
            // over all edges in the segment
            for(int k = 0; k < rst->nets.at(i).nroute.segments[j].numEdges; ++k){
                rst->nets.at(i).cost += rst->edgeWeight[indecies[k]];
            }
        }
    }
    //sort nets
    sort(rst->nets.begin(), rst->nets.end(), compare);
    return 1;
}

int calcRouteCost(route *newRoute, routingInst *rst){
    int cost, *indecies;
    cost = 0;
    // over all segments in the route
    for(int i = 0; i < newRoute->numSegs; ++i){
        indecies = getIndex(newRoute->segments[i].p1, newRoute->segments[i].p2, rst);
        // over all edges in the segment
        for(int j = 0; j < newRoute->segments[i].numEdges; ++j){
            cost += rst->edgeWeight[indecies[j]];
        }
    }
    return cost;
}

int writeOutput(const char *outRouteFile, routingInst *rst){
	try {
        ofstream output;
        output.open(outRouteFile);

        for (int i = 0; i < rst->numNets; i++) {

            net tmpNet = rst->nets.at(i);
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
        printf("net n%d, has %d pins\n", i, rst->nets.at(i).numPins);
       /* printf("the pin connections are:\n");
        for(int j = 0; j < rst->nets.at(i).numPins; ++j){
            printf("\t(%d,%d)\n", rst->nets.at(i).pins[j].x,rst->nets.at(i).pins[j].y);
        }*/
    }
    return 1;
}


int release(routingInst *rst){
	try {
        for (int i = 0; i < rst->numNets; i++) {
            net deleteNet = rst->nets.at(i);
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

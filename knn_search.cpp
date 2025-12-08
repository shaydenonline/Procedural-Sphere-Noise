#include <ANN/ANN.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

int main(){
    std::cout << sizeof(int*); 
    constexpr int numPoints {215}; 
    constexpr int dimension {3};
    constexpr ANNcoord c = 0;
    double point[3] = {.5, .3, .2};
    
    ANNidxArray nnIdx = new ANNidx[3];
    ANNdistArray dists = new ANNdist[3];

    ANNpointArray mainPointer = annAllocPts(numPoints, dimension); 
    mainPointer[0] = point;
    std::cout << "*mainPointer " << mainPointer[0][2] << '\n';
    std::cout << "sizeof(mainPointer)" << sizeof(mainPointer[1]) << '\n';
    std::ifstream PolarCoords("polar_coords_107.txt");


    int i {0};
    while(i < numPoints){
        double x, y, z;
        PolarCoords >> x; 
        PolarCoords.seekg(2, std::ifstream::cur);
        PolarCoords >> y;
        PolarCoords.seekg(2, std::ifstream::cur);
        PolarCoords >> z;
        PolarCoords.seekg(1, std::ifstream::cur);
        double point[3] = {x,y,z};
        std::cout << pow(x,2) + pow(y,2) + pow(z,2) << '\n';
        ANNpoint pointPtr = annCopyPt(3, point);
        mainPointer[i++] = pointPtr;
    }

    ANNkd_tree kd_tree(mainPointer, numPoints, dimension);
    
    double testX, testY, testZ;
    double phi = {M_PI/7.0};
    double theta = {M_PI/9.0};

    testX = sin(phi)*cos(theta);
    testY = sin(phi)*sin(theta);
    testZ = cos(theta);
    double queryArray[3] = {testX, testY, testZ};

    ANNpoint queryPoint = annCopyPt(3, queryArray);
    
    kd_tree.ANNkd_tree::annkSearch(queryPoint, 3, nnIdx, dists);
    for(int i {0}; i < 3; i++){
        for(int j {0}; j < 3; j++){
            std::cout << "nnIdx[" << i << ']' << '[' << j << "]: " << mainPointer[nnIdx[i]][j] << '\n';
        }
    }
    delete [] nnIdx;
    delete [] dists;
    //delete &kd_tree;
    annClose();
    return 0;
}



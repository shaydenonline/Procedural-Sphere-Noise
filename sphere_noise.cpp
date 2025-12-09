#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <algorithm>
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <chrono>
#include <random>
#include <ANN/ANN.h>

//Included for wall-clock timing certain methods
//Since all other namespaces are referred to specifically, naming conflicts aren't likely to occur
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;



constexpr double goldenRatio {(1 + pow(5, .5)) / 2};
constexpr double goldenAngle {2 * M_PI * (1 - (1/goldenRatio))};

constexpr double thetaRangeRadians {M_PI};
constexpr double phiRangeRadians {M_PI/2.0};

std::string convertIntToChar(int num);
std::vector<double> convertPolarToCartesian(const std::vector<double> &polarCoord);
std::vector<double> convertPolarToCartesian(const double theta, const double phi, const double rho = 1);
std::vector<double> crossProduct(const std::vector<double> & u, const std::vector<double> &v);
std::string createFibGrid(int points);
bool parseDoubles(const std::string& line, std::vector<double>& out);

double geodesic(const std::vector<double> &u, const std::vector<double> &v);
double dotProd(std::vector<double> point1, std::vector<double> point2, int entries = 3);
double distance(std::vector<double> point1, std::vector<double> point2);
double length(const std::vector<double>& point1);
double gaussianRBF(const std::vector<double>& u, const std::vector<double>& v, double epsilon = .5);
bool checkIfUnitVector(double x, double y, double z, double errorMargin = 1e-4);
//
//
//
//Use this to shuffle the permutation table entries, use something else for the gradient entries.
class SeededRandom {
private:
    long seed;
    static const long a = 48271;
    static const long c = 0;
    static const long m = 2147483647;

public:
    SeededRandom(int seed) : seed(seed) {}

    int Next() {
        seed = (a * seed + c) % m;
        return static_cast<int>(seed);
    }

    float NextFloat() {
        return Next() / static_cast<float>(m);
    }
};

class PerlinNoise {
private:
    std::vector<int> permutationTable;
    std::vector<std::vector<double>> gradients;
    static const int tableSize = 256;
    static const int tableSizeMask = tableSize - 1;
    SeededRandom seededRandom;
    ANNkd_tree* kd_tree;
    ANNpointArray mainPointer;
    //Mistake in this method we don't want it to modify point1 but ATM we're not using it.
    void slerp(std::vector<double>& point1, const std::vector<double>& point2, std::vector<double>& resPoint, double t){
        double dotProd {point1[0]*point2[0] + point1[1]*point2[1] + point1[2]*point2[2]};
        if(dotProd < 0.0){
            dotProd = -dotProd;
            point1[0] *= -1; point1[1] *= -1; point1[2] *= -1;
        }
        double omega = std::acos(dotProd);
        if (std::abs(omega) < 1e-8) {
            resPoint[0] = (1 - t)*point1[0] + t*point2[0];
            resPoint[1] = (1 - t)*point1[1] + t*point2[1];
            resPoint[2] = (1 - t)*point1[2] + t*point2[2];
            return;
        }
        double a {(std::sin((1-t)*omega))/std::sin(omega)};
        double b {std::sin(t*omega)/std::sin(omega)};
        resPoint[0] = a*point1[0] + b*point2[0];
        resPoint[1] = a*point1[1] + b*point2[1];
        resPoint[2] = a*point1[2] + b*point2[2];
    }
    
    //This method should be called once to initialize our gradients. The seed can either be non-deterministic ie std::random_device or deterministic (will be the same sequence each runtime)
    //For now we're only using one seed for both phi and theta, could alter to use two seeds if we suspect there's an issue
    std::vector<std::vector<double>> generateGradients(int seed, const int numGradients) {
        constexpr int cols {3};
        std::vector<std::vector<double>> gradients(numGradients, std::vector<double>(cols));


        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> theta_dist(-thetaRangeRadians, thetaRangeRadians);
        std::uniform_real_distribution<double> phi_dist(-phiRangeRadians, phiRangeRadians);
        
        for(int i {0}; i < numGradients; i++){
            double theta {theta_dist(gen)};
            double phi {phi_dist(gen)};
            gradients[i] = convertPolarToCartesian(theta, phi);
            assert(checkIfUnitVector(gradients[i][0], gradients[i][1], gradients[i][2]) && "We didn't generate gradients that are unit vectors");
        }

        return gradients;
    }

public:
    PerlinNoise(int seed) : seededRandom(seed) {
        gradients = generateGradients(12, 12);
        assert(tableSize == 256 && "we don't have access to table size\n");
        permutationTable.resize(tableSize * 2);
        std::vector<int> tempTable(tableSize);
        for (int i = 0; i < tableSize; i++)
            tempTable[i] = i;
        for (int i = 0; i < tableSize; i++) {
            int j = seededRandom.Next() % tableSize;
            std::swap(tempTable[i], tempTable[j]);
        }
        for (int i = 0; i < tableSize; i++) {
            permutationTable[i] = tempTable[i];
            permutationTable[tableSize + i] = tempTable[i];
        }
        static const int tableSizeMask = tableSize - 1;

        //Building our data structure to conduct k-nn searches we will generate noise for about 10,000 points on the sphere with 215 grid points (generated by the fibonacci grid)
        constexpr int numPoints {303};
        constexpr int dimension {3};
        constexpr ANNcoord c = 0;
        mainPointer = annAllocPts(numPoints, dimension);
        std::ifstream FibGridCoords("polar_coords_151.txt");
        int i {0};
       
        //Code to initialize the kd_tree of grid points from the file of fibonacci grid points. TESTED
        while(i < numPoints){
            std::vector<double> vals;
            std::string line;
            std::getline(FibGridCoords, line);
            
            if(!parseDoubles(line, vals)){
            std::cerr << "Error invalid input tokens \n";
            }
            assert(vals.size() == 3 && "vals.size() is not 3 \n");

            mainPointer[i][0] = vals[0];
            mainPointer[i][1] = vals[1];
            mainPointer[i][2] = vals[2];
            i++;
        }
        assert(mainPointer != nullptr && "mainPointer is not getting initialized properly");
        assert(mainPointer[0][0] != mainPointer[2][0] &&  "mainPointer is repeating initialization values");
        kd_tree = new ANNkd_tree(mainPointer, numPoints, dimension);
    }
    
    //Gets nearest K neighbors, default value is 3 neighbors. TESTED
    void getNearestKNeighbors(const std::vector<double>& targetPoint, std::vector<int>& indices, const int K = 3){
        //Dynamically allocate space inside method for the number of neighbors you want (distances which we won't be using and indices)
        ANNidxArray nnIdx = nullptr;
        ANNdistArray dists = nullptr;
        try {
            nnIdx = new ANNidx[K];
            dists = new ANNdist[K];
        }
        catch (const std::bad_alloc&) {
            std::cerr << "Memory allocation failed\n";
        }

        //We are only working with sphere, 3 coordinates and the special ANN library can only take in a vector so that's why we have to do all this complicated point rigamarole 
        double target[3] = {targetPoint[0], targetPoint[1], targetPoint[2]};
        assert(checkIfUnitVector(target[0], target[1], target[2]) && "This point doesn't have unit length, it's likely you're retrieving garbage memory");

        ANNpoint queryPoint = annCopyPt(3, target);

        kd_tree->ANNkd_tree::annkSearch(queryPoint, K, nnIdx, dists);

        
        //After we obtain the indices of the nearest neighbors put them in vector 
        indices.resize(K);
        for(int i {0}; i < K; i++) indices[i] = nnIdx[i];
        assert(indices[0] != indices[1] && "At least one pair of nearest neighbor indices aren't distinct");
       
        //cleanup
        annDeallocPt(queryPoint);
        delete [] nnIdx;
        delete [] dists;
    }

    double Noise(double x, double y, double z) {
        std::vector<double> targetPoint {x,y,z};
        std::vector<int> indices;
        getNearestKNeighbors(targetPoint, indices);
        std::vector<std::vector<double>> nearestNeighborGradients(indices.size(), std::vector<double>(3));
        std::vector<double> weightedDotProds(indices.size());
        
        double total {0};
        for(int i {0}; i < indices.size(); i++){
            std::vector<double> randomGrad = gradients[permutationTable[indices[i] % 256] % gradients.size()];
            //The default parameter for the gaussianRBF is set to 1. The points are already quite close to each other. a smaller epsilon will sphread the value of RBF closer to value 1, a larger value for epsilon will lead to values tending to 0. any value less than 1 risks a loss in precision for now but this would be good to try and adjust
            double weight {gaussianRBF(targetPoint, {mainPointer[i][0], mainPointer[i][1], mainPointer[i][2]})};
            double dot {dotProd({targetPoint[0] - mainPointer[i][0], targetPoint[1] - mainPointer[i][1], targetPoint[2] - mainPointer[i][2]}, randomGrad)};
            //double dot {dotProd(crossProduct(targetPoint, {mainPointer[i][0], mainPointer[i][1], mainPointer[i][2]}), randomGrad)};
            total += weight*dot;
        }
        return std::abs(total) > 1 ? std::abs(total) - 1 : std::abs(total); 
        //the most commonly used normalization method for perlin noise

    }
    void printTable(){
        int n = permutationTable.size();
        for(int i {0}; i < n; i++){
            std::cout << permutationTable[i] << '\n';
        }
    }
    //The ANN library has a lot of opportunities for memory leaks. Make sure you 1st delete the kd_tree then deallocate the points you put in it. Otherwise there'll be a memory leak.
    ~PerlinNoise() {
        delete kd_tree;
        annDeallocPts(mainPointer);
        annClose();
    }
};


double gaussianRBF(const std::vector<double>& u, const std::vector<double>& v, double epsilon){
    double r {geodesic(u,v)};
    double exponent {r * epsilon};
    exponent *= exponent;
    return std::exp(-exponent);
}


std::string createFibGrid(int points){

    std::vector<std::vector<double>> latLong {};
    double latitude {0}; double longitude {0};

    char pointsChar = points + '0';
    std::string fileName = "polar_coords_" + convertIntToChar(points) + ".txt";
    std::ofstream PolarCoords(fileName);

    for(int i = -points; i <= points; i++){
        latitude = asin((2.0*i)/(2.0*points + 1.0)) * (180 / M_PI);
        longitude = fmod(i * goldenAngle * (180.0 / M_PI) + 180.0, 360.0) - 180.0;
        if(longitude < -180) longitude += 360;
        if(longitude > 180) longitude -= 360;
        assert(("Latitude is not in the range [-90,90] degrees", latitude >= -90.0 && latitude <= 90.0));
        assert(("Longitude is not in the range [-180, 180] degrees", longitude >= -180.0 && longitude <=180.0));
        latLong.push_back({latitude, longitude});
    }

    std::vector<std::vector<double>> cartesianCoords(latLong.size());
    for(int i = 0; i < latLong.size(); i++){
        cartesianCoords[i] = convertPolarToCartesian(latLong[i]);
        assert(("There should be three coordinates", cartesianCoords[i].size() == 3));
    }
    for(auto &point: cartesianCoords){
        PolarCoords << point[0] << ' ' << point[1] << ' ' << point[2] << '\n';
    }
    PolarCoords.close();
    return fileName;
}

std::string convertIntToChar(int num){
    std::string stringRepresentation;
    while(num > 0){
        if(num % 10 != 0){
            stringRepresentation.push_back(num % 10 + '0');
            num -= (num % 10);
            num /= 10;
        }
        else{
            stringRepresentation.push_back('0');
            num /= 10;
        }
    }
    std::reverse(stringRepresentation.begin(), stringRepresentation.end());
    return stringRepresentation;
}

std::vector<double> convertPolarToCartesian(const std::vector<double> &polarCoord){
    double latitude = polarCoord[0] * (M_PI / 180);
    double longitude = polarCoord[1] * (M_PI / 180);
    double x = cos(latitude) * cos(longitude);
    double y = cos(latitude) * sin(longitude);
    double z = sin(latitude);
    return {x,y,z};
}

std::vector<double> convertPolarToCartesian(const double theta, const double phi, const double rho){
    assert(("Latitude is not in the range [-90,90] degrees", phi >= -90.0 && phi <= 90.0));
    assert(("Longitude is not in the range [-180, 180] degrees", theta >= -180.0 && theta <=180.0));

    double x = rho* cos(theta) * sin(phi);
    double y = rho* sin(theta) * sin(phi);
    double z = rho * cos(phi);
    return {x,y,z};
}
double dotProd(std::vector<double> point1, std::vector<double> point2, int entries){
    assert(point1.size() == point2.size() && point1.size() == entries && "Two vectors you're taking dot product of aren't the same size!\n");
    double res {0};
    for(int i {0}; i < entries; i++){
        res += point1[i]*point2[i];
    }
    return res;
}

double distance(std::vector<double> point1, std::vector<double> point2){
    assert(point1.size() == 3 && point2.size() == 3 && "We are only calculating distance between 3D vectors");
    double cosang {dotProd(point1, point2, 3)};
    double distance {acos(cosang)};
    return distance;
}

double length(const std::vector<double>& point1){
    return std::sqrt(point1[0]* point1[0] + point1[1]* point1[1] + point1[2]* point1[2]);
}

bool parseDoubles(const std::string& line, std::vector<double>& out){
    std::stringstream ss(line);
    double val;
    while(ss >> val){
        out.push_back(val);
    }
    if(!ss.eof()){
        return false;
    }
    return true;
}

std::vector<double> crossProduct(const std::vector<double> & u, const std::vector<double> &v){
    std::vector<double> result(3);
    result[0] = u[1]*v[2] - u[2]*v[1];
    result[1] = u[2]*v[0] - u[0]*v[2];
    result[2] = u[0]*v[1] - v[0]*u[1];
    return result;
}

double geodesic(const std::vector<double> &u, const std::vector<double> &v){
    return atan2(length(crossProduct(u,v)), dotProd(u,v));
}
bool checkIfUnitVector(double x, double y, double z, double errorMargin){
    double squaredNorm = x*x + y*y + z*z;
    if(!((1 - squaredNorm) < errorMargin)) std::cout << "not a unit vector: " << squaredNorm << '\n'; 
    return (1 - squaredNorm) < errorMargin;
}
int main() {
   
    std::cout << "Enter a number n to generate fibonacci grid: ";
    int points {0};
    std::cin >> points;
    if(points % 2 == 0) points++;

    std::string theFileName = createFibGrid(points);
     
    int seed = 12345;
    PerlinNoise perlinNoise(seed);
    
    std::ifstream NoisePoints(theFileName);
    std::ofstream NoiseValues("noise.txt");
    int i {0};
    double x, y, z;
    double value;
    while(i <= points*2){
        std::string line;
        std::getline(NoisePoints, line);
        std::stringstream ss(line);
        std::vector<double> vals;
        if(!parseDoubles(line, vals)) { 
            std::cerr << "Error: input contains invalid tokens. \n";
            return 1;
        }
        assert(vals.size() == 3 && "vals should be 3\n");
        x = vals[0]; y = vals[1]; z = vals[2];
        double myNoise = perlinNoise.Noise(x,y,z);
        NoiseValues << myNoise << '\n';
        i++;
    }
    NoisePoints.close(); NoiseValues.close();
    return 0;
}



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
#include <chrono>
#include <random>
#include <ANN/ANN.h>

    

constexpr double goldenRatio {(1 + pow(5, .5)) / 2};
constexpr double goldenAngle {2 * M_PI * (1 - (1/goldenRatio))};

std::string convertIntToChar(int num);
std::vector<double> convertPolarToCartesian(const std::vector<double> &polarCoord);
std::string createFibGrid();
double dotProd(std::vector<double> point1, std::vector<double> point2, int entries);
double distance(std::vector<double> point1, std::vector<double> point2);


class SeededRandom {
private:
  // Current seed value, evolves with each call to Next() or NextFloat().
  long seed;

  // The multiplier 'a'. This value is part of the parameters that define the quality and characteristics
  // of the linear congruential generator. 48271 is a commonly used value in LCG algorithms for its good properties.
  static const long a = 48271;

  // The increment 'c'. It is set to 0 in this implementation, defining it as a multiplicative LCG.
  static const long c = 0;

  // The modulus 'm'. Using 2^31 - 1 (a Mersenne prime) as the modulus helps in achieving a full period
  // for the generated pseudo-random sequence, maximizing the sequence's length before it repeats.
  static const long m = 2147483647;

public:
  /// Constructor that initializes the pseudo-random number generator with a seed value.
  /// @param seed Initial seed value for generating numbers.
  SeededRandom(int seed) : seed(seed) {}

  /// Generates the next number in the sequence of pseudo-random numbers.
  /// @return The next pseudo-random number as an int.
  int Next() {
    // Calculates the next seed value using the linear congruential generator formula.
    // The use of static_cast<int> ensures the result is properly typed as an integer.
    seed = (a * seed + c) % m;
    return static_cast<int>(seed);
  }

  /// Generates a floating-point number between 0 (inclusive) and 1 (exclusive) based on the pseudo-random sequence.
  /// @return A pseudo-random float between 0 and 1.
  float NextFloat() {
    // Utilizes the Next() function to obtain a pseudo-random integer, then divides it by 'm'
    // to normalize the result to a floating-point number in the range [0, 1).
    return Next() / static_cast<float>(m);
  }
};
// Existing implementation

// A class for generating Perlin noise, a procedural technique used in computer graphics to produce natural-looking textures.
class PerlinNoise {
private:
  // A table used for permutation in the Perlin noise algorithm. It's doubled to avoid overflow in indexing.
  std::vector<int> permutationTable;
 std::vector<std::vector<double>> gradients;
  // The size of the permutation table.
  static const int tableSize = 256;
  // A bitmask used for wrapping indices to remain within the permutation table's bounds.
  static const int tableSizeMask = tableSize - 1;
  // Seeded random number generator to shuffle the permutation table in a reproducible way.
  SeededRandom seededRandom;
    ANNkd_tree* kd_tree;
  // The Fade function smooths the input values to ease transitions, as described by Ken Perlin.
  float Fade(float t) { return t * t * t * (t * (t * 6 - 15) + 10); }

  // Linear interpolation between two values a and b using the blend factor t.
  float Lerp(float a, float b, float t) { return a + t * (b - a); }

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
    /*
    std::vector<Vec3> gradients;
    for (int i = 0; i < 64; i++) {
        double z = uniform(-1, 1);
        double t = uniform(0, 2*M_PI);
        double r = sqrt(1 - z*z);
        gradients.push_back({r*cos(t), r*sin(t), z});
    }
    */


std::vector<std::vector<double>> generate_gradients() {
    std::vector<std::vector<double>> gradients;
    gradients.reserve(64);

    std::random_device rd;               // seeds RNG from OS
    std::mt19937 gen(rd());              // Mersenne Twister RNG
    std::uniform_real_distribution<double> dist01(-1.0, 1.0);
    std::uniform_real_distribution<double> distAngle(0.0, 2.0 * M_PI);

    for (int i = 0; i < 64; ++i) {
        double z = dist01(gen);                 // uniform in [-1, 1]
        double theta = distAngle(gen);          // uniform angle in [0, 2π]
        double r = std::sqrt(1.0 - z*z);

        gradients.push_back({
            r * std::cos(theta),
            r * std::sin(theta),
            z
        });
    }

    return gradients;
}

public:
  // Constructor initializes the Perlin noise generator with a specific seed.
  PerlinNoise(int seed) : seededRandom(seed) {
    gradients = generate_gradients();
    permutationTable.resize(tableSize * 2);
    std::vector<int> tempTable(tableSize);
    for (int i = 0; i < tableSize; i++)
      tempTable[i] = i;

    // Shuffle the temporary permutation table using the seeded random generator.
    for (int i = 0; i < tableSize; i++) {
      int j = seededRandom.Next() % tableSize;
      std::swap(tempTable[i], tempTable[j]);
    }

    // Duplicate the shuffled permutation table to avoid overflow.
    for (int i = 0; i < tableSize; i++) {
      permutationTable[i] = tempTable[i];
      permutationTable[tableSize + i] = tempTable[i];
    }


    constexpr int numPoints {215}; 
    constexpr int dimension {3};
    constexpr ANNcoord c = 0;
    
    
    ANNpointArray mainPointer = annAllocPts(numPoints, dimension); 
    std::ifstream FibGridCoords("polar_coords_107.txt");
    
    
    int i {0};
    while(i < numPoints){
        double x, y, z;
        FibGridCoords >> x; 
        FibGridCoords.seekg(2, std::ifstream::cur);
        FibGridCoords >> y;
        FibGridCoords.seekg(2, std::ifstream::cur);
        FibGridCoords >> z;
        FibGridCoords.seekg(1, std::ifstream::cur);
        double point[3] = {x,y,z};
        std::cout << pow(x,2) + pow(y,2) + pow(z,2) << '\n';
        ANNpoint pointPtr = annCopyPt(3, point);
        mainPointer[i++] = pointPtr;
        annDeallocPt(pointPtr);
    }
    
    kd_tree = new ANNkd_tree(mainPointer, numPoints, dimension);
    
 
  }
    void getNearest3Neighbors(const std::vector<double>& targetPoint, std::vector<int>& indices){
       ANNidxArray nnIdx = new ANNidx[3];
       ANNdistArray dists = new ANNdist[3];

       
       double target[3] = {targetPoint[0], targetPoint[1], targetPoint[2]};
       ANNpoint queryPoint = annCopyPt(3, target);
       
       kd_tree->ANNkd_tree::annkSearch(queryPoint, 3, nnIdx, dists);
       indices.resize(3);
       indices[0] = nnIdx[0]; indices[1] = nnIdx[1]; indices[2] = nnIdx[2];
        annDeallocPt(queryPoint);
       delete [] nnIdx;
       delete [] dists;
    }  
  // Generates Perlin noise for a given point (x, y, z) in 3D space.
  float Noise(float x, float y, float z) {
        std::vector<double> targetPoint {x,y,z};
        std::vector<int> indices;
        getNearest3Neighbors(targetPoint, indices);
        std::vector<double> G1 = gradients[indices[0] % 64];
        std::vector<double> G2 = gradients[indices[1] % 64];
        std::vector<double> G3 = gradients[indices[2] % 64];
        double d1 = distance(G1, targetPoint);
        double d2 = distance(G2, targetPoint);
        double d3 = distance(G3, targetPoint);
        constexpr int alpha {3};
        double w1 = pow(1/d1, alpha);
        double w2 = pow(1/d2, alpha);
        double w3 = pow(1/d3, alpha);
        double norm = w1 + w2 + w3;
        w1 /= norm;
        w2 /= norm;
        w3 /= norm;
        return w1*d1 + w2*d2 + w3*d3;

        /*
        vector<double> G1 {kd_tree[indices[0]][0], kd_tree[indices[0]][1], kd_tree[indices[0]][2]}; 
 vector<double> G2 {kd_tree[indices[1]][0], kd_tree[indices[1]][1], kd_tree[indices[1]][2]}; 
 vector<double> G3 {kd_tree[indices[2]][0], kd_tree[indices[2]][1], kd_tree[indices[2]][2]}; 
*/
  }
};

int main() {

    /*
  	int seed = 12345;
  	PerlinNoise perlinNoise(seed);
	float x {0.0f}; float y {0.0f};
	float z {0.0f};
	std::ofstream PointsFile("test_points.txt");
	while(x < 1.0f){
		while(y < 1.0f){
			PointsFile << x << ',' << y  << ',' << perlinNoise.Noise(x,y,z+=0.1f)/10.0f << '\n';
			y += .01f;
		}
		y = 0.0f;
		x += .01f;
	}
	PointsFile.close();
*/
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();
    std::string theFileName = createFibGrid();
    auto t2 = high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cout << "Filename: " << theFileName << '\n';
    std::cout << "Time taken: " << ms_int.count() << '\n';

	/*
  float x = 1.5f;
  float y = 1.5f;
  float z = 0.0f;
	
  float noiseValue =
      perlinNoise.Noise(x, y, z);
  std::cout << "Noise Value: " << noiseValue << std::endl;
	*/
	/*
	float noiseValue;
	for(int y = 0; y < imageDimension* imageDimension; y++){
     		noiseValue = perlinNoise.Noise(x,y,0);	
		image[y] = (unsigned char)((noiseValue + 1.0f) * 127.5f);
	}
	*/
  return 0;
}


std::string createFibGrid(){
    std::cout << "Enter a number n to generate fibonacci grid: ";
    int points {0};
    std::cin >> points;
    if(points % 2 == 0) points++;

    std::vector<std::vector<double>> latLong {};
    double latitude {0}; double longitude {0};

    char pointsChar = points + '0';
    std::string fileName = "polar_coords_" + convertIntToChar(points) + ".txt";
    std::ofstream PolarCoords(fileName);
    std::cout << "goldenAngle: " << goldenAngle << '\n';

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
        PolarCoords << point[0] << ", " << point[1] << ", " << point[2] << '\n';
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

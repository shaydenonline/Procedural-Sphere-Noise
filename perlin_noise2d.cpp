#define STB_IMAGE_WRITE_IMPLEMENTATION  
#include "stb_image_write.h"  
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
/// A pseudo-random number generator class based on a linear congruential generator (LCG) algorithm.
/// This class allows generating a sequence of pseudo-random numbers from a seed value,
/// ensuring reproducible outcomes which are useful for simulations, testing, and gaming.
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
  // The size of the permutation table.
  static const int tableSize = 256;
  // A bitmask used for wrapping indices to remain within the permutation table's bounds.
  static const int tableSizeMask = tableSize - 1;
  // Seeded random number generator to shuffle the permutation table in a reproducible way.
  SeededRandom seededRandom;

  // The Fade function smooths the input values to ease transitions, as described by Ken Perlin.
  float Fade(float t) { return t * t * t * (t * (t * 6 - 15) + 10); }

  // Linear interpolation between two values a and b using the blend factor t.
  float Lerp(float a, float b, float t) { return a + t * (b - a); }

  // Calculates a gradient based on a hash value and the x, y, z coordinates. 
  // This contributes to the pseudo-randomness of the noise.
  float Grad(int hash, float x, float y, float z) {
    int h = hash & 15;
    float u = h < 8 ? x : y;
    float v = h < 4 ? y : h == 12 || h == 14 ? x : z;
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
  }

public:
  // Constructor initializes the Perlin noise generator with a specific seed.
  PerlinNoise(int seed) : seededRandom(seed) {
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
  }

  // Generates Perlin noise for a given point (x, y, z) in 3D space.
  float Noise(float x, float y, float z) {
    // Find the unit cube containing the point and wrap the integer parts to the table size.
    int X = static_cast<int>(std::floor(x)) & tableSizeMask;
    int Y = static_cast<int>(std::floor(y)) & tableSizeMask;
    int Z = static_cast<int>(std::floor(z)) & tableSizeMask;

    // Calculate the relative position of the point in the cube.
    x -= std::floor(x);
    y -= std::floor(y);
    z -= std::floor(z);

    // Compute fade curves for x, y, z.
    float u = Fade(x);
    float v = Fade(y);
    float w = Fade(z);

    // Hash coordinates of the cube's eight corners.
    int A = permutationTable[X] + Y;
    int AA = permutationTable[A] + Z;
    int AB = permutationTable[A + 1] + Z;
    int B = permutationTable[X + 1] + Y;
    int BA = permutationTable[B] + Z;
    int BB = permutationTable[B + 1] + Z;

    // Add blended results from the eight corners of the cube.
    float res =
        Lerp(w,
             Lerp(v,
                  Lerp(u, Grad(permutationTable[AA], x, y, z),
                       Grad(permutationTable[BA], x - 1, y, z)),
                  Lerp(u, Grad(permutationTable[AB], x, y - 1, z),
                       Grad(permutationTable[BB], x - 1, y - 1, z))),
             Lerp(v,
                  Lerp(u, Grad(permutationTable[AA + 1], x, y, z - 1),
                       Grad(permutationTable[BA + 1], x - 1, y, z - 1)),
                  Lerp(u, Grad(permutationTable[AB + 1], x, y - 1, z - 1),
                       Grad(permutationTable[BB + 1], x - 1, y - 1, z - 1))));

    // Normalize the result to be within the range [0, 1].
    return (res + 1.0f) / 2.0f;
  }
};

int main() {
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

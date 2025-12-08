#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cassert>


constexpr double goldenRatio {(1 + pow(5, .5)) / 2};
constexpr double goldenAngle {2 * M_PI * (1 - (1/goldenRatio))};
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
int main(){
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
        PolarCoords << point[0] << ' ' << point[1] << ' ' << point[2] << '\n';
    }
}

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

    std::vector<std::vector<double>> generate_gradients() {
        constexpr int rows {12};
        constexpr int cols {3};
        std::vector<std::vector<double>> gradients(rows, std::vector<double>(cols));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist01(-1.0, 1.0);
        std::uniform_real_distribution<double> distAngle(0.0, 2.0 * M_PI);
        for (int i = 0; i < 12; ++i) {
            double z = dist01(gen);
            double theta = distAngle(gen);
            double r = std::sqrt(1.0 - z*z);
            gradients[i] = {
                r * std::cos(theta),
                r * std::sin(theta),
                z
            };
            double norm = gradients[i][0]*gradients[i][0] + gradients[i][1]*gradients[i][1] + gradients[i][2]*gradients[i][2];
            assert(std::abs(1 - norm) < 1e-6 && "gradient vector doesn't have correct norm\n");
        }
        return gradients;
    }
int main(){
    std::vector<std::vector<double>> gradients = generate_gradients();
    for(int i {0}; i <gradients.size(); i++){
        for(int j {0}; j <gradients[0].size(); j++){
            std::cout << gradients[i][j] << ',';
        }
        std::cout << '\n';
    }
    return 0;
}

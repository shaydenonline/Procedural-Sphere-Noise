#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

void slerp(std::vector<double>& point1, const std::vector<double>& point2, std::vector<double>& resPoint, double t){
    double dotProd {point1[0]*point2[0] + point1[1]*point2[1] + point1[2]*point2[2]};
    if(dotProd < 0.0){
        dotProd = -dotProd;
        point1[0] *= -1; point1[1] *= -1; point1[2] *= -1;
    }
    
    double omega = std::acos(dotProd);
    
    if (std::abs(omega) < 1e-8) {
        out[0] = (1 - t)*p1[0] + t*p2[0];
        out[1] = (1 - t)*p1[1] + t*p2[1];
        out[2] = (1 - t)*p1[2] + t*p2[2];
        return;
    }

    double a {(std::sin((1-t)*omega))/std::sin(omega)};
    double b {std::sin(t*omega)/std::sin(omega)};

    resPoint[0] = a*point1[0] + b*point2[0];
    resPoint[1] = a*point1[1] + b*point2[1];
    resPoint[2] = a*point1[2] + b*point2[2];
}


void slerp(const std::vector<double>& p1_in,
           const std::vector<double>& p2_in,
           std::vector<double>& out,
double t)
{
    // ---- 1. Normalize inputs ----
    auto norm = [](const std::vector<double>& v) {
        double m = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return std::vector<double>{v[0]/m, v[1]/m, v[2]/m};
    };

    std::vector<double> p1 = norm(p1_in);
    std::vector<double> p2 = norm(p2_in);

    // ---- 2. Dot product, clamped ----
    double dot = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
    dot = std::min(1.0, std::max(-1.0, dot));

    // ---- 3. Angle between vectors ----
    double omega = std::acos(dot);

    // ---- 4. Handle nearly parallel vectors (avoid zero division) ----
    if (std::abs(omega) < 1e-8) {
        out[0] = (1 - t)*p1[0] + t*p2[0];
        out[1] = (1 - t)*p1[1] + t*p2[1];
        out[2] = (1 - t)*p1[2] + t*p2[2];
        return;
    }

    double sin_omega = std::sin(omega);

    double a = std::sin((1 - t) * omega) / sin_omega;
    double b = std::sin(t * omega) / sin_omega;

    out[0] = a*p1[0] + b*p2[0];
    out[1] = a*p1[1] + b*p2[1];
    out[2] = a*p2[2] + b*p2[2];
}
void simpleSlerp(std::vector<double>& point1, const std::vector<double>& point2, std::vector<double>& resPoint, const double t){
    const double a = 1 - t;
    const double b = t;
    resPoint[0] = a*point1[0] + b*point2[0];
    resPoint[1] = a*point1[1] + b*point2[1];
    resPoint[2] = a*point1[2] + b*point2[2];
}
int main(){
    
    std::vector<double> point1 {-0.00267475, -0.0544821, -0.998511};
    std::vector<double> point2 {-0.221592, 0.240907, -0.944913};
    std::vector<double> resPoint {0.0, 0.0, 0.0};
    constexpr double t {.5};
    simpleSlerp(point1, point2, resPoint, t);
    std::cout << resPoint[0] << ", " << resPoint[1] << ", " << resPoint[2] << '\n';
    std::cout << point1[0] << ", " << point1[1] << ", " << point1[2] << '\n';
    std::cout << point2[0] << ", " << point2[1] << ", " << point2[2] << '\n';
    slerp(point1, point2, resPoint, t);
    std::cout << resPoint[0] << ", " << resPoint[1] << ", " << resPoint[2] << '\n';
    std::cout << point1[0] << ", " << point1[1] << ", " << point1[2] << '\n';
    std::cout << point2[0] << ", " << point2[1] << ", " << point2[2] << '\n';

    return 0;
}

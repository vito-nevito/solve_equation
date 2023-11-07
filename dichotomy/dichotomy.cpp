#include <iostream>
#include <array>
#include <cmath>
#include <iomanip>
#include <fstream>

double func(double x){
    double res = exp(-pow(x-3,2))+log(1+x)-0.5*x;
    return res;
}

void change_points(double &a, double &b, double sep_point){
    if (func(a) * func(sep_point) <= 0){
        b = sep_point;
    }
    else if(func(sep_point) * func(b) <= 0){
        a = sep_point;
    }
}

int main(){
    double a, b;
    a = 0.0;
    b = 10.0;
    std::array<double,3> eps_array = {1.0e-3, 1.0e-6, 1.0e-9};
    double current_eps;
    double sep_point, value;
    std::cout<<std::setprecision(9);
    std::ofstream fout("dichotomy_data.txt");
    for (const auto eps : eps_array){
        current_eps = b - a;
        while (current_eps >= eps){
            sep_point = (a + b)/2;
            current_eps /= 2.0;
            change_points(a, b, sep_point);
        }
        value = func(sep_point);
        std::cout<<eps<<' '<<sep_point<< ' ' <<value<< std::endl;
        fout<<eps<<' '<<sep_point;
    }
    fout.close();
    return 0;
}
#include<iostream>
#include<math.h>
#include<array>
#include <windows.h>

double target(double x)
{
    return 4*sin(x/2) + cos(x) * tanh(x) - x + 2;
}

double Chorda(double a, double b, double eps, int& n)
{
    double x = a;
    double x_last = b;
    n = 0;
    while(abs(x - x_last) > eps)
    {
     n++;
     x_last = x;
     x = x - (x - b) * target(x) / (target(x) - target(b));
    }
    return x;
}

int main()
{
    SetConsoleCP (CP_UTF8) ;
    SetConsoleOutputCP (CP_UTF8) ;
    double a = 0.;
    double b = 10.;
    double eps = 0.000001;
    int n = 0;
    double x = Chorda(a, b, eps, n);
    std::cout << "Метод хорд"<< std::endl;
    std::cout << "Корень: "<< x << std::endl;
    std::cout << "Количество итераций: "<< n << std::endl;
    std::cout << "Выполнение второго критерия: ";
    if(abs(target(x)) < eps)
        std::cout << "True" << std::endl;
    else
        std::cout << "False" << std::endl;
}



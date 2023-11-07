#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#define A 0
#define B 10
#define EPS1 1e-3
#define EPS2 1e-6
#define EPS3 1e-9

double func(double x) { return  4*sin(x/2) + cos(x)*tanh(x) - x + 2; };
double dfunc(double x) { return  2*cos(x/2) - sin(x)*tanh(x) + cos(x) / (pow(cosh(x), 2)) - 1; };
double ddfunc(double x) { return  -sin(x/2) - cos(x)*tanh(x) - 2*sin(x) / (pow(cosh(x), 2)) - 2*cos(x)*sinh(x)/ (pow(cosh(x), 3)); };
typedef double(*function)(double x);

double getStartPosition(function f, function ddf, double xLeft, double xRight);
double newtonSolve(function f, function df, double epsilon, double x0, std::string filename);

bool breakCheck1(double xnew, double xold, double epsilon1);
bool breakCheck2(function f, double xnew, double epsilon2);

void setInterval(double* xl, double* xr);

int main()
{
    double xl(A), xr(A);
    setInterval(&xl, &xr);
    //xl, xr = A, B;
    double x0 = getStartPosition(func, ddfunc, xl, xr);
    //x0 = 1.1;
    std::ofstream out;

    double root(0);
    root = newtonSolve(func, dfunc, EPS1, x0, "1e-3.txt");
    std::cout << "x = " << root << std::endl;
    std::cout << std::setprecision(16) << "x0 = " << x0 << std::endl;
    std::cout << "[" << xl << "; " << xr << "]" << std::endl;
    std::cout << "e = " << EPS1 << std::endl << std::endl;

    root = newtonSolve(func, dfunc, EPS2, x0, "1e-6.txt");
    std::cout << "x = " << root << std::endl;
    std::cout << std::setprecision(16) << "x0 = " << x0 << std::endl;
    std::cout << "[" << xl << "; " << xr << "]" << std::endl;
    std::cout << "e = " << EPS2 << std::endl << std::endl;

    root = newtonSolve(func, dfunc, EPS3, x0, "1e-9.txt");
    std::cout << "x = " << root << std::endl;
    std::cout << std::setprecision(16) << "x0 = " << x0 << std::endl;
    std::cout << "[" << xl << "; " << xr << "]" << std::endl;
    std::cout << "e = " << EPS3 << std::endl << std::endl;
}

void setInterval(double* xl, double* xr)
{
    double h(2e0);
    for (*xr = A + h; *xr < B; *xr += h)
    {
        if ((func(*xl) * func(*xr) < 0) and (dfunc(*xl) * dfunc(*xr) > 0) and (ddfunc(*xl) * ddfunc(*xr) > 0))
        {
            break;
        }
        else if ((dfunc(*xl) * dfunc(*xr) <= 0))
        {
            //std::cout << *xl << " " << *xr << std::endl;
            *xl = *xr;
        }
    }
}

bool breakCheck1(double xnew, double xold, double epsilon)
{
    std::cout << "x = " << xnew << std::endl;
    if (abs(xnew - xold) <= epsilon)
    {
        std::cout << "|x(n+1) - x(n)| <= e" << std::endl;
        return true;
    }
    std::cout << "|x(n+1) - x(n)| > e" << std::endl;
    return false;
}

bool breakCheck2(function f, double xnew, double epsilon)
{
    if (abs(f(xnew)) <= epsilon)
    {
        std::cout << "|f(x(n+1))| <= e" << std::endl << std::endl;
        return true;
    }
    std::cout << "|f(x(n+1))| > e" << std::endl << std::endl;
    return false;
}

double getStartPosition(function f, function ddf, double xLeft, double xRight)
{
    double x0(xLeft);
    if ((ddf(xLeft) != 0) and (f(xLeft) * ddf(xLeft) > 0)) { x0 = xLeft; }
    else if ((ddf(xRight) != 0) and (f(xRight) * ddf(xRight) > 0)) { x0 = xRight; }
    return x0;
}

double newtonSolve(function f, function df, double epsilon, double x0, std::string filename)
{
    std::ofstream out;
    out.open(filename);

    double root(0), xnew(0), xold(x0);
    int iter(1);
    out << std::setprecision(16) << x0 << std::endl;
    std::cout << "Iteration: " << iter << std::endl;
    xnew = xold - f(xold) / df(xold);
    out << std::setprecision(16) << xnew << std::endl;
    bool a = breakCheck1(xnew, xold, epsilon);
    bool b = breakCheck2(func, xnew, epsilon);
    while (not(a or b))
    {
        iter += 1;
        std::cout << "Iteration: " << iter << std::endl;
        xold = xnew;
        xnew = xold - f(xold) / df(xold);
        out << std::setprecision(16) << xnew << std::endl;
        a = breakCheck1(xnew, xold, epsilon);
        b = breakCheck2(func, xnew, epsilon);
        if (iter > 30)
        {
            std::cout << "Too many" << std::endl;
            break;
        }
    }
    root = xnew;
    out.close();
    return root;
}

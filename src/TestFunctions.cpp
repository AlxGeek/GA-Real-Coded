#include "TestFunctions.hpp"

namespace tf
{

double sphere(std::vector<double> x)
{
    double res = 0;
    for (auto i : x)
        res += pow(i, 2);
    return res;
}

double ellipsoid(std::vector<double> x)
{
    double res = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        res += pow(10, 6 * i / (x.size() - 1.)) * pow(x[i], 2);
    }
    return res;
}

double zakharov(std::vector<double> x)
{
    double s1 = 0;
    double s2 = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        s1 += pow(x[i], 2);
        s2 += 0.5 * (i + 1) * x[i];
    }
    return s1 + pow(s2, 2) + pow(s2, 4);
}

double rosenbrock(std::vector<double> x)
{
    double res = 0;
    for (size_t i = 0; i < x.size() - 1; i++)
    {
        res += pow(1 - x[i], 2) + 100 * pow(x[i + 1] - pow(x[i], 2), 2);
    }
    return res;
}

double ackley(std::vector<double> x)
{
    double s1 = 0;
    double s2 = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        s1 += pow(x[i], 2);
        s2 += cos(2 * M_PI * x[i]);
    }
    return -20 * exp(-0.2 * sqrt(s1 / x.size())) - exp(s2 / x.size()) + 20 + exp(1);
}

double griewangk(std::vector<double> x)
{
    double s1 = 0;
    double m1 = 1;
    for (size_t i = 0; i < x.size(); i++)
    {
        s1 += pow(x[i], 2) / 4000.;
        m1 *= cos(x[i] / sqrt(i + 1));
    }
    return s1 - m1 + 1;
}

double rastrigin(std::vector<double> x)
{
    double s1 = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        s1 += pow(x[i], 2) - 10 * cos(2 * M_PI * x[i]);
    }
    return 10 * x.size() + s1;
}
}

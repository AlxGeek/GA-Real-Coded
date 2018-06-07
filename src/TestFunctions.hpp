#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include <cmath>
#include <vector>
#include <cstddef>

namespace tf
{
double sphere(std::vector<double> x);

double ellipsoid(std::vector<double> x);

double zakharov(std::vector<double> x);

double rosenbrock(std::vector<double> x);

double ackley(std::vector<double> x);

double griewangk(std::vector<double> x);

double rastrigin(std::vector<double> x);
}

#endif // TEST_FUNCTIONS_HPP
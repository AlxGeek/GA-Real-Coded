#ifndef INDIVIDUAL_FUNCTION_HPP
#define INDIVIDUAL_FUNCTION_HPP

#include <vector>
#include <iterator>
#include <random>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Individual.hpp"

class IndividualFunction : public Individual
{
private:
  double (*objetiveFunction)(std::vector<double> x);
  double minDomainValue;
  double maxDomainValue;
  size_t dimensions;

  std::vector<double> genotype;

  std::mt19937 gen;
  std::uniform_real_distribution<> rand0to1;
  std::uniform_real_distribution<> randGenotype;

public:
  IndividualFunction() = default;
  IndividualFunction(double (*objetiveFunction)(std::vector<double>), double minDomainValue, double maxDomainValue, size_t dimensions);

  void initRandom();
  
  void randomMutation(const double probability);
  void nonUniformMutation(const double probability, int t, int gmax, double b);
  void polynomialMutation(const double probability, double eta);

  void setFitness();

  std::vector<double> &getGenotype();

  bool toFile(const char *filename);

  bool toFile(std::ofstream &file) const;

  bool operator<(const IndividualFunction &ind) const;

  size_t getGenotypeLength() const;

  void printGenotype() const;

  double getDistance(IndividualFunction &ind);

  void setDCN(std::vector<IndividualFunction> &survivors);

};

#endif
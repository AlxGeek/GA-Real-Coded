#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <vector>
#include <memory>
#include <cstddef>
#include <random>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <set>
#include <type_traits>

enum class CrossoverOperator
{
  Simple,
  Blx,
  Sbx
};

enum class MutationOperator
{
  Random,
  NonUniform,
  Polynomial
};

enum class ReplacementOperator
{
  Clasic,
  Multi_Dym
};

template <class T>
class GeneticAlgorithm
{
private:
  size_t genotypeLength;
  double mutationProbability;
  double crossoverProbability;
  size_t eliteNumber;

  size_t populationSize;
  std::vector<T> population;
  std::vector<T> offspring;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_int_distribution<> randPopulation;
  std::uniform_int_distribution<> randGenotype;
  std::uniform_real_distribution<> randProb;

  CrossoverOperator crossoverOperator;
  MutationOperator mutationOperator;
  ReplacementOperator replacementOperator;

  void tournament(size_t n);
  void calcFitness();
  void elitism();
  void multiDynamic(double D);

  void randomMutation();
  void nonUniformMutation(int t, int gmax, double b);
  void polynomialMutation(double eta);

  void simpleCross();
  void blxCross(const double alpha);
  void sbxCross(const double u);

  std::vector<size_t> nonDominated(std::vector<T> individuals);

public:
  GeneticAlgorithm(const T &individual, size_t populationSize, double mutationProbability, double crossoverProbability, size_t eliteNumber,
                   CrossoverOperator crossoverOperator, MutationOperator mutationOperator, ReplacementOperator replacementOperator);
  void initPoblation();
  int run(size_t generations);

  const T &getBest();
};

#include "GeneticAlgorithm.tpp"

#endif //GENETIC_ALGORITHM_HPP
#include <iostream>
#include <ctime>
#include <chrono>
#include "string"

#include "GeneticAlgorithm.hpp"
#include "IndividualFunction.hpp"
#include "TestFunctions.hpp"

void printConfig(std::string filename, ReplacementOperator replacementOperator, CrossoverOperator crossoverOperator, MutationOperator mutationOperator)
{
    std::ofstream file(filename);

    file << "Configuracion del Algoritmo Genetico:" << std::endl;
    file << "Reemplazamiento: ";
    switch (replacementOperator)
    {
    case ReplacementOperator::Clasic:
        file << "Reemplazamiento Clasico" << std::endl;
        break;
    case ReplacementOperator::Multi_Dym:
        file << "Multi_Dynamic" << std::endl;
        break;
    default:
        file << "Reemplazamiento por defecto" << std::endl;
        break;
    }
    file << "Cruza: ";
    switch (crossoverOperator)
    {
    case CrossoverOperator::Simple:
        file << "Cruza Simple" << std::endl;
        break;
    case CrossoverOperator::Blx:
        file << "BLX" << std::endl;
        break;
    case CrossoverOperator::Sbx:
        file << "SBX" << std::endl;
        break;
    default:
        file << "Cruza por Defecto" << std::endl;
        break;
    }
    file << "Mutacion: ";
    switch (mutationOperator)
    {
    case MutationOperator::Random:
        file << "Mutacion Aleatoria" << std::endl;
        break;
    case MutationOperator::NonUniform:
        file << "Mutacion No Uniforme" << std::endl;
        break;
    case MutationOperator::Polynomial:
        file << "Mutacion Polinomial" << std::endl;
        break;
    default:
        file << "Mutacion por defecto" << std::endl;
        break;
    }
    file.close();
}

void runTests(std::ostream &output, IndividualFunction individual, int repetitions, size_t popSize, std::string fout, size_t dimension,
              ReplacementOperator replacementOperator, CrossoverOperator crossoverOperator, MutationOperator mutationOperator)
{
    int elite = replacementOperator == ReplacementOperator::Clasic ? 2 : 0;
    double crossoverProb = .8;
    double mutationProb = .01;
    GeneticAlgorithm<IndividualFunction> ga(individual, popSize, mutationProb, crossoverProb, elite, crossoverOperator, mutationOperator, replacementOperator);
    std::ofstream file(fout);
    int sumDuration = 0, durationGa;
    double average = 0;
    double best = 1000000;
    size_t maxItr = 1500;
    for (int i = 0; i < repetitions; i++)
    {
        double c;
        ga.initPoblation();
        output << "Test: " << i << " Initial fitness: " << ga.getBest().getFitness() << " ";
        auto start = std::chrono::steady_clock::now();
        output << "Generations: " << ga.run(maxItr) << " ";
        durationGa = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
        output << "Final fitness: " << (c = ga.getBest().getFitness()) << std::endl;
        file << durationGa << " ";
        ga.getBest().toFile(file);
        average += c;
        sumDuration += durationGa;
        if(c < best){
            best = c;
        }
    }
    output << "Average: " << average / (float)repetitions << " Best fitness: " << best << " Average Optimization duration(ms): " << sumDuration / (float)repetitions << std::endl;
    file.close();
}

void runAllTests(std::ostream &output, int repetitions, size_t popSize, size_t dimension, std::string directory,
                 ReplacementOperator replacementOperator, CrossoverOperator crossoverOperator, MutationOperator mutationOperator)
{
    std::cout << "Sphere" << std::endl;
    IndividualFunction individual(tf::sphere, -600, 600, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "sphere.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);

    std::cout << "ellipsoid" << std::endl;
    individual = IndividualFunction(tf::ellipsoid, -20, 20, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "ellipsoid.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);

    std::cout << "zakharov" << std::endl;
    individual = IndividualFunction(tf::zakharov, -20, 20, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "zakharov.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);

    std::cout << "rosenbrock" << std::endl;
    individual = IndividualFunction(tf::rosenbrock, -20, 20, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "rosenbrock.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);

    std::cout << "ackley" << std::endl;
    individual = IndividualFunction(tf::ackley, -20, 20, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "ackley.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);

    std::cout << "griewangk" << std::endl;
    individual = IndividualFunction(tf::griewangk, -600, 600, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "griewangk.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);

    std::cout << "rastrigin" << std::endl;
    individual = IndividualFunction(tf::rastrigin, -20, 20, dimension);
    runTests(std::cout, individual, repetitions, popSize, directory + "rastrigin.txt", dimension, replacementOperator, crossoverOperator, mutationOperator);
}

int main(int argc, char *argv[])
{

    if (argc != 3)
    {
        std::cout << "Uso: programa pruebas dimension" << std::endl;
        return -1;
    }
    std::srand(unsigned(std::time(0)));
    size_t popSize = 100;

    std::vector<ReplacementOperator> replacementOperators = {ReplacementOperator::Clasic, ReplacementOperator::Multi_Dym};
    std::vector<CrossoverOperator> crossoverOperators = {CrossoverOperator::Simple, CrossoverOperator::Blx, CrossoverOperator::Sbx};
    std::vector<MutationOperator> mutationOperators = {MutationOperator::Random, MutationOperator::NonUniform, MutationOperator::Polynomial};
    std::string directory;

    int configId = 1;
    for (auto ro : replacementOperators)
    {
        for (auto co : crossoverOperators)
        {
            for (auto mo : mutationOperators)
            {
                std::cout << "Configuracion: " << configId << std::endl;
                printConfig("results/c" + std::to_string(configId) + "/config.txt", ro, co, mo);
                directory = "results/c" + std::to_string(configId) + "/";
                runAllTests(std::cout, atoi(argv[1]), popSize, atoi(argv[2]), directory, ro, co, mo);
                configId++;
            }
        }
    }

    return 0;
}
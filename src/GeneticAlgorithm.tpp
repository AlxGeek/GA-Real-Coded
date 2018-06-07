#include "GeneticAlgorithm.hpp"

template <class T>
GeneticAlgorithm<T>::GeneticAlgorithm(
    const T &individual, size_t populationSize, double mutationProbability, double crossoverProbability, size_t eliteNumber,
    CrossoverOperator crossoverOperator, MutationOperator mutationOperator, ReplacementOperator replacementOperator)
    : genotypeLength(individual.getGenotypeLength()), mutationProbability(mutationProbability),
      crossoverProbability(crossoverProbability), eliteNumber(eliteNumber),
      populationSize(populationSize), population(populationSize, individual), offspring(populationSize),
      rd(), gen(rd()), randPopulation(0, populationSize - 1), randGenotype(0, genotypeLength - 1), randProb(0, 1),
      crossoverOperator(crossoverOperator), mutationOperator(mutationOperator), replacementOperator(replacementOperator)
{
}

template <class T>
void GeneticAlgorithm<T>::initPoblation()
{
    for (auto &i : population)
    {
        i.initRandom();
    }
}

template <class T>
void GeneticAlgorithm<T>::tournament(size_t n)
{
    for (size_t i = eliteNumber; i < populationSize; i++)
    {
        size_t selected = randPopulation(gen);
        double min = population[selected].getFitness();
        for (size_t j = 0; j < n - 1; j++)
        {
            size_t k = randPopulation(gen);
            if (population[k].getFitness() < min)
            {
                selected = k;
                min = population[k].getFitness();
            }
        }
        offspring[i] = population[selected];
    }
}

template <class T>
void GeneticAlgorithm<T>::simpleCross()
{
    T aux;
    for (size_t i = eliteNumber; i < populationSize - 1; i += 2)
    {

        if (randProb(gen) < crossoverProbability)
        {
            auto &genotype1 = offspring[i].getGenotype();
            auto &genotype2 = offspring[i + 1].getGenotype();
            aux = offspring[i + 1];
            size_t pos = randGenotype(gen);
            std::copy(std::next(genotype1.begin(), pos), genotype1.end(), std::next(genotype2.begin(), pos));
            std::copy(std::next(aux.getGenotype().begin(), pos), aux.getGenotype().end(), std::next(genotype1.begin(), pos));            
        }
    }
}

template <class T>
void GeneticAlgorithm<T>::blxCross(const double alpha)
{
    for (size_t i = eliteNumber; i < populationSize - 1; i += 2)
    {
        if (randProb(gen) < crossoverProbability)
        {
            auto &genotype1 = offspring[i].getGenotype();
            auto &genotype2 = offspring[i + 1].getGenotype();
            for (size_t i = 0; i < genotypeLength; i++)
            {
                double I = fabs(genotype1[i] - genotype2[i]);
                std::uniform_real_distribution<double> randNumber(std::min(genotype1[i], genotype2[i]) - I * alpha, std::max(genotype1[i], genotype2[i]) + I * alpha);
                genotype1[i] = randNumber(gen);
                genotype2[i] = randNumber(gen);
            }
        }
    }
}

template <class T>
void GeneticAlgorithm<T>::sbxCross(const double eta)
{
    for (size_t i = eliteNumber; i < populationSize - 1; i += 2)
    {
        if (randProb(gen) < crossoverProbability)
        {
            auto &genotype1 = offspring[i].getGenotype();
            auto &genotype2 = offspring[i + 1].getGenotype();
            for (size_t i = 0; i < genotypeLength; i++)
            {
                double u = randProb(gen);
                double betaQ = u <= 0.5 ? pow(2.0 * u, 1.0 / (eta + 1)) : pow(1 / (2 * (1.0 - u)), 1.0 / (eta + 1));
                double gen1 = genotype1[i];
                genotype1[i] = 0.5 * ((1 + betaQ) * gen1 + (1 - betaQ) * genotype2[i]);
                genotype2[i] = 0.5 * ((1 - betaQ) * gen1 + (1 + betaQ) * genotype2[i]);
            }
        }
    }
}

template <class T>
void GeneticAlgorithm<T>::randomMutation()
{
    for (size_t i = eliteNumber; i < populationSize; i++)
    {
        offspring[i].randomMutation(mutationProbability);
    }
}

template <class T>
void GeneticAlgorithm<T>::nonUniformMutation(int t, int gmax, double b)
{
    for (size_t i = eliteNumber; i < populationSize; i++)
    {
        offspring[i].nonUniformMutation(mutationProbability, t, gmax, b);
    }
}

template <class T>
void GeneticAlgorithm<T>::polynomialMutation(double eta)
{
    for (size_t i = eliteNumber; i < populationSize; i++)
    {
        offspring[i].polynomialMutation(mutationProbability, eta);
    }
}

template <class T>
void GeneticAlgorithm<T>::calcFitness()
{
    for (size_t i = eliteNumber; i < populationSize; i++)
    {
        offspring[i].setFitness();
    }
}

template <class T>
void GeneticAlgorithm<T>::elitism()
{
    std::sort(population.begin(), population.end());
    for (size_t i = 0; i < eliteNumber; i++)
    {
        offspring[i] = population[i];
    }
}

template <class T>
void GeneticAlgorithm<T>::multiDynamic(double D)
{
    size_t i;
    size_t c;
    std::vector<size_t> nd;
    std::vector<T> newPop;
    newPop.reserve(populationSize);
    std::vector<T> currentMembers = population;
    currentMembers.reserve(2 * populationSize);
    currentMembers.insert(currentMembers.end(), offspring.begin(), offspring.end());
    typename std::vector<T>::iterator best = std::min_element(currentMembers.begin(), currentMembers.end());
    newPop.push_back(*best);
    currentMembers.erase(best);
    while (newPop.size() < populationSize)
    {
        for (auto &cm : currentMembers)
        {
            cm.setDCN(newPop);
        }
        nd = nonDominated(currentMembers);
        c = 0;
        do
        {
            i = rand() % nd.size();
            c++;
        } while (currentMembers[nd[i]].getDCN() < D && c < nd.size());
        newPop.push_back(currentMembers[nd[i]]);
        currentMembers.erase(currentMembers.begin() + nd[i]);
    }
    population = newPop;
}

template <class T>
std::vector<size_t> GeneticAlgorithm<T>::nonDominated(std::vector<T> individuals)
{
    std::vector<size_t> nd;
    for (size_t i = 0; i < individuals.size(); i++)
    {
        bool nonDominated = true;
        for (size_t j = 0; j < individuals.size(); j++)
        {
            if (individuals[i].getFitness() > individuals[j].getFitness() && individuals[i].getDCN() < individuals[j].getDCN())
            {
                nonDominated = false;
                break;
            }
        }
        if (nonDominated)
        {
            nd.push_back(i);
        }
    }
    return nd;
}

template <class T>
int GeneticAlgorithm<T>::run(size_t generations)
{
    double DI = 10;
    size_t i = 0;
    do
    {
        tournament(2);
        switch (crossoverOperator)
        {
        case CrossoverOperator::Simple:
            simpleCross();
            break;
        case CrossoverOperator::Blx:
            blxCross(0.5);
            break;
        case CrossoverOperator::Sbx:
            sbxCross(1);
            break;
        default:
            simpleCross();
            break;
        }
        switch (mutationOperator)
        {
        case MutationOperator::Random:
            randomMutation();
            break;
        case MutationOperator::NonUniform:
            nonUniformMutation(i, generations, 30);
            break;
        case MutationOperator::Polynomial:
            polynomialMutation(60);
            break;
        default:
            randomMutation();
            break;
        }
        calcFitness();
        switch (replacementOperator)
        {
        case ReplacementOperator::Clasic:
            elitism();
            population = offspring;
            break;
        case ReplacementOperator::Multi_Dym:
            multiDynamic(DI - DI * (double)i / generations);
            break;
        default:
            elitism();
            population = offspring;
            break;
        }

        //std::cout << i << " " << getBest().getFitness() << std::endl;
        //getBest().printGenotype();
        i++;
    } while (i < generations);
    return i;
}

template <class T>
const T &GeneticAlgorithm<T>::getBest()
{
    return *std::min_element(population.begin(), population.end());
}
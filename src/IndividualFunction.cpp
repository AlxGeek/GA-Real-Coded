#include "IndividualFunction.hpp"

IndividualFunction::IndividualFunction(
    double (*objetiveFunction)(std::vector<double>), double minDomainValue, double maxDomainValue, size_t dimensions)
    : objetiveFunction(objetiveFunction), minDomainValue(minDomainValue), maxDomainValue(maxDomainValue),
      dimensions(dimensions), genotype(dimensions),
      rand0to1(0, 1), randGenotype(minDomainValue, maxDomainValue)
{
}

bool IndividualFunction::toFile(const char *filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        return false;
    }
    file << getFitness() << std::endl;
    for (size_t i = 0; i < genotype.size(); i++)
    {
        file << genotype[i] << " ";
    }
    file << std::endl;
    file.close();
    return true;
}

bool IndividualFunction::toFile(std::ofstream &file) const
{
    file << getFitness() << std::endl;
    return true;
}

void IndividualFunction::initRandom()
{
    std::random_device rd1;
    gen.seed(rd1());
    for (size_t i = 0; i < genotype.size(); i++)
    {
        genotype[i] = randGenotype(gen);
    }
    setFitness();
}

void IndividualFunction::randomMutation(const double probability)
{
    for (size_t i = 0; i < genotype.size(); i++)
    {
        if (rand0to1(gen) < probability)
        {
            genotype[i] = randGenotype(gen);
        }
    }
}

void IndividualFunction::nonUniformMutation(const double probability, int t, int gmax, double b)
{
    for (size_t i = 0; i < genotype.size(); i++)
    {
        if (rand0to1(gen) < probability)
        {
            if (rand0to1(gen) <= 0.5)
            {
                genotype[i] = genotype[i] + (maxDomainValue - genotype[i]) * (1 - pow(rand0to1(gen), pow(1 - (double)t / (double)gmax, b)));
            }
            else
            {
                genotype[i] = genotype[i] - (genotype[i] - minDomainValue) * (1 - pow(rand0to1(gen), pow(1 - (double)t / (double)gmax, b)));
            }
        }
    }
}

void IndividualFunction::polynomialMutation(const double probability, double eta)
{
    for (size_t i = 0; i < genotype.size(); i++)
    {
        if (rand0to1(gen) < probability)
        {
            double u = rand0to1(gen);
            if (u <= 0.5)
            {
                genotype[i] = genotype[i] + (pow(2 * u, 1 / (1 + eta)) - 1) * (genotype[i] - minDomainValue);
            }
            else
            {
                genotype[i] = genotype[i] + (1 - pow(2 * (1 - u), 1 / (1 + eta))) * (maxDomainValue - genotype[i]);
            }
        }
    }
}

void IndividualFunction::setFitness()
{
    fitness = objetiveFunction(genotype);
}

std::vector<double> &IndividualFunction::getGenotype()
{
    return genotype;
}

bool IndividualFunction::operator<(const IndividualFunction &ind) const
{
    return fitness < ind.getFitness();
}

size_t IndividualFunction::getGenotypeLength() const
{
    return genotype.size();
}

void IndividualFunction::printGenotype() const
{
    for (auto v : genotype)
    {
        std::cout << v << " ";
    }
    std::cout << std::endl;
}

void IndividualFunction::setDCN(std::vector<IndividualFunction> &survivors)
{
    dcn = getDistance(survivors[0]);
    for (size_t i = 1; i < survivors.size(); i++)
    {
        double temp = getDistance(survivors[i]);
        if (temp < dcn)
        {
            dcn = temp;
        }
    }
}

double IndividualFunction::getDistance(IndividualFunction &ind)
{
    auto genotype2 = ind.getGenotype();
    double distance = 0;
    for (size_t i = 0; i < genotype2.size(); i++)
    {
        distance += pow(genotype[i] - genotype2[i], 2);
    }
    return sqrt(distance);
}

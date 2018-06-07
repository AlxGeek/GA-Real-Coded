#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <cstdlib>
#include <vector>

class Individual
{
private:
protected:
  double fitness;
  double dcn;

public:
  double getFitness() const;

  virtual void initRandom() = 0;

  virtual void setFitness() = 0;

  double getDCN() const;  

  virtual size_t getGenotypeLength() const = 0;

};

#endif // INDIVIDUAL_HPP
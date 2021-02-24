#ifndef INDIVIDUAL_TYPE_H
#define INDIVIDUAL_TYPE_H
#include <iostream>

enum class phenotype
{
  active,
  sporulating,
  spore
};

std::string to_str(phenotype this_ind_type);

#endif // INDIVIDUAL_TYPE_H

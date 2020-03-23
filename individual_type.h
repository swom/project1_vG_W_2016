#ifndef INDIVIDUAL_TYPE_H
#define INDIVIDUAL_TYPE_H
#include <iostream>

enum class individual_type
{
  living,
  sporulating,
  spore
};

std::string to_str(individual_type this_ind_type);

void test_individual_type();
#endif // INDIVIDUAL_TYPE_H

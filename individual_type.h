#ifndef INDIVIDUAL_TYPE_H
#define INDIVIDUAL_TYPE_H
#include <iostream>

enum class ind_type
{
  active,
  sporulating,
  spore
};

std::string to_str(ind_type this_ind_type);

void test_individual_type();
#endif // INDIVIDUAL_TYPE_H

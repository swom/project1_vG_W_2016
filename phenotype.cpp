#include "phenotype.h"
#include <cassert>

std::string to_str(phenotype this_ind_type)
{
  switch (this_ind_type)
    {
    case phenotype::spore:
      return "spore";
    case phenotype::sporulating:
      return "sporulating";
    case phenotype::active:
      return "living";
    }
  return "[Unknown environment_type]";
}

void test_individual_type()
{
#ifndef NDEBUG
  // Conversion to string
  {
    assert(to_str(phenotype::spore) == "spore");
    assert(to_str(phenotype::sporulating) == "sporulating");
    assert(to_str(phenotype::active) == "living");
    assert(to_str(phenotype::active) != "spore");
    assert(to_str(phenotype::spore) != "sporulating");
    assert(to_str(phenotype::sporulating) != "living");
  }
#endif
}

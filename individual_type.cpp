#include "individual_type.h"
#include <cassert>

std::string to_str(individual_type this_ind_type)
{
  switch (this_ind_type)
    {
    case individual_type::spore:
      return "spore";
    case individual_type::sporulating:
      return "sporulating";
    case individual_type::living:
      return "living";
    }
  return "[Unknown environment_type]";
}

void test_individual_type()
{
  // Conversion to string
  {
    assert(to_str(individual_type::spore) == "spore");
    assert(to_str(individual_type::sporulating) == "sporulating");
    assert(to_str(individual_type::living) == "living");
    assert(to_str(individual_type::living) != "spore");
    assert(to_str(individual_type::spore) != "sporulating");
    assert(to_str(individual_type::sporulating) != "living");
  }
}

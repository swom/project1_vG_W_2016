#include "individual_type.h"
#include <cassert>

std::string to_str(ind_type this_ind_type)
{
  switch (this_ind_type)
    {
    case ind_type::spore:
      return "spore";
    case ind_type::sporulating:
      return "sporulating";
    case ind_type::active:
      return "living";
    }
  return "[Unknown environment_type]";
}

void test_individual_type()
{
  // Conversion to string
  {
    assert(to_str(ind_type::spore) == "spore");
    assert(to_str(ind_type::sporulating) == "sporulating");
    assert(to_str(ind_type::active) == "living");
    assert(to_str(ind_type::active) != "spore");
    assert(to_str(ind_type::spore) != "sporulating");
    assert(to_str(ind_type::sporulating) != "living");
  }
}

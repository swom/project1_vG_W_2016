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


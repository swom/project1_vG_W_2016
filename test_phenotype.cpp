#include"tests.h"

void test_phenotype()
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

#include <cassert>
#include "environment.h"
#include "simulation.h"
#include "individual.h"

void test() {
    test_environment();
    test_simulation();
    test_individual();
}

int main() //!OCLINT tests may be long
{
#ifndef NDEBUG
void    test();
#else
    // In release mode, all asserts are removed from the code
    assert(1 == 2);
#endif
    return 0;
}

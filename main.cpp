#include <cassert>
#include <string>
#include "environment.h"
#include "simulation.h"
#include "individual.h"

void test() {
    test_environment();
    test_simulation();
    test_individual();
}


int main(int argc, char ** argv) //!OCLINT tests may be long
{
#ifndef NDEBUG
    test();
#else
    // In release mode, all asserts are removed from the code
    assert(1 == 2);
#endif
    const std::vector<std::string> args(argv, argv + argc);

    //We've already tested, so the program is done
    if (args.size() > 1 && args[1] == "--test") return 0;



    return 0;
}


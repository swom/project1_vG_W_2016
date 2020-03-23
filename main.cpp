#include <cassert>
#include <string>
#include "simulation.h"

void test() {
    test_env_grid_cell();
    test_environment();
    test_individual_type();
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

//    simulation s(19,200);
//    auto time = 10;
//    exec(s, time);
    return 0;
}


#include "tests.h"
#include <cassert>
#include <iostream>
#include <string>



int main(int argc, char ** argv) //!OCLINT tests may be long
{

#ifndef NDEBUG
        test_all();
#endif

    return 0;
}




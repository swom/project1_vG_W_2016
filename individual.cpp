#include "individual.h"
#include <cassert>

individual::individual(double size):
    m_size(size)
{

}

void test_individual(){

    //An individual should be initialized with the defined starting size
    double starting_size = 14.0;
    individual i(starting_size);
    assert(i.get_size() - starting_size < 0.0000001);

}

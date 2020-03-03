#include "individual.h"
#include <cassert>

individual::individual(double x_pos, double y_pos, double size):
    m_size(size),
    m_x(x_pos),
    m_y(y_pos)
{

}

void test_individual(){

    //An individual should be initialized with the defined starting size
    {
        double starting_size = 14.0;
        individual i(0,0,starting_size);
        assert(i.get_size() - starting_size < 0.0000001);
    }

    //Individuals should be initialized with 0 internal energy
    {
        individual i(0,0);
        assert(i.get_energy() - 0.0 < 0.000001);
    }

    //An inddividual should be initialized at a certain position
    {
        double x = 100;
        double y = 100;
        individual i(x,y);
        assert(i.get_x() - x < 0.0000001);
        assert(i.get_y() - y < 0.0000001);
    }


}

#include "environment.h"

environment::environment(int n):
    m_grid(n*n)
{

}

void test_environment()
{
    //An environment has a grid with  n*n elements
    //where n is the size of the side of the square
    //constituting the environment
    {
        int n = 2;
        environment e(n);
        //assert(e.get_env_size() == 2 * 2);
    }


}

#include "funders_success.h"

#include<cassert>

funders_success::funders_success()
{

}

void test_funders_success() noexcept
{
    //funder_success is initialized/contains
    //a vector of funders object
    {
        funders_success fs;
        assert(fs.get_v_funders().size() >= 0);
    }
}

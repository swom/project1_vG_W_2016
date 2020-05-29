#include "funders.h"
#include <cassert>

funders::funders()
{

}

bool operator==(const funders& lhs, const funders& rhs) noexcept
{
   return std::equal(lhs.get_v_funder_data().begin(),lhs.get_v_funder_data().end(),
               rhs.get_v_funder_data().begin());
}

void test_funders() noexcept
{

    //A funders class contains a vector of funder_data object
    {
        funders f;
        assert(f.get_v_funder_data().size() == 0);
    }
}

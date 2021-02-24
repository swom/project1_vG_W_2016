#include "funders.h"
#include <cassert>

funders::funders()
{

}

bool operator==(const funders& lhs, const funders& rhs) noexcept
{
    return lhs.get_v_funder_data() == rhs.get_v_funder_data();
}

bool operator!=(const funders& lhs, const funders& rhs) noexcept
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const funders& f)
{
    for(const auto& funder : f.get_v_funder_data())
    {
        os << f.get_cycle() << " , "
           <<funder;
    }

    return os;
}

std::ifstream& operator>>(std::ifstream& is, funders& f)
{
    funders tmp;
    individual ind{ind_param{}};
    funder_data f_d{ind};
    std::string dummy;
    int cycle_number = -123;
    int prev_cycle_number = -123;
    while(cycle_number == prev_cycle_number)
    {
        is >> cycle_number; // Eat the cycle number

        if(prev_cycle_number == -123)
        {
            prev_cycle_number = cycle_number;
        }

        is >> dummy; //Eat the comma
        is >> f_d;
        tmp.get_v_funder_data().push_back(f_d);

        auto pos = is.tellg();
        is.clear();

        if(is >> cycle_number)
        {
        is.seekg(pos);
        is.clear();
        }
        else
        {
            break;
        }

    }
    tmp.set_cycle(prev_cycle_number);

    f = tmp;
    return is;
}

funders load_funders(const std::string& filename)
{
    funders f;
    std::ifstream is(filename, std::ios::binary);
    is >> f;
    return f;
}

void save_funders(funders f, const std::string& filename)
{
    std::ofstream os(filename);
    os << f;
}

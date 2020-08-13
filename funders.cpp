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

void test_funders() noexcept
{

    //A funders class contains a vector of funder_data object
    {
        funders f;
        assert(f.get_v_funder_data().size() == 0);
    }

    //A funder class is initialized with
    //a member variable meant to indicate to which cycle
    //The funders belong
    //by default = -1;
    {
        funders f;
        assert( f.get_cycle() == -1);
    }

    //Funders can be saved and loaded from a given filenmae
    {
        funders f;
        int number_of_funders = 3;
        for( int i = 0; i != number_of_funders; i++)
        {
            std::vector<int> ancestor_ID(1,i);
            GRN grn{1,2,1};
            funder_data f_d{ancestor_ID,
                        grn};
            f_d.set_success(i);
            f.get_v_funder_data().push_back(f_d);
        }
        std::string filename = "funders.csv";
        save_funders(f, filename);
        funders f1 = load_funders(filename);
        assert( f == f1);
    }
}

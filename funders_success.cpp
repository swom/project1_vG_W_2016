#include "funders_success.h"

#include<cassert>

funders_success::funders_success()
{

}

bool operator==(const funders_success& lhs, const funders_success& rhs) noexcept
{
    return std::equal(lhs.get_v_funders().begin(),lhs.get_v_funders().end(),
                      rhs.get_v_funders().begin());
}

std::ifstream& operator>>(std::ifstream& is, funders_success& f_s)
{
    funders f;
    while(is >> f)
    {
        f_s.get_v_funders().push_back(f);
    }
    f_s.get_v_funders().push_back(f);
    return is;
}

std::ostream& operator<<(std::ostream& os, const funders_success& f_s)
{

    for(const auto& funders : f_s.get_v_funders())
    {
        os << funders;
    }

    return os;
}

funders_success load_funders_success(const std::string& filename)
{
    funders_success f_s;
    std::ifstream is(filename, std::ios::binary);
    is >> f_s;
    return f_s;
}

void save_funders_success(const funders_success& f_s,const std::string& filename)
{
    std::ofstream os(filename);
    os << f_s;
}

void test_funders_success() noexcept
{
    //funder_success is initialized/contains
    //a vector of funders object
    {
        funders_success fs;
        assert(fs.get_v_funders().size() == 0);
    }

    //funder_success can be saved and loaded
    {
        funders_success f_s;
        int number_of_funders_cycles = 3;

        for (int i = 0; i != number_of_funders_cycles; i++)
        {
            funders f;
            int number_of_funders = 3;
            for( int j = 0; j != number_of_funders; j++)
            {
                std::vector<int> ancestor_ID{ i,j};
                GRN grn{1,2,1};
                funder_data f_d{ancestor_ID,
                               grn};
                f.get_v_funder_data().push_back(f_d);
            }
            f.set_cycle(i);
            f_s.get_v_funders().push_back(f);
        }

        std::string filename = "funders_success.csv";
        save_funders_success(f_s, filename);
        funders_success f_s1 = load_funders_success(filename);
        save_funders_success(f_s1, "funders_success1.csv");

        assert(f_s == f_s1);
    }
}

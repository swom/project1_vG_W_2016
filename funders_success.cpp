#include "funders_success.h"
#include<cassert>

funders_success::funders_success()
{

}

bool operator==(const funders_success& lhs, const funders_success& rhs) noexcept
{
    return lhs.get_v_funders() == rhs.get_v_funders();
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

std::string create_funders_success_name(int seed, int change_freq)
{
    return  std::string{
        "funders_success_s" +
        std::to_string(seed) +
                "change_" +
                std::to_string(change_freq) +
                ".csv"
    };
}

GRN find_last_gen_best_ind_grn(const funders_success& funders_success)
{
    auto funders_succ = funders_success.get_v_funders();
    auto last_pop = funders_succ[funders_succ.size() - 2] ;
    auto best_ind = *std::max_element(
                last_pop.get_v_funder_data().begin(),
                last_pop.get_v_funder_data().end(),
                [](const funder_data& lhs, const funder_data& rhs)
    {return lhs.get_success() < rhs.get_success();}
    );
    auto best_grn = best_ind.get_grn();
    return best_grn;
}

funders_success load_funders_success(const std::string& filename)
{
    funders_success f_s;
    std::ifstream is(filename, std::ios::binary);
    if(!is.is_open())
        {
            std::cout << "Could not find specified funders_success*.csv file. \n\n";
            abort();
        }
    is >> f_s;
    return f_s;
}

void save_funders_success(const funders_success& f_s,const std::string& filename)
{
    std::ofstream os(filename);
    os << f_s;
}


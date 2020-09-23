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

    /// The best network of a given generation is found based
    /// on the highest success of the BEFORE-LAST generation
    /// ATTENTION!!! this is done because by design the last
    /// funder object in the funders_success vector does not have
    /// success already calculated
    {
        funders_success funders_success;

        funder_data not_best{std::vector<int>{1},GRN{}};
        funder_data best{std::vector<int>{2},GRN{1,1,1,0.5}};
        not_best.set_success(0);
        best.set_success(10);

        auto n_not_best = 100;
        auto n_best = 2;

        //Create a first funders object that will be the one that will
        //actually be considered for finding the best network
        funders funders_before_last;
        funders_before_last.set_cycle(0);

        for(int i = 0; i != n_not_best; i++ )
            funders_before_last.get_v_funder_data().push_back(not_best);
        for(int i = 0; i != n_best; i++ )
            funders_before_last.get_v_funder_data().push_back(best);
        funders_success.get_v_funders().push_back(funders_before_last);

        //Create a mock last funder object which will not be considered
        funders funders_last;
        funders_last.set_cycle(1);

        for(int i = 0; i != n_not_best; i++ )
            funders_last.get_v_funder_data().push_back(not_best);
        funders_success.get_v_funders().push_back(funders_last);

        //save funders_success since find_best_ind_grn will need to load it
        assert(find_last_gen_best_ind_grn(funders_success) == best.get_grn());
    }
}

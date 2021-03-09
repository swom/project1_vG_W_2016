#include "tests.h"

void test_save_and_load()
{
    //Make a non default funders_success
    funders_success fs;
    funders fu;

    fu.get_v_funder_data().push_back(funder_data{individual{ind_param{12345}}});
    fs.get_v_funders().push_back(fu);

    std::string filename{"funders_success.csv"};
    save_funders_success(fs, filename);

    auto fs2 = load_funders_success(filename);
    assert(fs == fs2);
}

void test_save_and_load_json()
{
    //Make a non default funders_success
    funders_success fs;
    funders fu;

    fu.get_v_funder_data().push_back(funder_data{individual{ind_param{12345}}});
    fs.get_v_funders().push_back(fu);

    std::string filename{"funders_success.json"};
    save_funders_success_json(fs, filename);

    auto f2 = load_funders_success_json(filename);
    assert(fs == f2);
}


void test_funders_success() noexcept
{
    //test csv saving and loading
    test_save_and_load();

    //test json saving and loading
    test_save_and_load_json();

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

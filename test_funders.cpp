#include "tests.h"

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

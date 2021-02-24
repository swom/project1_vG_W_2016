#include "tests.h"


void test_funder_data() noexcept
{

    //A funder_data object can be initialized
    //with a ancestor_ID member
    //and with a GRN member
    {
        std::vector<int> ancestor_ID(2,2);
        GRN grn{1,2,1};
        funder_data f_d{ancestor_ID, grn};
        assert(ancestor_ID == f_d.get_ancestor_ID()
               &&
               grn == f_d.get_grn());
    }

    //A funder_data object can be initialized
    //with an individual reference
    {
        std::vector<int> ancestor_ID(2,2);
        GRN grn{1,2,1};
        individual i{ind_param{}};
        i.get_grn() = grn;
        i.get_ancestor() = ancestor_ID;
        funder_data f_d{i};
        assert(ancestor_ID == f_d.get_ancestor_ID()
               &&
               grn == f_d.get_grn());
    }

    //A funder object is initialized with a success member = 0
    {
        individual i{ind_param{}};
        funder_data f {i};
        assert(f.get_success() == 0);
    }

    //Funder data can be written and read from a streamfile
    {
        std::vector<int> ancestor_ID(2,2);
        GRN grn{1,2,1};
        individual i{ind_param{}};
        i.get_grn() = grn;
        i.get_ancestor() = ancestor_ID;
        funder_data f_d{i};
        const std::string filename = "funder_data.csv";
        save_funder_data(f_d, filename);

        std::ifstream i_s(filename);
        funder_data f_d1{i};
        i_s >> f_d1;
        assert(f_d == f_d1);
    }
}

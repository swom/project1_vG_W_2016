#include "funder_data.h"

#include<cassert>

funder_data::funder_data(std::vector<int> ancestro_ID,
                         GRN grn):
    m_ancestor_ID{ancestro_ID},
    m_grn{grn}
{

}

funder_data::funder_data(const individual& i):
    m_ancestor_ID{i.get_ancestor()},
    m_grn{i.get_grn()}
{

}

bool operator==(const funder_data& lhs, const funder_data& rhs)
{
    return
            lhs.get_grn() == rhs.get_grn()
            && lhs.get_ancestor_ID() == rhs.get_ancestor_ID();
}

bool operator!=(const funder_data& lhs, const funder_data& rhs)
{
    return !(lhs == rhs);
}

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
        individual i;
        i.get_grn() = grn;
        i.get_ancestor() = ancestor_ID;
        funder_data f_d{i};
        assert(ancestor_ID == f_d.get_ancestor_ID()
               &&
               grn == f_d.get_grn());
    }
}

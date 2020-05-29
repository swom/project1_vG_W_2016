#ifndef FUNDER_DATA_H
#define FUNDER_DATA_H

#include "individual.h"

#include <vector>

///Stores essential information about one of the funders of a population cycle
class funder_data
{
public:
    funder_data(std::vector<int> ancestro_ID,
                GRN grn);

    funder_data(const individual& i);

    ///Returns const ref to m_ancestor ID
    const std::vector<int>& get_ancestor_ID() const noexcept {return m_ancestor_ID;}

    ///Returns const ref to m_grn
    const GRN& get_grn() const noexcept {return m_grn;}

private:
    ///The ancestor_ID of the funder individual
    std::vector<int> m_ancestor_ID;

    ///The GRN of the funder individual
    GRN m_grn;
};

///Compares if two funders data have the same GRN and ancestor_ID
bool operator==(const funder_data& lhs, const funder_data& rhs);

///Checks if two funders data have different GRN and ancestor_ID
bool operator!=(const funder_data& lhs, const funder_data& rhs);

void test_funder_data() noexcept;

#endif // FUNDER_DATA_H

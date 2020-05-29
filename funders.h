#ifndef FUNDERS_H
#define FUNDERS_H

#include"funder_data.h"

#include<vector>

class funders
{
public:
    funders();

    ///Returns const ref to m_v_funder_data
    const std::vector<funder_data>& get_v_funder_data() const noexcept {return m_v_funder_data;}

    ///Returns ref to m_v_funder_data
    std::vector<funder_data>& get_v_funder_data() noexcept {return m_v_funder_data;}

private:

    ///The vector containing the data related to each funder of the population
    std::vector<funder_data> m_v_funder_data;
};

void test_funders() noexcept;
#endif // FUNDERS_H

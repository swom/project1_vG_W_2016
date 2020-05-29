#ifndef FUNDER_SUCCESS_H
#define FUNDER_SUCCESS_H

#include "funders.h"

#include<vector>

class funders_success
{
public:
    funders_success();

    ///Returns a const ref to the vector of founders
    /// of different cycles
    const std::vector<funders>& get_v_funders() const noexcept {return m_v_funders;}

    ///Returns a ref to the vector of founders
    /// of different cycles
    std::vector<funders>& get_v_funders() noexcept {return m_v_funders;}

private:
    ///Te vector containing the funders of different cycles
    std::vector<funders> m_v_funders;
};

void test_funders_success() noexcept;

#endif // FUNDER_SUCCESS_H

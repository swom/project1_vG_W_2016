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
    ///of different cycles
    std::vector<funders>& get_v_funders() noexcept {return m_v_funders;}

private:
    ///Te vector containing the funders of different cycles
    std::vector<funders> m_v_funders;
};

///Compares to objects of type funders_success to see that all their elements are the same
bool operator==(const funders_success& lhs, const funders_success& rhs) noexcept;

///Creates the name for the file of funders_success of an evolutionary run
/// given the seed and the frequency of environmentla change
std::string create_funder_success_name(int seed, int change_freq);

///Retruns the GRN of the individual with the highest success
/// of the funder_success
/// ATTENTION!!! by design funders_scuccess objects
/// will have as the last element a population without
/// a registered success, therefore to find the best network
/// we will look at the before-last element of the funders vector
GRN find_best_ind_grn(const funders_success& funders_success);

///Loads a funders?success object from a given filename
 funders_success load_funders_success(const std::string& filename);

///Saves a funders success object to a given filename
void save_funders_success(const funders_success& f_s,const std::string& filename);

void test_funders_success() noexcept;

#endif // FUNDER_SUCCESS_H

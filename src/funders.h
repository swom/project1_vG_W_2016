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

    ///Returns the cycle number to which funders belong
    int get_cycle() const noexcept { return  m_cycle;}

    ///Sets the value of m_cycle to indicate from which cycle the funders are from
    void set_cycle(int cycle) noexcept { m_cycle = cycle;}

private:

    ///The vector containing the data related to each funder of the population
    std::vector<funder_data> m_v_funder_data;

    ///This variable is used during saving to external file to indicate to
    /// which cycle the funders belong
    /// Its value will be changed/assigned at the moment of saving
    /// based on the position of the funder object in the funder
    /// vecto of funders_success
    int m_cycle = -1;
};

///Checks two vectors of funders to see if they are the same
///For ancestor_IDs, GRN and success
bool operator==(const funders& lhs, const funders& rhs) noexcept;


///Checks two vectors of funders to see if they are the NOT the same
///For ancestor_IDs, GRN and success
bool operator!=(const funders& lhs, const funders& rhs) noexcept;

///Prints a funders object to an ofstream
std::ostream& operator<<(std::ostream& os, const funders& f);

///Instantiates a funders object from a ifstream
std::ifstream& operator>>(std::ifstream& is, funders& f);

///Load funders from a given file_name
funders load_funders(const std::string& filename);

///Saves funders to a given filename
void save_funders(funders f, const std::string& filename);

void test_funders() noexcept;
#endif // FUNDERS_H

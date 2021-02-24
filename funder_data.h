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

    ///Returns the success of a funder
    double get_success() const noexcept {return m_success;}

    ///Sets the success of an individual
    void set_success(double success) noexcept {m_success = success;}

private:
    ///The ancestor_ID of the funder individual
    std::vector<int> m_ancestor_ID;

    ///The GRN of the funder individual
    GRN m_grn;

    ///The success fo an individual represented by the fraction of descendants
    /// composing the final population
    double m_success = 0;
};

///Compares if two funders data have the same GRN and ancestor_ID
bool operator==(const funder_data& lhs, const funder_data& rhs);

///Checks if two funders data have different GRN and ancestor_ID
bool operator!=(const funder_data& lhs, const funder_data& rhs);

///Instantiates a funder_data object from an ifstream
std::ifstream& operator>>(std::ifstream& is, funder_data& f_d);

///Prints funder_data object to an ofstream
std::ostream& operator<<(std::ostream& os, const funder_data& f_d);

///Covert the string representing the ancestor ID in the ifstream
/// into the ancestor ID int vector
std::vector<int> load_ancestor_ID(std::ifstream& is);

///Saves the funder_data to a given text_file
void save_funder_data(const funder_data& f_d, const std::string& filename);


#endif // FUNDER_DATA_H

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

            bool grn = lhs.get_grn() == rhs.get_grn();
            bool ancestor = lhs.get_ancestor_ID() == rhs.get_ancestor_ID();
            bool success = lhs.get_success() - rhs.get_success() < 0.001
            && lhs.get_success() - rhs.get_success() > -0.001;

            return grn && ancestor && success;
}

bool operator!=(const funder_data& lhs, const funder_data& rhs)
{
    return !(lhs == rhs);
}


std::ostream& operator<<(std::ostream& os, const funder_data& f_d)
{

    for(auto ancestor : f_d.get_ancestor_ID())
    {
        os << ancestor << " ";
    }
    os << "# , ";

    os << f_d.get_grn();

    os << f_d.get_success() << std::endl;

    return os;
}

std::ifstream& operator>>(std::ifstream& is, funder_data& f_d)
{

    double success;
    std::string dummy; // To remove the annotation in the file


    auto ancestor_ID = load_ancestor_ID(is);

    auto grn = f_d.get_grn();
    is >> grn;
    is >> success;
    is.get();//Eat the endl

            f_d = funder_data{
            ancestor_ID,
            grn};
            f_d.set_success(success);

    return is;
}

std::vector<int> load_ancestor_ID(std::ifstream& is)
{
    std::vector<int> ancestor_ID;
    std::string ancestor;
    std::getline(is, ancestor, '#');
    std::stringstream iss( ancestor );
    int number;
    while ( iss >> number )
        ancestor_ID.push_back( number );

    std::string dummy;
    is >> dummy;//Eat the comma after the hashtag
    return ancestor_ID;
}

void save_funder_data(const funder_data& f_d,
                      const std::string& filename)
{
    std::ofstream o_s(filename);
    o_s << f_d;
}

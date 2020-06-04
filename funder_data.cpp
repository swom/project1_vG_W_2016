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
            && lhs.get_ancestor_ID() == rhs.get_ancestor_ID()
            && lhs.get_success() == rhs.get_success();
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

    //A funder object is initialized with a success member = 0
    {
        individual i;
        funder_data f {i};
        assert(f.get_success() == 0);
    }

    //Funder data can be written and read from a streamfile
    {
        std::vector<int> ancestor_ID(2,2);
        GRN grn{1,2,1};
        individual i;
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

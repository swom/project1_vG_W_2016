#include "ind_param.h"
#include <cassert>

ind_param::ind_param(double radius,
                     double treshold_energy,
                     double uptake_rate,
                     double metabolic_rate,
                     double spor_metabolic_rate,
                     int transformation_time,
                     double wiggle_room,
                     double metab_secretion_rate):
    m_metabolic_rate{metabolic_rate},
    m_metab_secr_rate{metab_secretion_rate},
    m_radius{radius},
    m_spor_metabolic_rate{spor_metabolic_rate},
    m_transformation_time{transformation_time},
    m_treshold_energy{treshold_energy},
    m_uptake_rate{uptake_rate},
    m_wiggle_room{wiggle_room}
{
    assert(m_metabolic_rate > -0.000001);
    assert(m_metab_secr_rate > -0.000001);
    assert(m_radius > -0.000001);
    assert(m_transformation_time > 0);
    assert(m_treshold_energy > -0.000001);
    assert(m_uptake_rate > -0.000001);
    assert(m_wiggle_room > -0.000001);
}

std::ostream& operator<<(std::ostream& os, const ind_param& p)
{
    os << p.get_radius()
       << " , " << p.get_uptake_rate()
       << " , " << p.get_wiggle_room()
       << " , " << p.get_metabolic_rate()
       << " , " << p.get_metab_secr_rate()
       << " , " << p.get_treshold_energy()
       << " , " << p.get_spor_metabolic_rate()
       << " , " << p.get_transformation_time();
    return os;
}

std::ifstream& operator>>(std::ifstream& is, ind_param& p)
{

    double radius;
    double treshold_energy;
    double uptake_rate;
    double metabolic_rate;
    double spor_metabolic_rate;
    int transformation_time;
    double wiggle_room;
    double metab_secretion_rate;
    std::string dummy; // To remove the annotation in the file
    is
            >> radius >> dummy
            >> uptake_rate >> dummy
            >> wiggle_room >> dummy
            >> metabolic_rate >> dummy
            >> metab_secretion_rate >> dummy
            >> treshold_energy >> dummy
            >> spor_metabolic_rate >> dummy
            >> transformation_time;

    p = ind_param {radius,
            treshold_energy,
            uptake_rate,
            metabolic_rate,
            spor_metabolic_rate,
            transformation_time,
            wiggle_room,
            metab_secretion_rate
};

    return is;
}

bool operator==(const ind_param& lhs, const ind_param& rhs) noexcept
{
    return
            lhs.get_radius() == rhs.get_radius()
            && lhs.get_uptake_rate() == rhs.get_uptake_rate()
            && lhs.get_wiggle_room() == rhs.get_wiggle_room()
            && lhs.get_metabolic_rate() == rhs.get_metabolic_rate()
            && lhs.get_metab_secr_rate() == rhs.get_metab_secr_rate()
            && lhs.get_treshold_energy() == rhs.get_treshold_energy()
            && lhs.get_spor_metabolic_rate() == rhs.get_spor_metabolic_rate()
            && lhs.get_transformation_time() == rhs.get_transformation_time()
            ;
}

void save_ind_parameters(
        const ind_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p.get_radius()  << " , " <<
         p.get_uptake_rate() << " , " <<
         p.get_wiggle_room() << " , " <<
         p.get_metabolic_rate() << " , " <<
         p.get_metab_secr_rate() << " , " <<
         p.get_treshold_energy() << " , " <<
         p.get_spor_metabolic_rate() << " , " <<
         p.get_transformation_time();
}

ind_param load_ind_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    double radius;
    double treshold_energy;
    double uptake_rate;
    double metabolic_rate;
    double spor_metabolic_rate;
    int transformation_time;
    double wiggle_room;
    double metab_secretion_rate;
    std::string dummy; // To remove the annotation in the file

    f
            >> radius >> dummy
            >> uptake_rate >> dummy
            >> wiggle_room >> dummy
            >> metabolic_rate >> dummy
            >> metab_secretion_rate >> dummy
            >> treshold_energy >> dummy
            >> spor_metabolic_rate >> dummy
            >> transformation_time;

    ind_param p{
        radius,
                treshold_energy,
                uptake_rate,
                metabolic_rate,
                spor_metabolic_rate,
                transformation_time,
                wiggle_room,
                metab_secretion_rate

    };
    return p;
}

void test_ind_param() noexcept  //!OCLINT
{
    //It is possible to save and load ind parameters from a file
    {
        double radius  = 0.9;
        double treshold_energy = 15;
        double uptake_rate = 0.12;
        double metabolic_rate = 0.031;
        double spor_metabolic_rate = 0.8;
        int transformation_time = 7;
        double wiggle_room = 0.021;
        double metab_secretion_rate = 2;
        ind_param p{
            radius,
                    treshold_energy,
                    uptake_rate,
                    metabolic_rate,
                    spor_metabolic_rate,
                    transformation_time,
                    wiggle_room,
                    metab_secretion_rate
        };
        const std::string filename = "ind_param.csv";
        save_ind_parameters(p, filename);
        const ind_param q = load_ind_parameters(filename);
        assert(p == q);
        //Test >> operator overload
        std::ifstream f(filename);
        ind_param s;
        f >> s;
        assert(s == p);
    }
    {
        const ind_param p;
        std::ostringstream s;
        s << p;
    }
}

#include "ind_param.h"
#include <cassert>

ind_param::ind_param(double radius,
                     double treshold_energy,
                     double uptake_rate,
                     double uptake_range,
                     double metabolic_rate,
                     double reproduction_prob,
                     double reproduction_range,
                     double spor_metabolic_rate,
                     double spor_metab_range,
                     int transformation_time,
                     int transformation_range,
                     double metab_secretion_rate,
                     double range_on_change_ratio):
    m_metabolic_rate{metabolic_rate},
    m_metab_secr_rate{metab_secretion_rate},
    m_radius{radius},
    m_range_on_change_ratio{range_on_change_ratio},
    m_repr_prob{reproduction_prob},
    m_repr_prob_range{m_repr_prob - reproduction_range, m_repr_prob + reproduction_range},
    m_spor_metabolic_rate{spor_metabolic_rate},
    m_spor_metabolic_range{m_spor_metabolic_rate - spor_metab_range, m_spor_metabolic_rate + spor_metab_range},
    m_transformation_time{transformation_time},
    m_transformation_range{m_transformation_time - transformation_range, m_transformation_time + transformation_range},
    m_treshold_energy{treshold_energy},
    m_uptake_rate{uptake_rate},
    m_uptake_range{m_uptake_rate - uptake_range, m_uptake_rate + uptake_range}
{
    assert(m_metabolic_rate > -0.000001);
    assert(m_metab_secr_rate > -0.000001);
    assert(m_radius > -0.000001);
    assert(reproduction_prob < 1.0000001 && reproduction_prob > -0.0000000001 );
    assert(m_transformation_time > 0);
    assert(m_treshold_energy > -0.000001);
    assert(m_uptake_rate > -0.000001);
    assert(m_range_on_change_ratio > -0.0000001);
    assert(m_repr_prob_range.min() > -0.0000001);
    assert(m_spor_metabolic_range.min() > -0.0000001);
    assert(m_transformation_range.min() > -0.0000001);
    assert(m_uptake_range.min() > -0.0000001);
}

std::ostream& operator<<(std::ostream& os, const ind_param& p)
{
    os << p.get_radius()
       << " , " << p.get_uptake_rate()
       << " , " << (p.get_uptake_range().max() - p.get_uptake_range().min()) / 2
       << " , " << p.get_metabolic_rate()
       << " , " << p.get_metab_secr_rate()
       << " , " << p.get_repr_prob()
       << " , " << (p.get_repr_range().max() - p.get_repr_range().min()) / 2
       << " , " << p.get_treshold_energy()
       << " , " << p.get_spor_metabolic_rate()
       << " , " << (p.get_spor_metabolic_range().max() - p.get_spor_metabolic_range().min()) / 2
       << " , " << p.get_transformation_time()
       << " , " << (p.get_transformation_range().max() - p .get_transformation_range().min()) / 2
       << " , " << p.get_range_on_change_ratio()
          ;
    return os;
}

std::ifstream& operator>>(std::ifstream& is, ind_param& p)
{

    double radius;
    double treshold_energy;
    double uptake_rate;
    double uptake_range;
    double metabolic_rate;
    double reproduction_probability;
    double reproduction_range;
    double spor_metabolic_rate;
    double spor_metabolic_range;
    int transformation_time;
    int transformation_range;
    double metab_secretion_rate;
    double range_on_change_ratio;
    std::string dummy; // To remove the annotation in the file
    is
            >> radius >> dummy
            >> uptake_rate >> dummy
            >> uptake_range >> dummy
            >> metabolic_rate >> dummy
            >> metab_secretion_rate >> dummy
            >> reproduction_probability >> dummy
            >> reproduction_range >> dummy
            >> treshold_energy >> dummy
            >> spor_metabolic_rate >> dummy
            >> spor_metabolic_range >> dummy
            >> transformation_time >> dummy
            >> transformation_range >> dummy
            >> range_on_change_ratio;

    p = ind_param {radius,
            treshold_energy,
            uptake_rate,
            uptake_range,
            metabolic_rate,
            reproduction_probability,
            reproduction_range,
            spor_metabolic_rate,
            spor_metabolic_range,
            transformation_time,
            transformation_range,
            metab_secretion_rate,
            range_on_change_ratio
};

    return is;
}

bool operator==(const ind_param& lhs, const ind_param& rhs) noexcept
{
    return
            lhs.get_radius() == rhs.get_radius()
            && lhs.get_treshold_energy() == rhs.get_treshold_energy()
            && lhs.get_uptake_rate() == rhs.get_uptake_rate()
            && lhs.get_uptake_range() == rhs.get_uptake_range()
            && lhs.get_metabolic_rate() == rhs.get_metabolic_rate()
            && lhs.get_repr_prob() == rhs.get_repr_prob()
            && lhs.get_repr_range() == rhs.get_repr_range()
            && lhs.get_metab_secr_rate() == rhs.get_metab_secr_rate()
            && lhs.get_spor_metabolic_rate() == rhs.get_spor_metabolic_rate()
            && lhs.get_spor_metabolic_range() == rhs.get_spor_metabolic_range()
            && lhs.get_transformation_time() == rhs.get_transformation_time()
            && lhs.get_transformation_range() == rhs.get_transformation_range()
            && lhs.get_range_on_change_ratio() == rhs.get_range_on_change_ratio()
            ;
}

bool operator!=(const ind_param& lhs, const ind_param& rhs) noexcept
{
    return !(lhs==rhs);
}

double change_ind_uptake(ind_param& i, std::minstd_rand& rng)
{
    auto old_uptake = i.get_uptake_rate();
    auto new_uptake = i.get_uptake_rate();
    auto uptake_step =
            (i.get_uptake_range().max() - i.get_uptake_range().min())
            / i.get_range_on_change_ratio();

    while ( true ) {
        new_uptake = i.get_uptake_range()(rng);
        if(new_uptake > old_uptake + uptake_step
                || new_uptake < old_uptake - uptake_step)
        {
            return new_uptake;
        }
    }
}

double change_ind_repr_prob(ind_param& i, std::minstd_rand& rng)
{
    auto old_repr_prob = i.get_repr_prob();
    auto new_repr_prob = i.get_repr_prob();
    auto repr_step =
            (i.get_repr_range().max() - i.get_repr_range().min())
            / i.get_range_on_change_ratio();

    while ( true ) {
        new_repr_prob = i.get_repr_range()(rng);
        if(new_repr_prob > old_repr_prob + repr_step
                || new_repr_prob < old_repr_prob - repr_step )
        {
            return new_repr_prob;
        }
    }
}

double change_ind_spor_metab(ind_param& i, std::minstd_rand& rng)
{
    auto old_spor_metab = i.get_spor_metabolic_rate();
    auto new_spor_metab = i.get_spor_metabolic_rate();
    auto spor_metab_step =
            (i.get_spor_metabolic_range().max() - i.get_spor_metabolic_range().min())
            / i.get_range_on_change_ratio();

    while ( true ) {
        new_spor_metab = i.get_spor_metabolic_range()(rng);
        if(new_spor_metab > old_spor_metab + spor_metab_step
                ||new_spor_metab < old_spor_metab - spor_metab_step )
            return new_spor_metab;
    }
}



ind_param change_ind_param_unif( ind_param i,  std::minstd_rand& rng)
{

    i.set_uptake_rate(change_ind_uptake(i, rng));

    i.set_spor_metabolic_rate(change_ind_spor_metab(i, rng));

    i.set_repr_prob(change_ind_repr_prob(i,rng));

    i.set_transformation_time(i.get_transformation_range()(rng));

    return i;
}

void save_ind_parameters(
        const ind_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p;
}

ind_param load_ind_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    ind_param p;
    f >> p;
    return p;
}

void test_ind_param() noexcept  //!OCLINT
{
    //It is possible to save and load ind parameters from a file
    {
        double radius  = 0.9;
        double treshold_energy = 15;
        double uptake_rate = 0.12;
        double uptake_range = uptake_rate / 2;
        double metabolic_rate = 0.031;
        double reproduction_prob = 0.123456;
        double reproduction_range = reproduction_prob / 2;
        double spor_metabolic_rate = 0.8;
        double spor_metabolic_range  = reproduction_prob / 2;
        int transformation_time = 7;
        int transformation_range = transformation_time / 2;
        double metab_secretion_rate = 2;
        double range_on_change_ratio = 10;
        ind_param p{
            radius,
                    treshold_energy,
                    uptake_rate,
                    uptake_range,
                    metabolic_rate,
                    spor_metabolic_rate,
                    spor_metabolic_range,
                    reproduction_prob,
                    reproduction_range,
                    transformation_time,
                    transformation_range,
                    metab_secretion_rate,
                    range_on_change_ratio
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

    //Ind parameters can be changed based one the
    //range related to the parameters
    //and the minimum fraction of the range for which they change
    {
        std::minstd_rand rng;
        ind_param i{};
        auto new_i = change_ind_param_unif(i,rng);
        assert( i != new_i);
    }

    //Parameters will always change by a minimum amount dictated by the range_on_change_ratio
    // e.g. : if  range_on_change_ratio = 10, then
    // a parameter with range(max param possible value - min param possible value) = 40
    // will always change with a minimum step of 40/10 = 4
    {
        std::minstd_rand rng;
        double radius  = 0.9;
        double treshold_energy = 15;
        double uptake_rate = 0.12;
        double uptake_range = uptake_rate / 2;
        double metabolic_rate = 0.031;
        double reproduction_prob = 0.123456;
        double reproduction_range = reproduction_prob / 2;
        double spor_metabolic_rate = 0.8;
        double spor_metabolic_range  = reproduction_prob / 2;
        int transformation_time = 7;
        int transformation_range = transformation_time / 2;
        double metab_secretion_rate = 2;
        double range_on_change_ratio = 10;
        ind_param i_p {radius,
                    treshold_energy,
                    uptake_rate,
                    uptake_range,
                    metabolic_rate,
                    reproduction_prob,
                    reproduction_range,
                    spor_metabolic_rate,
                    spor_metabolic_range,
                    transformation_time,
                    transformation_range,
                    metab_secretion_rate,
                    range_on_change_ratio
                      };
        auto uptake_step = (i_p.get_uptake_range().max() - i_p.get_uptake_range().min()) / i_p.get_range_on_change_ratio();
        auto spor_metab_step =
                (i_p.get_spor_metabolic_range().max() - i_p.get_spor_metabolic_range().min()) / i_p.get_range_on_change_ratio();
        auto repr_step = (i_p.get_repr_range().max() - i_p.get_repr_range().min()) / i_p.get_range_on_change_ratio();
        auto new_i = i_p;
        int repeats = 100;
        for( int i = 0; i != repeats; i++)
        {
            new_i = change_ind_param_unif(i_p,rng);
            auto uptake_delta = new_i.get_uptake_rate() - i_p.get_uptake_rate();
            assert(uptake_delta * uptake_delta >= uptake_step  * uptake_step);
            auto repr_prob_delta = new_i.get_repr_prob() - i_p.get_repr_prob();
            assert( repr_prob_delta * repr_prob_delta >= repr_step * repr_step);
            auto spor_metabolic_rate_delta =
                    new_i.get_spor_metabolic_rate() - i_p.get_spor_metabolic_rate();
            assert( spor_metabolic_rate_delta * spor_metabolic_rate_delta >=
                    spor_metab_step * spor_metab_step );
        }
    }
}

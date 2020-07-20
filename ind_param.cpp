#include "ind_param.h"
#include <cassert>

ind_param::ind_param(double radius,
                     double treshold_energy,
                     double uptake_rate,
                     double uptake_rate_mean,
                     double uptake_rate_var,
                     double metabolic_rate,
                     double reproduction_prob,
                     double reproduction_prob_mean,
                     double reproduction_prob_var,
                     double spor_metabolic_rate,
                     double spor_metabolic_rate_mean,
                     double spor_metabolic_rate_var,
                     int transformation_time,
                     int transformation_time_mean,
                     int transformation_range,
                     double metab_secretion_rate):
    m_metabolic_rate{metabolic_rate},
    m_metab_secr_rate{metab_secretion_rate},
    m_radius{radius},
    m_repr_prob{reproduction_prob},
    m_mean_repr_prob{reproduction_prob_mean},
    m_var_repr_prob{reproduction_prob_var},
    m_spor_metabolic_rate{spor_metabolic_rate},
    m_mean_spor_metabolic_rate{spor_metabolic_rate_mean},
    m_var_spor_metabolic_rate{spor_metabolic_rate_var},
    m_transformation_time{transformation_time},
    m_mean_transformation_time{transformation_time_mean},
    m_transformation_range{transformation_range},
    m_treshold_energy{treshold_energy},
    m_uptake_rate{uptake_rate},
    m_mean_uptake_rate{uptake_rate_mean},
    m_var_uptake_rate{uptake_rate_var}
  //!OCLINT
{
    assert(m_metabolic_rate > -0.000001);
    assert(m_metab_secr_rate > -0.000001);
    assert(m_radius > -0.000001);
    assert(reproduction_prob < 1.0000001
           && reproduction_prob > -0.0000000001 );
    assert(m_mean_repr_prob > -0.00000001 && m_mean_repr_prob < 1.00001);
    assert(m_var_repr_prob > -0.00000001 && m_var_repr_prob < 1.00001);
    assert(m_var_repr_prob < m_mean_repr_prob / 3.0);
    assert(m_spor_metabolic_rate > -0.0001);
    assert(m_mean_spor_metabolic_rate > -0.000001);
    assert(m_var_spor_metabolic_rate > -0.000001);
    assert(m_var_spor_metabolic_rate < m_mean_spor_metabolic_rate / 3.0);
    assert(m_transformation_time > 0);
    assert(m_mean_transformation_time - m_transformation_range > -0.0000001);
    assert(m_treshold_energy > -0.000001);
    assert(m_uptake_rate > -0.000001);
    assert(m_mean_uptake_rate > -0.00001);
    assert(m_var_uptake_rate > -0.0001);
    assert(m_var_uptake_rate < m_mean_uptake_rate / 3.0);
}

std::ostream& operator<<(std::ostream& os, const ind_param& p)
{
    os << p.get_radius()
       << " , " << p.get_uptake_rate()
       << " , " << p.get_uptake_mean()
       << " , " << p.get_uptake_var()
       << " , " << p.get_metabolic_rate()
       << " , " << p.get_metab_secr_rate()
       << " , " << p.get_repr_prob()
       << " , " << p.get_repr_prob_mean()
       << " , " << p.get_repr_prob_var()
       << " , " << p.get_treshold_energy()
       << " , " << p.get_spor_metabolic_rate()
       << " , " << p.get_spor_metabolic_rate_mean()
       << " , " << p.get_spor_metabolic_rate_var()
       << " , " << p.get_transformation_time()
       << " , " << p.get_transformation_time_mean()
       << " , " << p.get_transformation_range()
       << " , " << p.get_metab_secr_rate()
          ;
    return os;
}

std::ifstream& operator>>(std::ifstream& is, ind_param& p)
{

    double radius;
    double treshold_energy;
    double uptake_rate;
    double uptake_rate_mean;
    double uptake_rate_var;
    double metabolic_rate;
    double reproduction_prob;
    double reproduction_prob_mean;
    double reproduction_prob_var;
    double spor_metabolic_rate;
    double spor_metabolic_rate_mean;
    double spor_metabolic_rate_var;
    int transformation_time;
    int transformation_time_mean;
    int transformation_range;
    double metab_secretion_rate;
    std::string dummy; // To remove the annotation in the file
    is
            >> radius >> dummy
            >> uptake_rate >> dummy
            >> uptake_rate_mean >> dummy
            >> uptake_rate_var >> dummy
            >> metabolic_rate >> dummy
            >> metab_secretion_rate >> dummy
            >> reproduction_prob >> dummy
            >> reproduction_prob_mean >> dummy
            >> reproduction_prob_var >> dummy
            >> treshold_energy >> dummy
            >> spor_metabolic_rate >> dummy
            >> spor_metabolic_rate_mean >> dummy
            >> spor_metabolic_rate_var >> dummy
            >> transformation_time >> dummy
            >> transformation_time_mean >> dummy
            >> transformation_range >> dummy
            >> metab_secretion_rate
            ;

    p = ind_param {radius,
            treshold_energy,
            uptake_rate,
            uptake_rate_mean,
            uptake_rate_var,
            metabolic_rate,
            reproduction_prob,
            reproduction_prob_mean,
            reproduction_prob_var,
            spor_metabolic_rate,
            spor_metabolic_rate_mean,
            spor_metabolic_rate_var,
            transformation_time,
            transformation_time_mean,
            transformation_range,
            metab_secretion_rate,
};

    return is;
}

bool operator==(const ind_param& lhs, const ind_param& rhs) noexcept
{
    return lhs.get_radius() - rhs.get_radius() < 0.00001
            && lhs.get_radius() - rhs.get_radius() > -0.00001
            && lhs.get_treshold_energy() - rhs.get_treshold_energy() < 0.00001
            && lhs.get_treshold_energy() - rhs.get_treshold_energy() > -0.00001
            && lhs.get_uptake_rate() - rhs.get_uptake_rate() < 0.0001
            && lhs.get_uptake_rate() - rhs.get_uptake_rate() > -0.0001
            && lhs.get_uptake_mean() - rhs.get_uptake_mean() < 0.0001
            && lhs.get_uptake_mean() - rhs.get_uptake_mean() > -0.0001
            && lhs.get_uptake_var() - rhs.get_uptake_var() < 0.00001
            && lhs.get_uptake_var() - rhs.get_uptake_var() > -0.00001
            && lhs.get_repr_prob() -  rhs.get_repr_prob() < 0.00001
            && lhs.get_repr_prob() -  rhs.get_repr_prob() > -0.00001
            && lhs.get_repr_prob_mean() -  rhs.get_repr_prob_mean() < 0.00001
            && lhs.get_repr_prob_mean() -  rhs.get_repr_prob_mean() > -0.00001
            && lhs.get_repr_prob_var() -  rhs.get_repr_prob_var() < 0.00001
            && lhs.get_repr_prob_var() -  rhs.get_repr_prob_var() > -0.00001
            && lhs.get_metabolic_rate() - rhs.get_metabolic_rate() < 0.00001
            && lhs.get_metabolic_rate() - rhs.get_metabolic_rate() > -0.00001
            && lhs.get_metab_secr_rate() - rhs.get_metab_secr_rate() < 0.00001
            && lhs.get_metab_secr_rate() - rhs.get_metab_secr_rate() > -0.00001
            && lhs.get_spor_metabolic_rate() - rhs.get_spor_metabolic_rate() < 0.00001
            && lhs.get_spor_metabolic_rate() - rhs.get_spor_metabolic_rate() > -0.00001
            && lhs.get_spor_metabolic_rate_mean() - rhs.get_spor_metabolic_rate_mean() < 0.0001
            && lhs.get_spor_metabolic_rate_mean() - rhs.get_spor_metabolic_rate_mean() > -0.0001
            && lhs.get_spor_metabolic_rate_var() - rhs.get_spor_metabolic_rate_var() < 0.00001
            && lhs.get_spor_metabolic_rate_var() - rhs.get_spor_metabolic_rate_var() >  -0.00001
            && lhs.get_transformation_time() == rhs.get_transformation_time()
            && lhs.get_transformation_time_mean() == rhs.get_transformation_time_mean()
            && lhs.get_transformation_range() == rhs.get_transformation_range();
}

bool operator!=(const ind_param& lhs, const ind_param& rhs) noexcept
{
    return !(lhs==rhs);
}

ind_param change_ind_param_norm( ind_param i,  std::minstd_rand& rng)
{

    i.set_uptake_rate( std::normal_distribution<double>{i.get_uptake_mean(),
                                                        i.get_uptake_var()}(rng));

    i.set_spor_metabolic_rate(std::normal_distribution<double>{i.get_spor_metabolic_rate_mean(),
                                                               i.get_spor_metabolic_rate_var()}(rng));

    i.set_repr_prob(std::normal_distribution<double>{i.get_repr_prob_mean(),
                                                     i.get_repr_prob_var()}(rng));

    ///This is a binomial distribution as it is the discrete equivalent of the normal distribution
    /// an offset is required as the binomial alway goes from 0 to the number of trials
    auto rnd_numb =
            std::binomial_distribution<>{(i.get_transformation_range() * 2) + 1, 0.5}(rng);
    auto offset = i.get_transformation_time_mean() - i.get_transformation_range();
    i.set_transformation_time(rnd_numb + offset);

    return i;
}

ind_param change_ind_param_unif( ind_param i,  std::minstd_rand& rng)
{

    i.set_uptake_rate( std::uniform_real_distribution<double>{
                           i.get_uptake_mean() - 3 * i.get_uptake_var(),
                           i.get_uptake_mean() + 3 * i.get_uptake_var()
                       }(rng));

    i.set_spor_metabolic_rate(std::uniform_real_distribution<double>{
                                  i.get_spor_metabolic_rate_mean() - 3 * i.get_spor_metabolic_rate_var(),
                                  i.get_spor_metabolic_rate_mean() + 3 * i.get_spor_metabolic_rate_var()
                              }(rng));

    i.set_repr_prob(std::uniform_real_distribution<double>{
                        i.get_repr_prob_mean() - 3 * i.get_repr_prob_var(),
                        i.get_repr_prob_mean() + 3 * i.get_repr_prob_var()
                    }(rng));

    i.set_transformation_time(std::uniform_int_distribution<int>{
                                  i.get_transformation_time_mean() - i.get_transformation_range(),
                                  i.get_transformation_time_mean() + i.get_transformation_range()}
                              (rng)
                              );

    return i;
}

ind_param resize_ind_param(const ind_param& i, double amplitude)
{
    ind_param ind{
        i.get_radius(),
                i.get_treshold_energy(),
                i.get_uptake_rate(),
                i.get_uptake_mean(),
                i.get_uptake_var() * amplitude,
                i.get_metabolic_rate(),
                i.get_repr_prob(),
                i.get_repr_prob_mean(),
                i.get_repr_prob_var() * amplitude,
                i.get_spor_metabolic_rate(),
                i.get_spor_metabolic_rate_mean(),
                i.get_spor_metabolic_rate_var() * amplitude,
                i.get_transformation_time(),
                i.get_transformation_time_mean(),
                static_cast<int>( i.get_transformation_range() * amplitude),
                i.get_metab_secr_rate()
    };
    return ind;
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
        double radius  = 0.8;
        double treshold_energy = 10;
        double uptake_rate = 0.1;
        double uptake_rate_mean = 0.1;
        double uptake_rate_var = 0.02;
        double metabolic_rate = 0.01;
        double reproduction_prob = 0.5;
        double reproduction_prob_mean = 0.5;
        double reproduction_prob_var = 0.07;
        double spor_metabolic_rate = 0.5;
        double spor_metabolic_rate_mean = 0.5;
        double spor_metabolic_rate_var = 0.07;
        int transformation_time = 5;
        int transformation_time_mean = 5;
        int transformation_range = 2;
        double metab_secretion_rate = 1;
        ind_param p{
            radius,
                    treshold_energy,
                    uptake_rate,
                    uptake_rate_mean,
                    uptake_rate_var,
                    metabolic_rate,
                    reproduction_prob,
                    reproduction_prob_mean,
                    reproduction_prob_var,
                    spor_metabolic_rate,
                    spor_metabolic_rate_mean,
                    spor_metabolic_rate_var,
                    transformation_time,
                    transformation_time_mean,
                    transformation_range,
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

    //Ind parameters can be changed based one the
    //range related to the parameters
    //and the minimum fraction of the range for which they change
    {
        std::minstd_rand rng;
        ind_param i{};
        auto new_i = change_ind_param_norm(i,rng);
        assert( i != new_i);
    }

    ///Parameters can be changed both according a normal or uniform distribution
    ///If normal the variance members (m_var*) are use as the variance of the normal
    /// If uniform the variance members (m_var*) are used to determine the range
    /// range = [ m_mean - 3 * m_var, m_mean + 3 * m_var]
    {
        std::minstd_rand rng;
        double radius  = 0.8;
        double treshold_energy = 10;
        double uptake_rate = 0.1;
        double uptake_rate_mean = 0.1;
        double uptake_rate_var = 0.02;
        double metabolic_rate = 0.01;
        double reproduction_prob = 0.5;
        double reproduction_prob_mean = 0.5;
        double reproduction_prob_var = 0.07;
        double spor_metabolic_rate = 0.5;
        double spor_metabolic_rate_mean = 0.5;
        double spor_metabolic_rate_var = 0.07;
        int transformation_time = 5;
        int transformation_time_mean = 5;
        int transformation_range = 2;
        double metab_secretion_rate = 1;
        ind_param i_p_n{
            radius,
                    treshold_energy,
                    uptake_rate,
                    uptake_rate_mean,
                    uptake_rate_var,
                    metabolic_rate,
                    reproduction_prob,
                    reproduction_prob_mean,
                    reproduction_prob_var,
                    spor_metabolic_rate,
                    spor_metabolic_rate_mean,
                    spor_metabolic_rate_var,
                    transformation_time,
                    transformation_time_mean,
                    transformation_range,
                    metab_secretion_rate
        };

        int repeats = 1000;

        double normal_uptake_rate_mean = uptake_rate_mean;
        double normal_repr_prob_mean = reproduction_prob_mean;
        double normal_spor_metab_mean = spor_metabolic_rate_mean;
        double normal_transformation_time_mean = transformation_time;
        double unif_uptake_rate_mean = uptake_rate_mean;
        double unif_repr_prob_mean = reproduction_prob_mean;
        double unif_spor_metab_mean = spor_metabolic_rate_mean;
        double unif_transformation_time_mean = transformation_time;

        for( int i = 0; i != repeats; i++)
        {
            i_p_n = change_ind_param_norm(i_p_n,rng);
            normal_repr_prob_mean += i_p_n.get_repr_prob();
            normal_spor_metab_mean += i_p_n.get_spor_metabolic_rate();
            normal_uptake_rate_mean += i_p_n.get_uptake_rate();
            normal_transformation_time_mean += i_p_n.get_transformation_time();

            i_p_n = change_ind_param_unif(i_p_n,rng);
            unif_repr_prob_mean += i_p_n.get_repr_prob();
            unif_spor_metab_mean += i_p_n.get_spor_metabolic_rate();
            unif_uptake_rate_mean += i_p_n.get_uptake_rate();
            unif_transformation_time_mean += i_p_n.get_transformation_time();
        }

        normal_repr_prob_mean /= repeats;
        normal_spor_metab_mean /= repeats;
        normal_uptake_rate_mean /= repeats;
        normal_transformation_time_mean /= repeats;

        unif_repr_prob_mean /= repeats;
        unif_spor_metab_mean /= repeats;
        unif_uptake_rate_mean /= repeats;
        unif_transformation_time_mean /= repeats;

        assert(normal_repr_prob_mean - reproduction_prob_mean < 0.01
               && normal_repr_prob_mean - reproduction_prob_mean > -0.01);
        assert(normal_spor_metab_mean - spor_metabolic_rate < 0.01
               && normal_spor_metab_mean - spor_metabolic_rate > -0.01);
        assert( normal_uptake_rate_mean - uptake_rate_mean < 0.01
                &&  normal_uptake_rate_mean - uptake_rate_mean > -0.01);
        //Very coarse grained
        assert(normal_transformation_time_mean - transformation_time < 1
               && normal_transformation_time_mean - transformation_time > -1);

        assert(unif_repr_prob_mean - reproduction_prob_mean < 0.01
               && unif_repr_prob_mean - reproduction_prob_mean > -0.01);
        assert(unif_spor_metab_mean - spor_metabolic_rate < 0.01
               && unif_spor_metab_mean - spor_metabolic_rate > -0.01);
        assert( unif_uptake_rate_mean - uptake_rate_mean < 0.01
                &&  unif_uptake_rate_mean - uptake_rate_mean > -0.01);
        //Very coarse grained
        assert(unif_transformation_time_mean - transformation_time < 1
               && unif_transformation_time_mean - transformation_time > -1);
    }

    /// It is possible to create ind_param object
    /// starting from initial ones,
    /// that have a wider range of possible values
    {
        ind_param i;
        double amplitude = 1.5;

        auto i2 = resize_ind_param(i,amplitude);

        assert(i.get_repr_prob_var() - (i2.get_repr_prob_var() / amplitude) < 0.00001
               && i.get_repr_prob_var() - (i2.get_repr_prob_var() / amplitude) > -0.00001);
        assert(i.get_uptake_var() - (i2.get_uptake_var() / amplitude) < 0.00001
               && i.get_uptake_var() - (i2.get_uptake_var() / amplitude) > -0.00001);
        assert(i.get_spor_metabolic_rate_var() - (i2.get_spor_metabolic_rate_var() / amplitude) < 0.00001
               && i.get_spor_metabolic_rate_var() - (i2.get_spor_metabolic_rate_var() / amplitude) > -0.00001);
        assert(i.get_transformation_range() - (i2.get_transformation_range() / amplitude) < 0.00001
               && i.get_transformation_range() - (i2.get_transformation_range() / amplitude) > -0.00001);
    }
}

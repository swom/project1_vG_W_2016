#include"tests.h"

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

        auto i2 = change_range_ind_param(i,amplitude);

        assert(i.get_repr_prob_var() - (i2.get_repr_prob_var() / amplitude) < 0.00001
               && i.get_repr_prob_var() - (i2.get_repr_prob_var() / amplitude) > -0.00001);
        assert(i.get_uptake_var() - (i2.get_uptake_var() / amplitude) < 0.00001
               && i.get_uptake_var() - (i2.get_uptake_var() / amplitude) > -0.00001);
        assert(i.get_spor_metabolic_rate_var() - (i2.get_spor_metabolic_rate_var() / amplitude) < 0.00001
               && i.get_spor_metabolic_rate_var() - (i2.get_spor_metabolic_rate_var() / amplitude) > -0.00001);
        assert(i.get_transformation_range() - ( i2.get_transformation_range() / static_cast<int>(amplitude)) < 0.00001
               && i.get_transformation_range() - (i2.get_transformation_range() / static_cast<int>(amplitude)) > -0.00001);
    }
}

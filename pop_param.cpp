#include "pop_param.h"

#include <cassert>
#include <iostream>
#include <string>

pop_param::pop_param(unsigned int start_pop_size,
                     unsigned int exp_new_pop_size,
                     double min_dist,
                     double mutation_prob,
                     double mutation_step,
                     double base_disp_prob,
                     double spore_advantage,
                     double death_rate,
                     class ind_param ind_parameters):
    m_ind_param{ind_parameters},
    m_base_disp_prob{base_disp_prob},
    m_death_prob{death_rate},
    m_exp_new_pop_size{exp_new_pop_size},
    m_min_init_dist_btw_inds{min_dist},
    m_mutation_prob{mutation_prob},
    m_mutation_step{mutation_step},
    m_spore_advantage{spore_advantage},
    m_start_pop_size{start_pop_size}
{
#ifndef IS_ON_TRAVIS
    try {
        if(base_disp_prob * spore_advantage > 1)
        {throw std::string{"base dispersal probability * spore advantage > 1, too high!\n"};}
    }
    catch (std::string& e) {
        std::cout << e;
#ifdef NDEBUG
        abort();
#endif
    }
    assert(m_base_disp_prob > -0.000001);
    assert(m_exp_new_pop_size >= 0);
    assert(m_min_init_dist_btw_inds > -0.000001);
    assert(m_mutation_prob > -0.000001 && m_mutation_prob < 1.000001);
    assert(m_mutation_step > -0.000001 && m_mutation_step < 1.000001);
    assert(m_spore_advantage > -0.000001);
    assert(m_start_pop_size >= 0);
#endif

}


std::ostream& operator<<(std::ostream& os, const pop_param& p)
{
    os << p.get_mu_p() << " , "
       << p.get_mu_st() << " , "
       << p.get_spo_adv() << " , "
       << p.get_min_dist() << " , "
       << p.get_death_rate() << " , "
       << p.get_base_disp_prob() << " , "
       << p.get_pop_start_size() << " , "
       << p.get_exp_new_pop_size();
    return os;
}

std::ifstream& operator >>(std::ifstream& is, pop_param& p)
{
    unsigned int start_pop_size;
    unsigned int exp_new_pop_size;
    double min_dist;
    double mutation_prob;
    double mutation_step;
    double base_disp_prob;
    double spore_advantage;
    double death_rate;
    std::string dummy; // To remove the annotation in the file

    is      >> mutation_prob >> dummy
            >> mutation_step >> dummy
            >> spore_advantage >> dummy
            >> min_dist >> dummy
            >> death_rate >> dummy
            >> base_disp_prob >> dummy
            >> start_pop_size >> dummy
            >> exp_new_pop_size;

    p =     pop_param {start_pop_size,
            exp_new_pop_size,
            min_dist,
            mutation_prob,
            mutation_step,
            base_disp_prob,
            spore_advantage,
            death_rate
};

    return is;
}

bool operator==(const pop_param& lhs, const pop_param& rhs) noexcept
{
    return
            lhs.get_mu_p() == rhs.get_mu_p()
            && lhs.get_mu_st() == rhs.get_mu_st()
            && lhs.get_spo_adv() == rhs.get_spo_adv()
            && lhs.get_min_dist() == rhs.get_min_dist()
            && lhs.get_death_rate() == rhs.get_death_rate()
            && lhs.get_base_disp_prob() == rhs.get_base_disp_prob()
            && lhs.get_pop_start_size() == rhs.get_pop_start_size()
            && lhs.get_exp_new_pop_size() == rhs.get_exp_new_pop_size()
            ;
}

bool operator!=(const pop_param& lhs, const pop_param& rhs) noexcept
{
    return
            !(lhs == rhs);
}

void save_pop_parameters(
        const pop_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p;
}

pop_param load_pop_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    pop_param p;
    f >> p;
    return p;
}

void test_pop_param() noexcept  //!OCLINT
{

    //At initialization a pop checks that base_disp_dist * spore advantage is not > 1
    //--------> constructor throws exception. Tested directly in constructor when in release mode
    {
#ifndef IS_ON_TRAVIS
        try {
            pop_param(0,0,0,0,0,1,10);
        } catch (std::string& e) {
            assert(e == "base dispersal probability * spore advantage > 1, too high!\n" );
        }
#endif
    }

    //A pop_param has a member variable m_new_pop_size that states the max number of
    //individuals that will start a new population

    {
        unsigned int exp_pop_size = 100;
        pop_param p{1, exp_pop_size};
        assert(p.get_exp_new_pop_size() == exp_pop_size);
    }

    //It is possible to save and load ind parameters from a file
    {
        unsigned int start_pop_size = 1;
        unsigned int exp_new_pop_size = 1;
        double min_dist = 0.1;
        double mutation_prob = 0.0015;
        double mutation_step = 0.1;
        double base_disp_prob = 0.01;
        double spore_advantage = 10.0;
        double reproduction_prob = 0.5;
        double death_rate = 0.0;
        pop_param p{
            start_pop_size,
                    exp_new_pop_size,
                    min_dist,
                    mutation_prob,
                    mutation_step,
                    base_disp_prob,
                    spore_advantage,
                    reproduction_prob,
                    death_rate
        };

        //Test load and save
        const std::string filename = "pop_param.csv";
        save_pop_parameters(p, filename);
        const pop_param q = load_pop_parameters(filename);
        assert(p == q);
        //Test >> operator overload
        std::ifstream f(filename);
        pop_param s;
        f >> s;
        assert(s == p);
    }

    {
        const ind_param p;
        std::ostringstream s;
        s << p;
    }
}

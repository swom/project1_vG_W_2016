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
                     double death_rate):
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
        std::cout << e << std::endl;
#ifdef NDEBUG
        abort();
#endif
    }
    assert(m_base_disp_prob > -0.000001);
    assert(m_min_init_dist_btw_inds > -0.000001);
    assert(m_mutation_prob > -0.000001 && m_mutation_prob < 1.000001);
    assert(m_mutation_step > -0.000001 && m_mutation_step < 1.000001);
    assert(m_spore_advantage > -0.000001);

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

void save_pop_parameters_json(
        const pop_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    nlohmann::json json_out;
    json_out = p;
    f << json_out;
}

pop_param load_pop_parameters_json(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    pop_param p;
    nlohmann::json json_in;
    f >> json_in;
    p = json_in;
    return p;
}


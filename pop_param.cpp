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
                     double reproduction_prob ,
                     class ind_param ind_parameters):
    m_ind_param{ind_parameters},
    m_base_disp_prob{base_disp_prob},
    m_exp_new_pop_size{exp_new_pop_size},
    m_min_init_dist_btw_inds{min_dist},
    m_mutation_prob{mutation_prob},
    m_mutation_step{mutation_step},
    m_repr_prob{reproduction_prob},
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
#endif
    assert(m_base_disp_prob > 0);
    assert(m_exp_new_pop_size > 0);
    assert(m_min_init_dist_btw_inds > 0);
    assert(m_mutation_prob > 0 && m_mutation_prob < 1.000001);
    assert(m_mutation_step > 0 && m_mutation_step < 1.000001);
    assert(m_repr_prob > 0 && m_repr_prob < 1.0000001);
    assert(m_spore_advantage > 0);
    assert(m_start_pop_size >0);
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
}

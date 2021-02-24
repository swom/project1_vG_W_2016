#include"tests.h"
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
        double death_rate = 0.0;
        pop_param p{
            start_pop_size,
                    exp_new_pop_size,
                    min_dist,
                    mutation_prob,
                    mutation_step,
                    base_disp_prob,
                    spore_advantage,
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

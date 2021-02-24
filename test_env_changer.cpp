#include "test_env_changer.h"
#include <cassert>

///An env_changer is constructed with env_param and an amplitude
void test_e_c_constructor(){

    env_param intial_mean_params{};
    double amplitude = 1.2345456;
    env_changer ec{intial_mean_params, amplitude};
    assert(ec.get_amplitude() == amplitude);
    assert(ec.get_mean_params() == intial_mean_params);
}

///An env_changer can produce an env_param object that
/// has parametes value drawn from an uniform distribution
/// determined by the mean and the amplitude
void test_change() {

    double amplitude = 0.2;
    env_changer ec{env_param{}, amplitude};

//    change
}

///The maximum amplitude can be 1
void test_legal_amplitude(){
    double amplitude = 2;
    try {
        env_changer{env_param{}, amplitude};
    } catch (std::string& e) {
        assert(e == "amplitude above max_limit" );
    }
}

/// amplitude 1 means that the lowest limit of the distribution is 0
void test_max_amplitude(){
    double amplitude = 1;
    env_changer ec{env_param{}, amplitude};
    int repeats = 10000;
    for(int i = 0; i!= repeats; i++)
    {
//    change_env(ec);
    }
}
void test_env_changer(){

    test_e_c_constructor();
    test_legal_amplitude();
    test_max_amplitude();
    test_change();

}

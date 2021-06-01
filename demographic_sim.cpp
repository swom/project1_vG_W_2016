#include "demographic_sim.h"
#include<cassert>

demographic_sim::demographic_sim()
{

}

bool operator==(const demographic_sim& lhs, const demographic_sim& rhs) noexcept
{
    return lhs.get_demo_cycles() == rhs.get_demo_cycles();
}

bool operator!=(const demographic_sim& lhs, const demographic_sim& rhs) noexcept
{
    return !(lhs == rhs);
}


std::string create_sim_demo_name(int seed, int change_freq, std::string prefix, std::string suffix)
{
    return  std::string{
        prefix +
                "sim_demographic_s" +
                std::to_string(seed) +
                "change_" +
                std::to_string(change_freq) +
                suffix +
                ".csv"
    };
}
demographic_sim load_demographic_sim(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    if(!f.is_open())
    {
        std::cout << "Could not find specified demographic_sim*.csv file: "
        << filename
        << "n\n";
        abort();
    }
    demographic_sim d_s;
    demographic_cycle d_c{0,
                          0,
                          0,
                          0,
                          env_param{}, ind_param{}};
    std::string dummy;
    while(f >> dummy)//Skips the cycle number
    {
        f >> dummy; //skips the comma
        f >> d_c;
        d_s.get_demo_cycles().push_back(d_c);
    }

    return d_s;
}

void save_demographic_sim(
        const demographic_sim& d_s,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    int cycle_number = 0;
    for(const auto& demo_cycle : d_s.get_demo_cycles())
    {
        f << cycle_number << " , " << demo_cycle ;
        cycle_number++;
    }
}

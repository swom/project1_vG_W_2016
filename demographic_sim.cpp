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


std::string create_sim_demo_name(int seed, int change_freq)
{
    return  std::string{
        "sim_demographic_s" +
        std::to_string(seed) +
                "change_" +
                std::to_string(change_freq) +
                ".json"
    };
}

demographic_sim load_demographic_sim(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    if(!f.is_open())
    {
        std::cout << "Could not find specified demographic_sim*.json file. \n\n";
        abort();
    }
    demographic_sim d_s;
    nlohmann::json json_in;
    f >> json_in;
    d_s = json_in;

    return d_s;
}

void save_demographic_sim(
        const demographic_sim& d_s,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    nlohmann::json json_out;
    json_out = d_s;
    f << json_out;
}

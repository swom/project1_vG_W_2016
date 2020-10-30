#include "utilities.h"
#include <sys/stat.h>

bool compare_with_tolerance(const std::vector<double>& lhs,const std::vector<double>& rhs)
{
    return std::equal(lhs.begin(), lhs.end(), rhs.begin(),
                      [](const double& lhs_w, const double& rhs_w)
    {return lhs_w - rhs_w < 0.0001 && lhs_w - rhs_w > -0.00001;});
}

void check_for_cmd_param(const std::vector<std::string>& args,
                         int& seed,
                         int& change_freq,
                         int& n_conditions,
                         int& replay_cycle,
                         double& amplitude,
                         bool& overwrite)
{
    if ((args.size() > 3
            && (args[1] == "--sim"
                || args[1] == "--rand"
                || args[1] == "--rand_best"
                || args[1] == "--reac_norm"
                || args[1] == "--replay")
            )
            ||
            (args.size() > 1
             &&
             args[1] == "--create_rand_cond_vec"))
    {
        take_amplitude_arg(args, amplitude);
        take_change_freq_arg(args,change_freq);
        take_n_conditions_arg(args, n_conditions);
        take_seed_arg(args, seed);
        take_overwrite_arg(args, overwrite);
        take_replay_cycle_arg(args,replay_cycle);
    }
}


std::uniform_real_distribution<double> create_unif_dist(double a, double b) noexcept
{
    return std::uniform_real_distribution<double>(a,b);
}

std::bernoulli_distribution create_bernoulli_dist(double p) noexcept
{
    return std::bernoulli_distribution{p};
}

std::normal_distribution<double> create_normal_dist(double m, double v)
{
    return std::normal_distribution<double>{m, v};
}

bool exists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

void take_amplitude_arg(const std::vector<std::string>& args, double& amplitude)
{
    if (args.size() > 4
            &&  (args[1] == "--rand"
                 || args[1] == "--rand_best")
            &&  args[4][0] == 'a'         
            )
    {
        std::string s_amplitude;
        for(size_t i = 1; i != args[4].size(); i++)
        {
            s_amplitude += args[4][i];
        }
        amplitude = std::stod(s_amplitude);
    }
    else if (args.size() > 1
              &&
              args[1] == "--create_rand_cond_vec"
              &&
              args[2][0] == 'a'
              )
    {
        std::string s_amplitude;
        for(size_t i = 1; i != args[2].size(); i++)
        {
            s_amplitude += args[2][i];
        }
        amplitude = std::stod(s_amplitude);
    }

}

void take_n_conditions_arg(const std::vector<std::string>& args, int& n_conditions)
{
    if ((args.size() > 5
            &&  (args[1] == "--rand"
                 || args[1] == "--rand_best")
            &&  args[5][0] == 'n')
            ||
            (args.size() > 2
             &&
             args[1] == "--create_rand_cond_vec"
             &&
             args[3][0] == 'n')
            )
    {
        std::string s_n_random_conditions;
        for(size_t i = 1; i != args[5].size(); i++)
        {
            s_n_random_conditions += args[5][i];
        }
        n_conditions = std::stoi(s_n_random_conditions);
    }
}

void take_change_freq_arg(const std::vector<std::string>& args, int& change_freq)
{

    if(args[3][0] == 'f' && std::isdigit(args[3][1]))
    {
        std::string s_change_freq;
        for(size_t i = 1; i != args[3].size(); i++)
        {
            s_change_freq += args[3][i];
        }
        change_freq = std::stoi(s_change_freq);
    }
    else
    {
        std::cout << "Invalid fourth argument: it has to be an fN\n"
                     " Where N is the number of cycles after which\n"
                     "the environment changes";
        abort();
    }

}

void take_replay_cycle_arg(const std::vector<std::string>& args, int& replay_cycle)
{
    if (args.size() > 3 &&
            (args[1] == "--replay")
            )
    {
        if(args[4][0] == 'c' && std::isdigit(args[4][1]))
        {
            std::string s_replay_cycle;
            for(size_t i = 1; i != args[4].size(); i++)
            {
                s_replay_cycle += args[4][i];
            }
            replay_cycle = std::stoi(s_replay_cycle);
        }
        else
        {
            std::cout << "Invalid fifth argument: if cmd line starts with --replay\n "
                         "this argument has to be cN:Where N is the cycle which is going\n"
                         "to be replayed";
            abort();
        }
    }
}

void take_overwrite_arg(const std::vector<std::string>& args, bool& overwrite)
{
    if (args.back() == "--overwrite")
        overwrite = true;
    else
        overwrite = false;
}

void take_seed_arg(const std::vector<std::string>& args, int& seed)
{
    if(args[2][0] == 's' && std::isdigit(args[2][1]))
    {
        std::string s_seed;
        for(size_t i = 1; i != args[2].size(); i++)
        {
            s_seed += args[2][i];
        }
        seed = std::stoi(s_seed);
    }
    else
    {
        std::cout << "Invalid third argument: it has to be an sN\n"
                     " Where N is the seed number";
        abort();
    }
}

void test_utilities() noexcept //!OCLINT
{

}

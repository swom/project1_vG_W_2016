#include "utilities.h"
#include <sys/stat.h>
#include <chrono>
#include <cassert>

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
                         bool& overwrite,
                         int& seed_rand_cond,
                         int& rand_cond_n)
{
    if (args.size() > 3
            && (args[1] == "--sim"
                || args[1] == "--rand"
                || args[1] == "--rand_best"
                || args[1] == "--reac_norm"
                || args[1] == "--replay"
                || args[1] == "--continue_evo"
                || args[1] ==  "--replay_rand_cond"
                || args[1] == "--replay_best_rand_cond")
            )
    {
        take_amplitude_arg(args, amplitude);
        take_change_freq_arg(args,change_freq);
        take_n_conditions_arg(args, n_conditions);
        take_seed_arg(args, seed);
        take_overwrite_arg(args, overwrite);
        take_replay_cycle_arg(args,replay_cycle);
        take_seed_rand_cond(args,seed_rand_cond);
        take_rand_cond_n(args, rand_cond_n);

    }
    else if(args.size() > 1
            &&
            args[1] == "--create_rand_cond_vec")
    {
        take_amplitude_arg(args, amplitude);
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
                 || args[1] == "--rand_best"
                 || args[1] == "--replay_rand_cond"
                 || args[1] == "--replay_best_rand_cond")
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
              || args[1] == "--rand_best"
              || args[1] ==  "--replay_rand_cond"
              || args[1] == "--replay_best_rand_cond")
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

void take_rand_cond_n(const std::vector<std::string>& args, int& rand_cond_n)
{
    if((args[1] == "--replay_rand_cond"
            || args[1] == "--replay_best_rand_cond")&&
            args.size() > 6 &&
            args[7][0] == 'r' &&
            args[7][1] == 'n'
            )
    {
        std::string s_rand_cond_n;
        for(size_t i = 2; i != args[7].size(); i++)
        {
            s_rand_cond_n += args[7][i];
        }
        rand_cond_n = std::stoi(s_rand_cond_n);
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

void take_seed_rand_cond(const std::vector<std::string>& args, int& seed_rand_cond)
{
    if((args[1] == "--replay_rand_cond"
            || args[1] == "--replay_best_rand_cond") &&
            args.size() > 5 &&
            args[6][0] == 's' &&
            args[6][1] == 'r')
    {
        std::string s_seed_rand_cond;
        for(size_t i = 2; i != args[6].size(); i++)
        {
            s_seed_rand_cond += args[6][i];
        }
        seed_rand_cond = std::stoi(s_seed_rand_cond);
    }
}

void test_utilities() noexcept //!OCLINT
{

    ///The seed of the random random number generator
    /// of the simulation can be taken as an argument
    /// from a commnad line
    {

        int seed = 123;
        int expected_seed = 0;
        assert(seed != expected_seed);

        //create commnad line where seed of random condition is 0
        auto seed_s = "s" + std::to_string(expected_seed);
        const std::vector<std::string>& args{
                "executable's_folder",
                "--replay_rand_cond",
                seed_s,
                "freq",
                "ampl",
                "n_rand_cond",
                "seed_random_conditions",
                "rand_cond_n"};

        take_seed_arg(args, seed);
        assert(seed == expected_seed);
    }

    ///The frequency with which the enviornment changes
    /// can be taken as an argument
    /// from a commnad line
    {

        int change_freq = 123;
        int expected_change_freq = 0;
        assert(change_freq != expected_change_freq);

        //create commnad line where seed of random condition is 0
        auto change_freq_s = "f" + std::to_string(expected_change_freq);
        const std::vector<std::string>& args{
                "executable's_folder",
                "--replay_rand_cond",
                "seed",
                change_freq_s,
                "ampl",
                "n_rand_cond",
                "seed_random_conditions",
                "rand_cond_n"};

        take_change_freq_arg(args, change_freq);
        assert(change_freq == expected_change_freq);
    }

    ///The frequency with which the enviornment changes
    /// can be taken as an argument
    /// from a commnad line
    {

        double amplitude = 123;
        double expected_amplitude = 0;
        assert(amplitude != expected_amplitude);

        //create commnad line where seed of random condition is 0
        auto amplitude_s = "a" + std::to_string(expected_amplitude);
        const std::vector<std::string>& args{
                "executable's_folder",
                "--replay_rand_cond",
                "seed",
                "change_freq",
                amplitude_s,
                "n_rand_cond",
                "seed_random_conditions",
                "rand_cond_n"};

        take_amplitude_arg(args, amplitude);
        assert(amplitude == expected_amplitude);
    }

    ///The frequency with which the enviornment changes
    /// can be taken as an argument
    /// from a commnad line
    {

        int n_rand_cond = 123;
        int expected_n_rand_cond = 0;
        assert(n_rand_cond != expected_n_rand_cond);

        //create commnad line where seed of random condition is 0
        auto n_rand_cond_s = "n" + std::to_string(expected_n_rand_cond);
        const std::vector<std::string>& args{
                "executable's_folder",
                "--replay_rand_cond",
                "seed",
                "change_freq",
                "ampl",
                n_rand_cond_s,
                "seed_random_conditions",
                "rand_cond_n"};

        take_n_conditions_arg(args, n_rand_cond);
        assert(n_rand_cond == expected_n_rand_cond);
    }

    ///The seed of the reandom conditions can be takenas an argument
    /// from a commnad line
    {

        int seed_rand_cond = 123;
        int expected_seed_rand_cond = 0;
        assert(seed_rand_cond != expected_seed_rand_cond);

        //create commnad line where seed of random condition is 0
        auto seed_rand_cond_s = "sr" + std::to_string(expected_seed_rand_cond);
        const std::vector<std::string>& args{
                "executable's_folder",
                "--replay_rand_cond",
                "seed",
                "freq",
                "ampl",
                "n_rand_cond",
                seed_rand_cond_s,
                "rand_cond_n"};

        take_seed_rand_cond(args, seed_rand_cond);
        assert(seed_rand_cond == expected_seed_rand_cond);
    }

    ///The index of the random condition that need to be selcted
    /// from the random condiiton vector can be taken as an argument
    /// from a commnad line
    {

        int rand_cond_n = 123;
        int expected_rand_cond_n = 0;
        assert(rand_cond_n != expected_rand_cond_n);

        //create commnad line where seed of random condition is 0
        auto rand_cond_n_s = "rn" + std::to_string(expected_rand_cond_n);
        const std::vector<std::string>& args{
                "executable's_folder",
                "--replay_rand_cond",
                "seed",
                "freq",
                "ampl",
                "n_rand_cond",
                "seed_random_conditions",
                rand_cond_n_s};

        take_rand_cond_n(args, rand_cond_n);
        assert(rand_cond_n == expected_rand_cond_n);
    }



}

#include "tests.h"

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
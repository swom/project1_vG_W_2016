#include "tests.h"

void test_GRN()//!OCLINT , tests may be long
{
#ifndef NDEBUG
    //A grn is initialized with connection weights
    //all to a given value, default  = 1
    {
        double weight_value = 3.14;
        GRN g{3,3,2,weight_value};
        for(const auto& node : g.get_I2H())
        {
            for(const auto& weight : node)
                assert(weight - weight_value < 0.00001 && weight - weight_value > -0.00001);
        }
        for(const auto& node : g.get_H2H())
        {
            for(const auto& weight : node)
                assert(weight - weight_value < 0.00001 && weight - weight_value > -0.00001);
        }
        for(const auto& node : g.get_H2O())
        {
            for(const auto& weight : node)
                assert(weight - weight_value < 0.00001 && weight - weight_value > -0.00001);
        }
    }
    //A GRN is initialized with three layers
    // and a certain amount of nodes in each layer
    //These hyperparameters cannot be altered
    //during the simulation therefore no set() functions.
    //By default 3-3-2
    {
        size_t n_input = 5;
        size_t n_hidden = 5;
        size_t n_output = 5;
        GRN g(n_input,n_hidden,n_output);
        assert(g.get_I2H().size() == n_input);
        assert(g.get_H2H().size() == n_hidden);
        assert(g.get_H2O().size() == n_output);
    }

    //It is possible to get the values of the nodes states
    {
        GRN g;
        for(auto input : g.get_input_nodes())
            assert(input < 0.00001 && input > -0.00001);
        for(auto hidden : g.get_hidden_nodes())
            assert(hidden);
        for(auto output : g.get_output_nodes())
            assert(output);
    }

    //All weights in H2O connection can be set to a value
    {
        GRN g;
        double weight_value = 3.14;
        g.set_all_H2O(weight_value);
        for(const auto& node : g.get_H2O())
            for(const auto& weight : node)
                assert(weight - weight_value < 0.00001 &&
                       weight - weight_value > -0.00001);
    }

    //All weights in H2H connection can be set to a value
    {
        GRN g;
        double weight_value = 3.14;
        g.set_all_H2H(weight_value);
        for(const auto& node : g.get_H2H())
            for(const auto& weight : node)
                assert(weight - weight_value < 0.00001 &&
                       weight - weight_value > -0.00001);
    }

    //All output nodes can be set to true or false
    {
        GRN g;
        bool nodes_state = false;
        g.set_all_out_nodes(nodes_state);
        for(auto node : g.get_output_nodes())
            assert(node == nodes_state);
    }
    //All output nodes treshold can be set to a value
    {
        GRN g;
        double treshold_value = 3.14;
        g.set_all_out_tresh(treshold_value);
        for (const auto& treshold : g.get_out_tresh())
        {
            assert(treshold - treshold_value < 0.00001 &&
                   treshold - treshold_value > -0.00001);
        }
    }
    //All hidden nodes can be set to true or false
    {
        GRN g;
        bool nodes_state = false;
        g.set_all_hid_nodes(nodes_state);
        for(auto node : g.get_hidden_nodes())
            assert(node == nodes_state);
    }

    //Output nodes can be updated based on values of hidden nodes
    {
        GRN g;
        bool hidden_nodes_value = true;
        //all hidden nodes are 1
        g.set_all_hid_nodes(hidden_nodes_value);
        //all weights H20 are
        //hidden_nodes_value / number_of_connections_per_output_node
        double weight_value = hidden_nodes_value;
        weight_value /= g.get_H2O().size();

        g.set_all_H2O(weight_value);
        //Therefore each output node will receive a signal of hidden_nodes_value
        g.set_all_out_tresh(hidden_nodes_value - 0.00001);
        hid_updates_out(g);
        for (size_t i = 0; i != g.get_output_nodes().size(); i++)
            assert(g.get_output_nodes()[i]);
    }

    //Given the initial inputs the hidden nodes will recieve a signal
    {
        GRN g;
        double inputs_value = 3.14; //in this case all inputs are the same
        auto inputs = std::vector<double>(g.get_input_nodes().size(),inputs_value);
        g.set_inputs(inputs);
        g.set_all_hid_tresh(inputs_value - 0.00001);
        //all weight I2H are
        //1 / number_of_connections_of each_hidden_to_input
        g.set_all_I2H( 1.0 / g.get_I2H().size() );
        auto inp_signals = inp_updates_hid(g);
        //Therefore each hidden node will receive a signal of hidden_nodes_value
        for (size_t i = 0; i != inp_signals.size(); i++)
            assert(inp_signals[i] - inputs_value < 0.00001 &&
                   inp_signals[i] - inputs_value > -0.00001);

    }

    //  //The values of the hidden nodes are stored in a vector to be used in
    //  //the successive timestep to auto-update the hidden nodes
    ////At initialization the values in this vector are all 0
    //  {
    //    GRN g;
    //    for (const auto& value : g.get_past_hid_val())
    //      {
    //        assert(value < 0.000001 && value > -0.000001)
    //      }
    //    assert(store_hid_values(g) != g.set_all_hid_nodes(1));
    //  }
    //The hidden nodes updated themselves toghether with input
    {
        GRN g;
        bool hidden_value = true; //in this case all nodes are the same
        g.set_all_hid_nodes(hidden_value);
        //all weight H2H are
        //1 / number_of_connection_to_each_hidden_node
        g.set_all_H2H( 1.0 / g.get_H2H().size());
        auto hid_signals = hid_updates_hid(g);
        //Therefore each hidden node will receive a signal == hidden_value
        for (size_t i = 0; i != hid_signals.size(); i++)
            assert(hid_signals[i] - static_cast<double>(hidden_value) < 0.00001 &&
                   hid_signals[i] - static_cast<double>(hidden_value) > -0.00001);
        //if instead the weights are smaller signal != hidden_value
        g.set_all_H2H(1.0 / g.get_H2H().size() * 2);
        hid_signals = hid_updates_hid(g);
        for (size_t i = 0; i != hid_signals.size(); i++)
            assert( hid_signals[i] - static_cast<double>(hidden_value) > 0.00001 ||
                    hid_signals[i] - static_cast<double>(hidden_value) < -0.00001);

    }

    //Given three initial input a GRN gives back an output
    {
        GRN g;//the weights and connection and nodes are all to 0, so output should be 0
        std::vector<double> inputs{100,100,100};
        g.set_inputs(inputs);
        jordi_response_mech(g);
        for(auto output : g.get_output_nodes())
        {
            assert(output);
        }
    }

    //The number of connections of a network can be counted
    {
        GRN g;
        assert(n_connections(g) == static_cast<int> (
                   g.get_H2H().size() * g.get_H2H()[0].size() +
               g.get_H2O().size() * g.get_H2O()[0].size() +
                g.get_I2H().size() * g.get_I2H()[0].size()
                )
                );

    }

    //The weights of a network can be summed toghether
    {
        GRN g;
        g.set_all_H2O(1);
        g.set_all_I2H(1);
        g.set_all_H2H(1);
        auto sum = weights_sum(g);
        auto n_con = static_cast<double>(n_connections(g));
        assert(sum - n_con < 0.00001 &&
               sum - n_con > -0.00001);
    }

    //The mean of the weights of the network can be calculated
    {
        GRN g;
        double weights_value = 1;
        g.set_all_H2O(weights_value);
        g.set_all_I2H(weights_value);
        g.set_all_H2H(weights_value);
        assert(weights_mean(g) - weights_value < 0.00001 &&
               weights_mean(g) - weights_value > -0.00001);
    }

    //The variance of weights of a network can be calculated
    {
        GRN g;
        double weights_value = 1;
        //All weights are the same so variance should be 0
        g.set_all_H2O(weights_value);
        g.set_all_I2H(weights_value);
        g.set_all_H2H(weights_value);
        assert(weights_var(g) < 0.00001 &&
               weights_var(g) > -0.00001);
    }

    //A GRN can be printed and loaded from an ifstream and ofstream file
    {
        GRN g{1,2,1};
        g.set_all_H2O(0.12);
        g.set_all_I2H(0.123);
        g.set_all_H2H(0.124);
        const std::string filename = "grn.csv";
        save_grn(g, filename);
        GRN g1;
        g1 = load_grn(filename);
        assert( g == g1);
    }

    ///The reaction norm of a grn can be calculated
    /// given the range of values for the inputs
    /// (always starting from 0 to a max value)
    /// and the step of change/resolution of points
    /// It will return all points for whihc the network would start sporulation
    {
        GRN g;
        double size = 2;
        double max_en = size;
        double max_food = size;
        double max_metab = size;
        double step = 0.2;
        int seed = 123;
        int change_freq = 12;

        auto reac_norm = calc_reaction_norm(g,
                                            max_en,
                                            max_food,
                                            max_metab,
                                            step);

        //The size of the reaction norm will
        //always be less or equal to all possible combinations
        assert(reac_norm.size() <= pow(size / step, 3));

        std::string filename = create_reaction_norm_name( seed,  change_freq);

        save_reaction_norm(reac_norm, filename);
        auto loaded_reac_norm = load_reaction_norm(filename);


        for(size_t reaction = 0; reaction  != reac_norm.size(); reaction++)
            for(size_t value = 0; value != reac_norm[reaction].size(); value++)
                assert(reac_norm[reaction][value] - loaded_reac_norm[reaction][value] < 0.000001
                       && reac_norm[reaction][value] - loaded_reac_norm[reaction][value] > -0.000001);
    }

    ///A reaction norm can store multiple responses to the smae set of inputs
    /// the multiple responses are the output of the network after
    /// it has responded n times to a certain set of input
    /// excluding the first response that is not dictated by the inputs,
    /// but by the starting internal state of the input
    {

        GRN g;
        double size = 2;
        double max_en = size;
        double max_food = size;
        double max_metab = size;
        double step = 0.2;
        int seed = 1234;
        int change_freq = 12;
        int n_responses = 2;

        auto reac_norm = calc_reaction_norm(g,
                                            max_en,
                                            max_food,
                                            max_metab,
                                            step,
                                            n_responses);

        auto filename = create_reaction_norm_name(seed, change_freq);
        save_reaction_norm(reac_norm, filename);
        auto loaded_reac_norm = load_reaction_norm(filename);

        for(size_t condition = 0; condition  != reac_norm.size(); condition++)
            for(size_t value = 0; value != reac_norm[condition].size(); value++)
            {
                //assert each condition vector size is equal
                //to the number of inputs taken (3 HARDCODED!!!!) + the number of responses
                assert(static_cast<int>(reac_norm[condition].size()) == 3 + n_responses);
                assert(reac_norm[condition][value] - loaded_reac_norm[condition][value] < 0.000001
                       && reac_norm[condition][value] - loaded_reac_norm[condition][value] > -0.000001);
            }
    }
#endif

}

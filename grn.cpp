#include "grn.h"
#include<cassert>


GRN::GRN(size_t n_input, size_t n_hidden, size_t n_output):
  m_ConI2H(n_input,std::vector<double>(n_hidden,0.0)),
  m_ConH2H(n_hidden,std::vector<double>(n_hidden,0.0)),
  m_ConH2O(n_hidden,std::vector<double>(n_output,0.0)),
  m_THidden(n_hidden,0),
  m_TOutput(n_output,0),
  m_ExInput(n_input,0),
  m_ExHidden(n_hidden,0),
  m_ExOutput(n_output,0)
{

}

void GRN::set_hid_node(int index_node, bool state)
{m_ExHidden[static_cast<unsigned int>(index_node)] = state;}

void GRN::set_all_H2O(double value) noexcept
{
  for (auto& node : get_H2O())
    {
      for(auto& connection_weight : node)
      connection_weight = value;
    }
}

void GRN::set_all_I2H(double value) noexcept
{
  for (auto& node : get_I2H())
    {
      for(auto& connection_weight : node)
      connection_weight = value;
    }
}

void GRN::set_all_out_nodes(bool state)
{
  for(int i = 0; i != static_cast<int>(get_output_nodes().size()); i++)
    set_out_node(i, state);
}


void GRN::set_all_hid_nodes(bool state)
{
  for(int i = 0; i != static_cast<int>(get_hidden_nodes().size()); i++)
    set_hid_node(i,state);
}

void GRN::set_all_hid_tresh(double value)
{
  for(auto& treshold : get_hid_tresh())
     treshold = value;
}

void GRN::set_all_out_tresh(double value)
{
  for(auto& treshold : get_out_tresh())
     treshold = value;
}

void GRN::set_inputs(std::vector<double> inputs)
{
  assert(inputs.size() == get_input_nodes().size());
  get_input_nodes().swap(inputs);
}

void GRN::set_out_node(int index_node, bool state)
{m_ExOutput[static_cast<unsigned int>(index_node)] = state;}

void inp_updates_hid(GRN& g) noexcept
{
  for(unsigned int i = 0; i != g.get_hidden_nodes().size(); i++)
    {
      double signal_hid = 0;
      for (unsigned int j = 0 ; j != g.get_input_nodes().size(); j++)
        {
          signal_hid += g.get_input_nodes()[i] * g.get_I2H()[j][i];
        }
      if(signal_hid > g.get_hid_tresh()[i]){g.set_hid_node(static_cast<int>(i),signal_hid);}
    }
}

void hid_updates_out(GRN& g) noexcept
{
  for(unsigned int i = 0; i != g.get_output_nodes().size(); i++)
    {
      double signal_out = 0;
      for (unsigned int j = 0 ; j != g.get_hidden_nodes().size(); j++)
        {
          signal_out  += g.get_hidden_nodes()[j] * g.get_H2O()[j][i];

        }
      if(signal_out > g.get_out_tresh()[i]){g.set_out_node(static_cast<int>(i),signal_out);}
    }
}

void jordi_response_mech(GRN& g) //Jordi style
{
  //First he updates output from hidden (is not the first thing he does
  //in the code but is the first thing that actually changes something in the
  //GRN
  hid_updates_out(g);
  //Update the state of the hidden nodes
  inp_updates_hid(g);
}


void test_GRN()//!OCLINT , tests may be long
{
  //A grn is initialized with connection weights
  //all to 0
  {
    GRN g;
    for(const auto& node : g.get_I2H())
      {
        for(const auto& weight : node)
          assert(weight < 0.00001 && weight > -0.00001);
      }
    for(const auto& node : g.get_H2H())
      {
        for(const auto& weight : node)
          assert(weight < 0.00001 && weight > -0.00001);
      }
    for(const auto& node : g.get_H2O())
      {
        for(const auto& weight : node)
          assert(weight < 0.00001 && weight > -0.00001);
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
      assert(hidden < 0.00001 && hidden > -0.00001);
    for(auto output : g.get_output_nodes())
      assert(output < 0.00001 && output > -0.00001);
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
    bool hidden_nodes_value = 1;
    //all hidden nodes are 1
    g.set_all_hid_nodes(hidden_nodes_value);
    //all weight H20 are
    //hidden_nodes_value / number_of_output_nodes / number_of_connections
    g.set_all_H2O(hidden_nodes_value /
              g.get_output_nodes().size() /
              (g.get_H2O().size() * g.get_H2O()[0].size()));
    //Therefore each output node will receive a signal of hidden_nodes_value
    g.set_all_out_tresh(hidden_nodes_value - 0.00001);
    hid_updates_out(g);
    for (size_t i = 0; i != g.get_output_nodes().size(); i++)
    assert(g.get_output_nodes()[i]);
  }

  //Given the initial inputs the hidden nodes can be updated
  {
    GRN g;
    double inputs_value = 3.14; //in this case all inputs are the same
    auto inputs = std::vector<double>(g.get_input_nodes().size(),inputs_value);
    g.set_inputs(inputs);
    g.set_all_hid_tresh(inputs_value - 0.00001);
    //all weight I2H are
    //inputs_value / number_of_output_nodes / number_of_connections
    g.set_all_I2H( inputs_value /
          g.get_hidden_nodes().size() /
             g.get_I2H().size() * g.get_I2H()[0].size());
    inp_updates_hid(g);
    //Therefore each output node will receive a signal of hidden_nodes_value
    for (size_t i = 0; i != g.get_hidden_nodes().size(); i++)
    assert(g.get_hidden_nodes()[i]);

  }
  //Given three initial input a GRN gives back an output
  {
    GRN g;//the weights and connection and nodes are all to 0, so output should be 0
    std::vector<double> inputs{100,100,100};
    g.set_inputs(inputs);
    jordi_response_mech(g);
    for(auto output : g.get_output_nodes())
      {
        assert(output < 0.00001
               && output > -0.00001);
      }
  }

}

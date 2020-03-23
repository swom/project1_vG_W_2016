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

void GRN::set_out_node(int index_node, bool state)
{m_ExOutput[static_cast<unsigned int>(index_node)] = state;}

void GRN::set_hid_node(int index_node, bool state)
{m_ExHidden[static_cast<unsigned int>(index_node)] = state;}

void jordi_response(GRN& g, std::vector<double> inputs) //Jordi style
{
  assert(inputs.size() == g.get_input_states().size());
  //First he updates output from hidden (is not the first thing he does
  //in the code but is the first thing that actually changes something in the
  //GRN
  for(unsigned int i = 0; i != g.get_output_states().size(); i++)
    {
      double signal_out = 0;
      for (unsigned int j = 0 ; j != g.get_hidden_states().size(); j++)
        {
          signal_out  += g.get_hidden_states()[j] * g.get_H2O()[j][i];

        }
      if(signal_out > g.get_out_tresh()[i]){g.set_out_node(static_cast<int>(i),signal_out);}
    }

  //Update the state of the hidden nodes
  for(unsigned int i = 0; i != g.get_hidden_states().size(); i++)
    {
      double signal_hid = 0;
      for (unsigned int j = 0 ; j != g.get_input_states().size(); j++)
        {
          signal_hid += inputs[i] * g.get_I2H()[j][i];
        }
      if(signal_hid > g.get_hid_tresh()[i]){g.set_hid_node(static_cast<int>(i),signal_hid);}
    }
}


void test_GRN()//!OClint, tests may be long
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
    for(auto input : g.get_input_states())
      assert(input < 0.00001 && input > -0.00001);
    for(auto hidden : g.get_hidden_states())
      assert(hidden < 0.00001 && hidden > -0.00001);
    for(auto output : g.get_output_states())
      assert(output < 0.00001 && output > -0.00001);
  }

  //Given three initial input a GRN gives back an output
  {
    GRN g;//the weights and connection and nodes are all to 0, so output should be 0
    std::vector<double> inputs{100,100,100};
    jordi_response(g, inputs);
    for(auto output : g.get_output_states())
      {
        assert(output < 0.00001
               && output > -0.00001);
      }
  }

}

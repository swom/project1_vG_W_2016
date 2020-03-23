#ifndef GRN_H
#define GRN_H
#include <vector>

class GRN
{
public:
  GRN(size_t n_input = 3, size_t n_inner = 3, size_t n_output = 2);

  ///Gets the connection from input to hidden layer
  const std::vector<std::vector<double> >& get_I2H() const noexcept {return m_ConI2H;}

  ///Gets the connection from hidden to hidden layer
  const std::vector<std::vector<double> >& get_H2H() const noexcept {return m_ConH2H;}

  ///Gets the connection from hidden to output layer
  const std::vector<std::vector<double> >& get_H2O() const noexcept {return m_ConH2O;}

  ///Gets const reference to input layer states
  const std::vector<double>& get_input_states() const noexcept {return  m_ExInput;}

  ///Gets const reference hidden layer states
  const std::vector<bool>& get_hidden_states() const noexcept {return  m_ExHidden;}

  ///Gets const ref to tresholds of hidden nodes
  const std::vector<double>& get_hid_tresh() const noexcept {return m_THidden;}

  ///Gets const reference output layer states
  const std::vector<bool>& get_output_states() const noexcept {return  m_ExOutput;}

  ///Gets reference output layer states
  std::vector<bool>& get_output_states() noexcept {return  m_ExOutput;}

  ///Gets const ref to tresholds of output nodes
  const std::vector<double>& get_out_tresh() const noexcept {return m_TOutput;}

  ///Sets an output node to true or false
  void set_hid_node(int index_node, bool state);

  ///Sets an output node to true or false
  void set_out_node(int index_node, bool state);


private:
  std::vector<std::vector<double> > m_ConI2H;	// Connections from Input to Hidden layer
  std::vector<std::vector<double> > m_ConH2H;	// Connections from Hidden to Hidden layer
  std::vector<std::vector<double> > m_ConH2O;	// Connections from Hidden to Output layer

  std::vector<double> m_THidden;			// Threshold hidden nodes
  std::vector<double> m_TOutput;			// Threshold output nodes

  std::vector<double> m_ExInput;		// Expression during 'development'
  std::vector<bool> m_ExHidden;			// Expression during 'development'
  std::vector<bool> m_ExOutput;			// Expression during 'development'
};

///The GRN reads inputs and gives back the outputs
///!!!!****Implemented as in Jordi's model***!!! I do not like it
void jordi_response(GRN& g, std::vector<double> inputs);

void test_GRN();
#endif // GRN_H

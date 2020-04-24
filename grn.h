#ifndef GRN_H
#define GRN_H
#include <vector>
#include <random>

class GRN
{
public:
  GRN(size_t n_input = 3, size_t n_hidden = 3, size_t n_output = 2, double weights = 0.1);

  ///Gets the constant reference to connection from input to hidden layer
  const std::vector<std::vector<double> >& get_I2H() const noexcept {return m_ConI2H;}

  ///Gets the reference to connection from input to hidden layer
  std::vector<std::vector<double> >& get_I2H() noexcept {return m_ConI2H;}

  ///Gets the const reference to connection from hidden to hidden layer
  const std::vector<std::vector<double> >& get_H2H() const noexcept {return m_ConH2H;}

  ///Gets the reference to connection from hidden to hidden layer
  std::vector<std::vector<double> >& get_H2H() noexcept {return m_ConH2H;}

  ///Gets the const ref to connection from hidden to output layer
  const std::vector<std::vector<double> >& get_H2O() const noexcept {return m_ConH2O;}

  ///Gets the reference connection from hidden to output layer
  std::vector<std::vector<double> >& get_H2O() noexcept {return m_ConH2O;}

  ///Gets const reference to input layer states
  const std::vector<double>& get_input_nodes() const noexcept {return  m_ExInput;}

  ///Gets const reference to input layer states
  std::vector<double>& get_input_nodes() noexcept {return  m_ExInput;}

  ///Gets const reference hidden layer states
  const std::vector<bool>& get_hidden_nodes() const noexcept {return  m_ExHidden;}

  ///Gets reference hidden layer states
  std::vector<bool>& get_hidden_nodes() noexcept {return  m_ExHidden;}

  ///Gets const ref to tresholds of hidden nodes
  const std::vector<double>& get_hid_tresh() const noexcept {return m_THidden;}

  ///Gets ref to tresholds of hidden nodes
  std::vector<double>& get_hid_tresh() noexcept {return m_THidden;}

  ///Gets const reference output layer states
  const std::vector<bool>& get_output_nodes() const noexcept {return  m_ExOutput;}

  ///Gets reference output layer states
  std::vector<bool>& get_output_nodes() noexcept {return  m_ExOutput;}

  ///Gets const reference to sporulation output (m_ExOutput[0])
  bool get_output_spo() const noexcept {return m_ExOutput[0];}

  ///Gets const ref to tresholds of output nodes
  const std::vector<double>& get_out_tresh() const noexcept {return m_TOutput;}

  ///Gets ref to tresholds of output nodes
  std::vector<double>& get_out_tresh() noexcept {return m_TOutput;}

  ///Sets all weights of H2O to a given value
  void set_all_H2H(double value) noexcept;

  ///Sets all weights of H2O to a given value
  void set_all_H2O(double value) noexcept;

  ///Sets all weights of I2H to a given value
  void set_all_I2H(double value) noexcept;

  ///Sets an output node to true or false
  void set_hid_node(int index_node, bool state);

  ///Sets an output node to true or false
  void set_out_node(int index_node, bool state);

  ///Sets all output nodes to a state: T/F
  void set_all_hid_nodes(bool state);

  ///Sets all output nodes to a state: T/F
  void set_all_out_nodes(bool state);

  ///Sets all the tresholds of the output nodes to a certain value
  void set_all_out_tresh(double value);

  ///Sets all the tresholds of the hidden nodes to a certain value
  void set_all_hid_tresh(double value);

  ///Sets inputs to the values given by a vector
  ///(that is created in individual in working circumstances)
  void set_inputs(std::vector<double> inputs);

//  ///Sets the values of the m_prev_hid_val vector
//  void set_prev_hid_val(std::vector<bool> hid_val);

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

///Input nodes update the states of output nodes
std::vector<double> hid_updates_hid(GRN& g) noexcept;

///Hidden nodes update the states of output nodes
void hid_updates_out(GRN& g) noexcept;

///Input nodes update the states of output nodes
std::vector<double> inp_updates_hid(GRN& g) noexcept;

///The GRN reads inputs and gives back the outputs
///!!!!****Implemented as in Jordi's model***!!! I do not like it
void jordi_response_mech(GRN& g);

///Mutates weights of connections between input and hidden layer
///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutation_I2H(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept;

///Mutates weights of connections between input and hidden layer
///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutation_H2H(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept;

///Mutates weights of connections between input and hidden layer
///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutation_H2O(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept;

///Mutates weights of thresholds of hidden layer
///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutation_hid_tr(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept;

///Mutates weights of thresholds of output layer
///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutation_hid_tr(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept;


///Applies random mutation to weights and thresholds of an individual
void mutation(GRN& g, std::minstd_rand& rng,
              std::bernoulli_distribution& mu_p,
              std::normal_distribution<double> mu_st) noexcept;

///Counts the number of connections in a network
int n_connections(const GRN& g) noexcept;

// ///Stores the values of hidden layer in m_prev_hid_val vector
// void store_hid_val(GRN& g) noexcept;

///Sums the weights of hidden to hidden connections
double sum_I2H(const GRN& g) noexcept;

///Sums the weights of hidden to hidden connections
double sum_H2H(const GRN& g) noexcept;

///Sums the weights of hidden to ouput connections
double sum_H2O(const GRN& g) noexcept;

///Updates hidden nodes states based on inut and hidden nodes themselves
void update_hid(GRN& g) noexcept;

///Calculates the sum of all weights
double weights_sum (const GRN& g) noexcept;

///Calculates the mean value of all weights
double weights_mean (const GRN& g) noexcept;

///Calculates the variance of weights in a GRN
double weights_var(const GRN& g) noexcept;



void test_GRN();
#endif // GRN_H

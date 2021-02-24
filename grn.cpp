#include "grn.h"
#include "../project1_vG_W_2016/utilities.h"
#include "algorithm"
#include <cassert>
#include <numeric>



GRN::GRN(size_t n_input, size_t n_hidden, size_t n_output, double weights):
    m_ConI2H(n_input,std::vector<double>(n_hidden, weights)),
    m_ConH2H(n_hidden,std::vector<double>(n_hidden, weights)),
    m_ConH2O(n_hidden,std::vector<double>(n_output, weights)),
    m_THidden(n_hidden,0),
    m_TOutput(n_output,0),
    m_ExInput(n_input,0),
    m_ExHidden(n_hidden,true),
    m_ExOutput(n_output,true)
{

}

GRN::GRN(std::vector<std::vector<double> > ConI2H,
         std::vector<std::vector<double> > ConH2H,
         std::vector<std::vector<double> > ConH2O,
         std::vector<double> THidden,
         std::vector<double> TOutput,
         std::vector<bool>ExHidden,
         int n_input,
         int n_output):
    m_ConI2H{ ConI2H},
    m_ConH2H{ ConH2H},
    m_ConH2O{ ConH2O},
    m_THidden{ THidden},
    m_TOutput{ TOutput},
    m_ExInput(n_input, 0),
    m_ExHidden{ExHidden},
    m_ExOutput(n_output, 0)
{

}


bool operator==(const GRN& lhs, const GRN& rhs)
{

    bool h2h = compare_weights_with_tolerance(lhs.get_H2H() , rhs.get_H2H());
    bool h2o = compare_weights_with_tolerance(lhs.get_H2O(), rhs.get_H2O());
    bool i2h = compare_weights_with_tolerance(lhs.get_I2H(), rhs.get_I2H());
    bool hid_t = compare_with_tolerance(lhs.get_hid_tresh(), rhs.get_hid_tresh());
    bool out_t = compare_with_tolerance(lhs.get_out_tresh(), rhs.get_out_tresh());
    bool hid_n = lhs.get_hidden_nodes() == rhs.get_hidden_nodes();
    //We do not look at input and output nodes since they might vary
    //depending on the status of the simulation
    //this might be true also for hidden nodes,
    //but not before the first time step of a cycle
    return h2h && h2o && i2h && hid_n && hid_t && out_t;
}

bool operator!=(const GRN& lhs, const GRN& rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const GRN& p)
{
    save_I2H_weights(os,p);
    save_H2H_weights(os,p);
    save_H2O_weights(os,p);
    save_THidden(os,p);
    save_TOutput(os,p);
    save_ExHidden(os,p);
    save_n_input_nodes(os, p);
    save_n_output_nodes(os, p);

    return os;
}

std::ifstream& operator>>(std::ifstream& is, GRN& grn)
{
    auto ConI2H = load_I2H_weights( is);
    auto ConH2H = load_H2H_weights( is);
    auto ConH2O = load_H2O_weights( is);
    auto THidden = load_hid_tresh( is);
    auto TOutput = load_out_tresh( is);
    auto ExHidden = load_hid_state( is);
    int n_input = load_n_input_nodes( is);
    int n_output = load_n_output_nodes( is);

    grn = GRN{ConI2H,
            ConH2H,
            ConH2O,
            THidden,
            TOutput,
            ExHidden,
            n_input,
            n_output};

    return is;
}

void GRN::set_hid_node(int index_node, bool state)
{m_ExHidden[static_cast<unsigned int>(index_node)] = state;}

void GRN::set_all_H2H(double value) noexcept
{
    for (auto& node : get_H2H())
    {
        for(auto& connection_weight : node)
            connection_weight = value;
    }
}

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

std::vector<std::vector<double>> calc_reaction_norm(const GRN& g,
                                                    double max_energy,
                                                    double max_food,
                                                    double max_metabolite,
                                                    double step,
                                                    int n_responses)
{
    std::vector<std::vector<double>> reaction_norm;
    for(int i = 0; i * step < max_energy; i++)
        for(int j = 0; j * step < max_food; j++)
            for(int z = 0; z * step < max_metabolite; z++)
            {
                auto grn = g;
                auto energy = i * step;
                auto food = j * step;
                auto metabolite = z * step;
                grn.set_inputs({energy,
                                food,
                                metabolite});
                ///Always run response one
                /// to actually see the response of the grn
                /// to the inputs,starting from the next response
                jordi_response_mech(grn);


                auto response = std::vector<double>{energy, food, metabolite};
                for(int r = 0; r != n_responses; r++)
                {
                    jordi_response_mech(grn);
                    response.push_back(static_cast<double>(grn.get_output_spo()));

                }
                reaction_norm.push_back(response);
            }
    return reaction_norm;
}

bool compare_weights_with_tolerance(const std::vector<std::vector<double>>& lhs,
                                    const std::vector<std::vector<double>>& rhs)
{
    return std::equal(lhs.begin(), lhs.end(), rhs.begin(), compare_with_tolerance);
}

std::string create_reaction_norm_name(int seed, int change_freq)
{
    return std::string{
        "reaction_norm_best_ind_s"
        + std::to_string(seed)
                + "_f"
                + std::to_string(change_freq)
                +".csv"
    };
}

GRN load_grn( const std::string& filename)
{
    GRN g;
    std::ifstream  i_f(filename);
    i_f >> g;
    return g;
}

std::vector<std::vector<double>> load_I2H_weights(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;
    std::string subline;
    std::vector<std::vector<double> > ConI2H;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');

    std::stringstream iss(line);
    while (std::getline(iss, subline ,'!').good())
    {
        std::stringstream isss(subline);
        std::vector<double> node_weights;
        double tmp;
        while(isss >> tmp)
        {
            node_weights.push_back(tmp);
            isss >> dummy;
        }
        ConI2H.push_back(node_weights);
        iss >> dummy; //Eat the " , " after the " ! "
    }
    is >> dummy; //Eat the " , " after the " | "

    return ConI2H;
}

std::vector<std::vector<double>> load_H2H_weights(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;
    std::string subline;
    std::vector<std::vector<double> > ConH2H;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');

    std::stringstream iss(line);
    while (std::getline(iss, subline ,'!').good())
    {
        std::stringstream isss(subline);
        std::vector<double> node_weights;
        double tmp;
        while(isss >> tmp)
        {
            node_weights.push_back(tmp);
            isss >> dummy;
        }
        ConH2H.push_back(node_weights);
        iss >> dummy; //Eat the " , " after the " ! "
    }
    is >> dummy; //Eat the " , " after the " | "

    return ConH2H;
}

std::vector<std::vector<double>> load_H2O_weights(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;
    std::string subline;
    std::vector<std::vector<double> > ConH2O;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');

    std::stringstream iss(line);
    while (std::getline(iss, subline ,'!').good())
    {
        std::stringstream isss(subline);
        std::vector<double> node_weights;
        double tmp;
        while(isss >> tmp)
        {
            node_weights.push_back(tmp);
            isss >> dummy;
        }
        ConH2O.push_back(node_weights);
        iss >> dummy; //Eat the " , " after the " ! "
    }
    is >> dummy; //Eat the " , " after the " | "

    return ConH2O;
}

std::vector<double> load_hid_tresh(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;
    std::vector<double> THidden;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');
    std::stringstream iss(line);
    double tmp;
    while(iss >> tmp)
    {
        THidden.push_back(tmp);
        iss >> dummy;
    }
    is >> dummy; //Eat the " , " after the " | "

    return THidden;
}

std::vector<bool> load_hid_state(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;
    std::vector<bool> ExHidden;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');
    std::stringstream iss(line);
    bool temporary;

    while(iss >> temporary)
    {
        ExHidden.push_back(temporary);
        iss >> dummy;
    }
    is >> dummy; //Eat the " , " after the " | "

    return ExHidden;
}


std::vector<double> load_out_tresh(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;
    std::vector<double> TOutput;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');
    std::stringstream iss(line);
    double tmp;
    while(iss >> tmp)
    {
        TOutput.push_back(tmp);
        iss >> dummy;
    }
    is >> dummy; //Eat the " , " after the " | "
    return TOutput;
}

int load_n_input_nodes(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');
    std::stringstream iss(line);
    int n_input_nodes;
    while(iss >> n_input_nodes)
    {
        iss >> dummy;
    }
    is >> dummy; //Eat the " , " after the " | "
    return n_input_nodes;
}

int load_n_output_nodes(std::ifstream& is)
{
    std::string dummy; // To remove the annotation in the file
    std::string line;

    //Get everything until character that signals the end of the layer
    std::getline(is, line, '|');
    std::stringstream iss(line);
    int n_output_nodes;
    while(iss >> n_output_nodes)
    {
        iss >> dummy;
    }
    is >> dummy; //Eat the " , " after the " | "
    return n_output_nodes;
}


std::vector<double> hid_updates_hid(GRN& g) noexcept
{
    std::vector<double> signal_from_hid;
    for(unsigned int i = 0; i != g.get_hidden_nodes().size(); i++)
    {
        double signal_hid = 0;
        for (unsigned int j = 0 ; j != g.get_hidden_nodes().size(); j++)
        {
            signal_hid += g.get_hidden_nodes()[j] * g.get_H2H()[j][i];
        }
        signal_from_hid.push_back(signal_hid);
    }
    return signal_from_hid;
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
        if(signal_out > g.get_out_tresh()[i]) {g.set_out_node(static_cast<int>(i),true);}
        else {g.set_out_node(static_cast<int>(i),false);}
    }
}

std::vector<double> inp_updates_hid(GRN& g) noexcept
{
    std::vector<double> signal_from_input;
    for(unsigned int i = 0; i != g.get_hidden_nodes().size(); i++)
    {
        double signal_hid = 0;
        for (unsigned int j = 0 ; j != g.get_input_nodes().size(); j++)
        {
            signal_hid += g.get_input_nodes()[j] * g.get_I2H()[j][i];
        }
        signal_from_input.push_back(signal_hid);
    }
    return signal_from_input;
}

void jordi_response_mech(GRN& g) //Jordi style
{
    //First he updates output from hidden (is not the first thing he does
    //in the code but is the first thing that actually changes something in the
    //GRN
    hid_updates_out(g);
    //Update the state of the hidden nodes
    update_hid(g);
}

std::vector<std::vector<double>> load_reaction_norm(std::string filename)
{
    std::ifstream f(filename);
    std::vector<std::vector<double>> reaction_norm;
    std::string line;
    std::string dummy;
    double value = 0.0;
    while (std::getline(f, line))
    {
        std::istringstream ss(line);
        std::vector<double> reaction;
        while (ss >> value)
        {
            reaction.push_back(value);
            ss >> dummy;
        }
        reaction_norm.push_back(reaction);
    }
    return reaction_norm;
}

void mutation_I2H(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept
{
    for(auto & node : g.get_I2H())
        for(auto & weight : node)
        {
            if(mu_p(rng)){weight += mu_st(rng);}
        }
}

void mutation_H2H(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept
{
    for(auto & node : g.get_H2H())
        for(auto & weight : node)
        {
            if(mu_p(rng)){weight += mu_st(rng);}
        }
}

void mutation_H2O(GRN& g, std::minstd_rand& rng,
                  std::bernoulli_distribution& mu_p,
                  std::normal_distribution<double> mu_st) noexcept
{
    for(auto & node : g.get_H2O())
        for(auto & weight : node)
        {
            if(mu_p(rng)){weight += mu_st(rng);}
        }
}

void mutation_hid_tr(GRN& g, std::minstd_rand& rng,
                     std::bernoulli_distribution& mu_p,
                     std::normal_distribution<double> mu_st) noexcept
{
    for(auto & treshold : g.get_hid_tresh())
    {
        if(mu_p(rng)){treshold += mu_st(rng);}
    }
}

void mutation_out_tr(GRN& g, std::minstd_rand& rng,
                     std::bernoulli_distribution& mu_p,
                     std::normal_distribution<double> mu_st) noexcept
{
    for(auto & treshold : g.get_out_tresh())
    {
        if(mu_p(rng)){treshold += mu_st(rng);}
    }
}

void mutation(GRN& g, std::minstd_rand& rng,
              std::bernoulli_distribution& mu_p,
              std::normal_distribution<double> mu_st) noexcept
{
    mutation_I2H(g,rng,mu_p,mu_st);
    mutation_H2H(g,rng,mu_p,mu_st);
    mutation_H2O(g,rng,mu_p,mu_st);
    mutation_hid_tr(g,rng,mu_p,mu_st);
    mutation_out_tr(g,rng,mu_p,mu_st);
}

int n_connections(const GRN& g) noexcept
{
    return static_cast<int> (
                g.get_H2H().size() * g.get_H2H()[0].size() +
            g.get_H2O().size() * g.get_H2O()[0].size() +
            g.get_I2H().size() * g.get_I2H()[0].size());
}

void save_grn(const GRN& grn, const std::string& filename)
{
    std::ofstream f(filename);
    f << grn;
}

std::ostream& save_I2H_weights(std::ostream& os, const GRN& p)
{
    for(const auto& node : p.get_I2H())
    {
        for(const auto& weight : node)
        {
            os << weight << " , ";
        }
        os << " ! , ";
    }
    os << " | , ";
    return os;
}

std::ostream& save_H2H_weights(std::ostream& os, const GRN& p)
{
    for(const auto& node : p.get_H2H())
    {
        for(const auto& weight : node)
        {
            os << weight << " , ";
        }
        os << " ! , ";
    }
    os << " | , ";
    return os;
}

std::ostream& save_H2O_weights(std::ostream& os, const GRN& p)
{
    for(const auto& node : p.get_H2O())
    {
        for(const auto& weight : node)
        {
            os << weight << " , ";
        }
        os << " ! , ";
    }
    os << " | , ";
    return os;
}

std::ostream& save_THidden(std::ostream& os, const GRN& p)
{
    for(const auto& treshold : p.get_hid_tresh())
    {
        os << treshold << " , ";
    }
    os << " | , ";
    return os;
}

std::ostream& save_TOutput(std::ostream& os, const GRN& p)
{
    for(const auto& treshold : p.get_out_tresh())
    {
        os << treshold << " , ";
    }
    os << " | , ";
    return os;
}

std::ostream& save_ExHidden(std::ostream& os, const GRN& p)
{
    int state;
    for(const auto& boolstate : p.get_hidden_nodes())
    {
        state = boolstate;
        os << state << " , ";
    }
    os << " | , ";
    return os;
}

std::ostream& save_n_input_nodes(std::ostream& os, const GRN& grn)
{
    int n_inputs = grn.get_input_nodes().size();
    os << n_inputs << " , "
       << " | " << " , ";
    return  os;
}

std::ostream& save_n_output_nodes(std::ostream& os, const GRN& grn)
{
    int n_outputs = grn.get_output_nodes().size();
    os << n_outputs << " , "
       << " | " << " , ";
    return os;
}

void save_reaction_norm(const std::vector<std::vector<double>> reaction_norm,
                        const std::string& filename)
{
    std::ofstream f(filename);
    for (const auto& reaction : reaction_norm)
    {
        for (const auto& value : reaction)
        {
            f << value << " , ";
        }
        f << std::endl;
    }
}

double sum_I2H(const GRN& g) noexcept
{
    double sum_I2H = 0;
    for (const auto& node : g.get_I2H()) {
        sum_I2H += std::accumulate(node.begin(), node.end(), 0.0);
    }
    return sum_I2H;
}

double sum_H2H(const GRN& g) noexcept
{
    double sum_H2H = 0;
    for (const auto& node : g.get_H2H()) {
        sum_H2H += std::accumulate(node.begin(), node.end(), 0.0);
    }
    return sum_H2H;
}

double sum_H2O(const GRN& g) noexcept
{
    double sum_H2O = 0;
    for (const auto& node : g.get_H2O()) {
        sum_H2O += std::accumulate(node.begin(), node.end(), 0.0);
    }
    return sum_H2O;
}

void update_hid(GRN& g) noexcept
{
    auto signal_in = inp_updates_hid(g);
    auto signal_hid = hid_updates_hid(g);
    assert(signal_in.size() == signal_hid.size());
    for(size_t i = 0; i != signal_in.size(); i++)
    {
        auto signal = signal_in[i] + signal_hid[i];
        if(signal > g.get_hid_tresh()[i])
        {g.set_hid_node(static_cast<int>(i),true);}
        else
        {g.set_hid_node(static_cast<int>(i),false);}
    }
}

double weights_sum (const GRN& g) noexcept
{
    return sum_H2H(g) + sum_H2O(g) + sum_I2H(g);
}

double weights_mean (const GRN& g) noexcept
{
    return weights_sum(g)/n_connections(g);
}

double weights_var(const GRN& g) noexcept
{
    double var = 0;
    double mean = weights_mean(g);
    for(const auto& node : g.get_H2H())
        for(const auto& weight : node)
        {
            var += (weight - mean) * (weight - mean);
        }
    for(const auto& node : g.get_H2O())
        for(const auto& weight : node)
        {
            var += (weight - mean) * (weight - mean);
        }
    for(const auto& node : g.get_I2H())
        for(const auto& weight : node)
        {
            var += (weight - mean) * (weight - mean);
        }
    return var /= n_connections(g);
}


#ifndef SIMULATION_H
#define SIMULATION_H
#include "vector"
#include"individual.h"

class simulation
{
public:
    simulation(int pop_size = 1, double min_dist = 0.1);

    ///finds the indexes of the individuals ready to divide
    std::vector<int> get_dividing_individuals() const noexcept;

    ///Gets an inidividual at a certain index in the vector
    const individual get_ind(int i) const
    {
      return population[static_cast<unsigned int>(i)];
    }

    ///Gets an inidividual at a certain index in the vector and allows to change it
    individual& get_ind(int i)
    {
      return population[static_cast<unsigned int>(i)];
    }

    ///Gets the excess energy from a referenced vector of index of individuals in population vector
    std::vector<double> get_excess_energies(std::vector<int> &v_ind) const noexcept;

    ///Gets the excess energy for a vector of index of individuals in population vector
    std::vector<double> get_excess_energies(std::vector<int> v_ind) const noexcept;

    ///Gets an individual's energy
    double get_ind_en(int i) const
    {
      return population[static_cast<unsigned int>(i)].get_energy();
    }

    ///Gets an individual's treshold energy
    double get_ind_tr_en(int i) const
    {
      return population[static_cast<unsigned int>(i)].get_treshold_energy();
    }

    ///Gets the position of an individual as a vector x,y
    const std::pair<double,double> get_ind_pos(int i);

    ///Gets distance in vector elements between two duaghters of same mother cell
    std::vector<std::pair<int, int> > get_sisters_index_offset() const noexcept;

    ///Get minimum distance between individuals at the start of the simulation
    double get_min_dist() const noexcept {return m_min_init_dist_btw_cells;}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return static_cast<int>(population.size());}

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_pop() const noexcept {return population;}

    ///Gets the vector containing all the individuals of the population
    std::vector<individual>& get_pop() noexcept {return population;}

    ///Divides an individual at a given index
    void divide_ind(individual &i);


    ///The individuals in the vector are copied at the end of the pop vector
    // To be chaged to take vec of references as an argument
    void do_division(std::vector<int> dividing_individuals);

    ///Individuals with enough energy divide and redistribute energy
    void do_reprduction() noexcept;

    /// Places cells in an hexagonal grid
    /// https://stackoverflow.com/questions/14280831/algorithm-to-generate-2d-magic-hexagon-lattice
    void place_start_cells() noexcept;

    ///Places an individual of index i at position x,y
    void set_ind_pos(individual& i, double x, double y);

    ///Places an individual of index i at position x,y
    void set_ind_pos(individual& i, const std::pair<double, double> pos);

    ///Sets and individual's energy
    void set_ind_en(int i, double en)
    {
      population[static_cast<unsigned int>(i)].set_energy(en);
    }


private:
    std::vector<individual> population;
    double m_min_init_dist_btw_cells;
};

///Counts the number of hexagonal layers necessary to place all individuals in hex pattern
int count_hex_layers(int pop_size)  noexcept ;

///Calculate the distance between two cells given their positions
double calculate_distance(std::pair<double, double> lhs, std::pair<double, double> rhs) noexcept;

///Checks if any two cells are colliding
bool has_collision(const simulation& s);

/// Displaces colliding cells so that they do not collide anymore
void manage_static_collisions(simulation& s);

void test_simulation();

#endif // SIMULATION_H

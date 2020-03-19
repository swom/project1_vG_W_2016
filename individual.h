#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include "vector"
#include <utility>
#include <algorithm>

class individual
{
public:

    individual(double x_pos = 0, double y_pos = 0, double size  = 1,
               double energy = 0, double treshold_energy = 2,
               double uptake_rate = 0.1);

     ///Changes x of an individual
     void change_x(double x) noexcept {m_x += x;}

     ///Changes y of an individual
     void change_y(double y) noexcept {m_y += y;}

     ///Changes the energy level of an individual
     void change_en(double en_change) noexcept {m_energy += en_change;}

     ///gets size of individual
     double get_size() const noexcept {return m_size;}

     ///gets energy of individual
     double get_energy() const noexcept {return m_energy;}

     ///gets the food uptake_rate of an individual
     double get_uptake_rate() const noexcept {return m_uptake_rate;}

     ///gets x position of the individual
     double get_x() const noexcept {return m_x;}

     ///gets y position of the individual
     double get_y() const noexcept {return m_y;}

     ///gets the treshold energy at which an individual reproduces
     double get_treshold_energy() const noexcept {return m_treshold_energy;}

     ///sets the energy of an individual
     void set_energy(double new_energy) {m_energy = new_energy;}

     ///Sets the x of an individual
     void set_x(double x) noexcept {m_x = x;}

     ///Sets the y of an individual
     void set_y(double y) noexcept {m_y = y;}

     ///Splits the excess energy not required for division in two9to be then
     ///(to be then assigned to the two daughter cells by
     /// simulation::reproduce/cells_divide)
     double split_excess_energy() const noexcept {return (m_energy - m_treshold_energy)/2;}

private:

     double m_x;
     double m_y;
     double m_size;
     double m_energy;
     double m_treshold_energy;
     double m_uptake_rate;

};

///returns the position of the second daughter cell just outside the mother
const std::pair<double,double> get_d2_pos(individual& i) noexcept;

///Gets the x,y coordinates as a pair
const std::pair<double,double> get_pos(individual& i) noexcept;

///Sets x,y given a pair of doubles(x,y)
void set_pos(individual& i, std::pair<double, double> pos)  noexcept;

/// Finds the distance between two individuals
double distance(const individual& lhs, const individual& rhs) noexcept;

///Finds the overlap between two individuals
double overlap(const individual& lhs, const individual& rhs) noexcept;

///Checks if two individuals are colliding
bool are_colliding(const individual &lhs, const individual &rhs) noexcept;

///Calculates how much individuals need to be displaced to not overlap
std::pair<double, double> get_displacement(const individual &lhs, const individual &rhs) noexcept;

///Dislpaces two individuals so that they do not overlap
void displace(individual& lhs, individual& rhs) noexcept;

///Finds the grid cell where an individual is on,
///given the side of the environment grid
int find_grid_index(const individual& i, double  grid_side);

///An individual increases its energy depleting food
void feed(individual& i, double& food) noexcept;

void test_individual();

#endif // INDIVIDUAL_H

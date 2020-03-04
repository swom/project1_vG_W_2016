#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H


class individual
{
public:

    individual(double x_pos, double y_pos, double size  = 0,
               double energy = 0, double treshold_energy = 2);

     ///gets size of individual
     double get_size() const noexcept {return m_size;}

     ///gets energy of individual
     double get_energy() const noexcept {return m_energy;}

     ///sets the energy of an individual
     void set_energy(double new_energy) {m_energy = new_energy;}

     ///gets x position of the individual
     double get_x() const noexcept {return m_x;}

     ///gets y position of the individual
     double get_y() const noexcept {return m_y;}

     ///gets the treshold energy at which an individual reproduces
     double get_treshold_energy() const noexcept {return m_treshold_energy;}

private:

     double m_x;
     double m_y;
     double m_size;
     double m_energy;
     double m_treshold_energy;

};

void test_individual();

#endif // INDIVIDUAL_H

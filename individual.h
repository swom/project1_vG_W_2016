#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H


class individual
{
public:

    individual(double x_pos, double y_pos, double size  = 0);

     ///gets size of individual
     double get_size() const noexcept {return m_size;}

     ///gets energy of individual
     double get_energy() const noexcept {return m_energy;}

     ///gets x position of the individual
     double get_x() const noexcept {return m_x;}

     ///gets y position of the individual
     double get_y() const noexcept {return m_y;}

private:

    double m_size;
    double m_energy = 0;
    double m_x;
    double m_y;
};

void test_individual();

#endif // INDIVIDUAL_H

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H


class individual
{
public:
    individual(double size);

    ///gets size of individual
     double get_size() const noexcept {return m_size;}

private:
    double m_size;
};

void test_individual();

#endif // INDIVIDUAL_H

#include "ind_param.h"
#include <cassert>

ind_param::ind_param(double radius,
                     double treshold_energy,
                     double uptake_rate,
                     double metabolic_rate,
                     int transformation_time,
                     double wiggle_room,
                     double metab_secretion_rate):
    m_metabolic_rate{metabolic_rate},
    m_metab_secr_rate{metab_secretion_rate},
    m_radius{radius},
    m_transformation_time{transformation_time},
    m_treshold_energy{treshold_energy},
    m_uptake_rate{uptake_rate},
    m_wiggle_room{wiggle_room}
{
    assert(m_metabolic_rate > -0.000001);
    assert(m_metab_secr_rate > -0.000001);
    assert(m_radius > -0.000001);
    assert(m_transformation_time > 0);
    assert(m_treshold_energy > -0.000001);
    assert(m_uptake_rate > -0.000001);
    assert(m_wiggle_room > -0.000001);
}

void test_ind_param() noexcept  //!OCLINT
{

}

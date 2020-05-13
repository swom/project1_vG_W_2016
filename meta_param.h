#ifndef META_PARAM_H
#define META_PARAM_H


class meta_param
{
public:
    meta_param(int n_cycles = 1, int cycle_duration = 2000);

    ///Returns number of cycles for which the simulation will last
    int get_n_cycles() const noexcept {return m_n_cycles;}

    ///Returns number of ticks per cycle
    int get_cycle_duration() const noexcept {return m_cycle_duration;}

private:

    ///The number of timesteps executed in one colony cycle
    int m_cycle_duration;

    ///The number of cycles of funding and growing
    /// of individual colonies that will be simulated
    int m_n_cycles;

};

void test_meta_param() noexcept;
#endif // META_PARAM_H

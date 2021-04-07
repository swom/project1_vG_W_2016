#ifndef LOGIC_ONLY
#include "sim_view.h"
#include <SFML/Graphics.hpp>
#include <cmath>
#include <cassert>


sim_view::sim_view(float start_zoom,
                   float zoom_step,
                   float pan_step,
                   float scale) :
    m_grid_view{scale},
    m_max_zoom{start_zoom},
    m_pan_step{pan_step},
    m_scale{scale},
    m_window{sf::VideoMode(1280, 720), "Sim_W_V_G_2016"},
    m_view{
        sf::Vector2f{0,0},
        sf::Vector2f(m_window.getSize().x / m_max_zoom * m_scale,
                     m_window.getSize().y / m_max_zoom * m_scale)
        },
    m_zoom_step{zoom_step}
{
#ifndef IS_ON_TRAVIS
    try {
        if(zoom_step > 1)
        {throw std::string{"zoom_step > 1... too high!\n"};}
    }
    catch (std::string e) {
        std::cout << e;
#ifdef NDEBUG
        abort();
#endif
    }
#endif
}

sim_view::~sim_view()
{

}

void sim_view::draw_background(const simulation& s) noexcept
{
    m_grid_view.update_grid_quads(s.get_env().get_grid());
    m_window.draw(m_grid_view);
}

void sim_view::draw_food() noexcept
{
    // Get position of food
    // Position in landscape
}


void sim_view::draw_inds(const simulation& s) noexcept
{
    update_pop(s);
    for(const auto& ind : m_pop_shapes)
    {
        m_window.draw(ind);
    }
}

void sim_view::exec(simulation& s) noexcept
{
    prepare(s);

    while(s.get_cycle() != s.get_meta_param().get_n_cycles() &&
          s.get_pop().get_pop_size() < s.get_meta_param().get_pop_max())
    {
        if(exec_cycle_visual(s))
            break;
        if(s.get_cycle() != 0
                && s.get_meta_param().get_change_freq() != 0
                && s.get_cycle() % s.get_meta_param().get_change_freq() == 0
                )
        {
            change_env(s);
            change_pop(s);
        }
        s.tick_cycles();
    }
    return;
}

void sim_view::exec_change(simulation& s,
                 const std::vector<std::pair<env_param, ind_param>>& rand_conditions)
{
    try {
        if(rand_conditions.size() > s.get_meta_param().get_n_cycles())
        {
            throw std::string{"More random conditions than cycles!"};
        }
        else if(s.get_meta_param().get_n_cycles() % rand_conditions.size())
        {
            throw std::string{"conditions are not a dividend of the number of cycles!"};
        }
    }
    catch (std::string& e) {
        std::cout << e << std::endl;
#ifdef NDEBUG
        abort();
#endif
    }

    size_t condition_index = 0;
    auto condition_duration = s.get_meta_param().get_n_cycles() / rand_conditions.size();

    while(s.get_cycle() != s.get_meta_param().get_n_cycles())
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto modulus =  s.get_cycle() % condition_duration;

        if( modulus == 0 && condition_index < rand_conditions.size())
        {
            reproduce_rand_cond(s, rand_conditions, condition_index);
            condition_index++;
        }

        exec_cycle_visual(s);

        s.reset_timesteps();
        s.tick_cycles();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout<< "cycle n " << s.get_cycle() << ":" << std::endl <<
                    "time: " << duration.count() << "s" << std::endl <<
                    "n_individuals: " << s.get_pop().get_pop_size() << std::endl << std::endl;
    }
}

bool sim_view::exec_cycle_visual(simulation& s) noexcept
{
    auto rand_start = std::chrono::high_resolution_clock::now();

    add_new_funders(s);
    while(s.get_timestep() != s.get_meta_param().get_cycle_duration() &&
          s.get_pop().get_pop_size() < s.get_meta_param().get_pop_max())
    {
        show(s);
        bool must_quit{process_events()};
        if (must_quit)
            return true;
        if(!m_stop)
        {
            tick(s, s.get_meta_param().get_collision_check_interval());
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random_cond_evo cycle n " << s.get_cycle() << ":" << std::endl <<
                "time: " << duration.count() << "s" << std::endl <<
                "n_individuals: " << s.get_pop().get_pop_size() << std::endl <<
                "n_timesteps: " << s.get_timestep() << std::endl <<
                std::endl;

    add_success_funders(s);
    store_demographics(s);
    dispersal(s);
    s.reset_timesteps();

    return false;
}

void sim_view::k_pan() noexcept
{
    m_view.move(0,m_pan_step * m_pan_up);
    m_view.move(0,- m_pan_step * m_pan_down);
    m_view.move(-m_pan_step * m_pan_left,0);
    m_view.move(m_pan_step * m_pan_right,0);
}

void sim_view::k_zoom() noexcept
{
    m_view.zoom( 1 - m_zoom_step * m_zoom_in);
    m_view.zoom(1 + m_zoom_step * m_zoom_out);
}

void sim_view::pan_k_input_starts(const sf::Event &event) noexcept
{
    switch (event.key.code) {
    case sf::Keyboard::Up:
        m_pan_up = true;
        break;
    case sf::Keyboard::Down:
        m_pan_down = true;
        break;
    case sf::Keyboard::Left:
        m_pan_left = true;
        break;
    case sf::Keyboard::Right:
        m_pan_right = true;
        break;
    default:
        break;
    }
}

void sim_view::pan_k_input_ends(const sf::Event &event) noexcept
{
    switch (event.key.code) {
    case sf::Keyboard::Up:
        m_pan_up = false;
        break;
    case sf::Keyboard::Down:
        m_pan_down = false;
        break;
    case sf::Keyboard::Left:
        m_pan_left = false;
        break;
    case sf::Keyboard::Right:
        m_pan_right = false;
        break;
    default:
        break;
    }
}

void sim_view::prepare (const simulation& s) noexcept
{
    prepare_pop(s);
    m_grid_view.prepare_grid(s.get_env().get_grid());
}

void sim_view::prepare_pop(const simulation& s) noexcept
{
    if(s.get_pop().get_v_ind().size() == 0)
    {
        auto r = individual{}.get_radius() * m_scale;
        // Create the sprite
        sf::CircleShape circle;
        circle.setRadius(r);
        circle.setFillColor(sf::Color::Blue);
        circle.setOutlineColor(sf::Color::Red);
        circle.setOutlineThickness(0.1f);
        circle.setOrigin(r, r);
        circle.setPosition(0, 0);
        circle.setPointCount(40);
        m_pop_shapes.push_back(circle);
        return;
    }

    for (const auto& ind :s.get_pop().get_v_ind())
    {
        // Type conversions that simplify notation
        const float r{static_cast<float>(ind.get_radius()) * m_scale};
        const float x{static_cast<float>(ind.get_x()) * m_scale};
        const float y{static_cast<float>(ind.get_y()) * m_scale};

        // Create the individual sprite
        sf::CircleShape circle;
        circle.setRadius(r);
        circle.setFillColor(sf::Color::Blue);
        circle.setOutlineColor(sf::Color::Red);
        circle.setOutlineThickness(0.1f);
        circle.setOrigin(r, r);
        circle.setPosition(x, y);
        circle.setPointCount(40);
        m_pop_shapes.push_back(circle);
    }
}

bool sim_view::process_events()
{

    // User interaction
    sf::Event event;
    while (m_window.pollEvent(event))
    {
        switch (event.type) {
        case sf::Event::Closed:
            m_window.close();
            return true; // Sim is done
        case sf::Event::KeyPressed:
            zoom_k_input_starts(event);
            pan_k_input_starts(event);
            show_food_input(event);
            show_metab_input(event);
            start_stop_input(event);
            break;
        case sf::Event::KeyReleased:
            zoom_k_input_ends(event);
            pan_k_input_ends(event);
            break;
        default:
            break;
        }
    }
    //Once all inputs are processed apply changes
    k_zoom();
    k_pan();
    return false; // if no events proceed with tick
}

void sim_view::show(const simulation& s) noexcept
{
    // Start drawing the new frame, by clearing the screen
    m_window.clear();
    m_window.setView(m_view);

    draw_background(s);

    draw_inds(s);

    // Display all shapes
    m_window.display();
}

void sim_view::update_pop(const simulation& s) noexcept
{
    if(m_pop_shapes.size() != s.get_pop().get_v_ind().size())
    {
        m_pop_shapes.resize(s.get_pop().get_v_ind().size(), m_pop_shapes[0]);
        for( size_t i = 0; i != m_pop_shapes.size(); i++)
        {
            float x = m_scale * static_cast<float>(s.get_pop().get_v_ind()[i].get_x());
            float y = m_scale * static_cast<float>(s.get_pop().get_v_ind()[i].get_y());
            m_pop_shapes[i].setPosition(x, y);
        }
    }

    for( size_t i = 0; i != m_pop_shapes.size(); i++)
    {
        float radius = m_scale * static_cast<float>(s.get_pop().get_v_ind()[i].get_radius());
        m_pop_shapes[i].setRadius(radius);
    }

    for( size_t i = 0; i != m_pop_shapes.size(); i++)
    {
        if(s.get_pop().get_v_ind()[i].get_phen() == phenotype::active)
        {
            m_pop_shapes[i].setFillColor(sf::Color::Blue);
            continue;
        }
        else if(s.get_pop().get_v_ind()[i].get_phen() == phenotype::sporulating)
        {
            m_pop_shapes[i].setFillColor(sf::Color::Magenta);
            continue;
        }
        else
            m_pop_shapes[i].setFillColor(sf::Color::Yellow);

    }
}


void sim_view::update_indexes() noexcept
{
    if(m_pop_shapes.size() != m_pop_indexes.size())
    {
        m_pop_indexes.resize(m_pop_shapes.size(), m_pop_indexes[0]);
        for( size_t i = 0; i != m_pop_shapes.size(); i++)
        {
            m_pop_indexes[i].setPosition(m_pop_shapes[i].getPosition());
            m_pop_indexes[i].setString(std::to_string(i));
        }
    }
}

void sim_view::show_food_input(sf::Event& event)
{
    if(event.key.code == sf::Keyboard::F)
    {
        m_grid_view.show_food_conc();
    }
}

void sim_view::show_metab_input(sf::Event& event)
{
    if(event.key.code == sf::Keyboard::M)
    {
        m_grid_view.show_metab_conc();
    }
}
void sim_view::start_stop_input(const sf::Event& event) noexcept
{
    if (event.key.code == sf::Keyboard::S)
    {
        if(m_stop)
            m_stop = false;
        else
            m_stop = true;
    }
}


void sim_view::zoom_k_input_ends(const sf::Event &event) noexcept
{
    if (event.key.code == sf::Keyboard::Add)
    {
        m_zoom_in = false;
    }
    else if(event.key.code == sf::Keyboard::Subtract)
    {
        m_zoom_out = false;
    }
}

void sim_view::zoom_k_input_starts(const sf::Event &event) noexcept
{
    if (event.key.code == sf::Keyboard::Add)
    {
        m_zoom_in = true;
    }
    else if(event.key.code == sf::Keyboard::Subtract)
    {
        m_zoom_out = true;
    }
}

#ifndef LOGIC_ONLY
int run_visual_evo (const env_param& e,
                    const ind_param& i,
                    const meta_param& m,
                    const pop_param& p)
{

    simulation s{sim_param{e, i, m, p}};
    sim_view v;
    v.exec(s);
    save_data(s);
    return 0;
}

int replay_cycle_from_evo (
        int change_frequency,
        int seed,
        int cycle)
{

    simulation s = load_sim_no_pop(seed, change_frequency);
    auto f_el = std::find_if(s.get_funders_success().get_v_funders().begin(),
                             s.get_funders_success().get_v_funders().end(),
                             [& cycle](const funders& f){return f.get_cycle() == cycle;});
    auto last_pop = *f_el;
    s.get_pop().get_v_ind().resize(last_pop.get_v_funder_data().size());
    for(size_t i = 0; i != last_pop.get_v_funder_data().size(); i++)
    {
        s.get_pop().get_v_ind()[i].get_grn()
                = last_pop.get_v_funder_data()[i].get_grn();
    }

    sim_view v;
    v.prepare(s);
    v.exec_cycle_visual(s);
    return 0;
}

int  replay_rand_cond_evo (double change_freq,
                           int seed_sim,
                           int n_conditions,
                           double amplitude,
                           int seed_rand_cond,
                           int rand_cond_n,
                           int pop_max)
{
    auto rand_cond = load_random_conditions(create_name_vec_rand_cond(n_conditions, amplitude, seed_rand_cond));
    auto rand_s = load_sim_last_pop(seed_sim,change_freq);

    rand_s.get_pop().get_param() = pop_param{1,100,0.01,0.0015,0.1,0.01,10};
    rand_s.get_meta_param().get_pop_max() = pop_max;

    sim_view v;
    reproduce_rand_cond(rand_s,rand_cond, rand_cond_n);

    //just change uptake rate for testing
    auto uptake_rate = 0.115;
    auto i = rand_s.get_pop().get_v_ind()[0].get_param();
    i.set_uptake_rate(uptake_rate);
    rand_s.get_pop().set_pop_inds(set_new_ind_par(rand_s.get_pop().get_v_ind(),i));

    std::cout << "uptake rate: " << uptake_rate << std::endl;

    ///keep this for repeatability
    ///(it is useless but if taken away will change rng behaviour in replicates)
    ///I am an idiot
    {
        std::vector<double> rnd_n;
        for(int i = 0; i !=5; i++)
            rnd_n.push_back( rand_s.get_rng()());
    }

    std::cout<< "base_fitness: " << rand_s.get_pop().get_param().get_base_disp_prob()
             <<" spore advantage: " << rand_s.get_pop().get_param().get_spo_adv() << std::endl;

    auto rand_start = std::chrono::high_resolution_clock::now();

    v.exec(rand_s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition n" << rand_cond_n <<": " << duration.count() << "s" << std::endl;

    save_demographic_sim(rand_s.get_demo_sim(),"test_death_selection" + create_sim_demo_name(rand_s));
    auto tot_inds = rand_s.get_demo_sim().get_demo_cycles().back().get_n_actives() +
            rand_s.get_demo_sim().get_demo_cycles().back().get_n_spores() +
            rand_s.get_demo_sim().get_demo_cycles().back().get_n_sporulating();
    std::cout<< "n individuals:" << tot_inds << std::endl;

    return 0;
}
int  replay_rand_cond_test (double change_freq,
                            int seed_sim,
                            int n_conditions,
                            double amplitude,
                            int seed_rand_cond,
                            int rand_cond_n,
                            int pop_max)
{
    auto rand_cond = load_random_conditions(create_name_vec_rand_cond(n_conditions, amplitude, seed_rand_cond));
    auto rand_s = load_sim_last_pop(seed_sim,change_freq);
    rand_s.get_meta_param().get_pop_max() = pop_max;

    sim_view v;
    reproduce_rand_cond(rand_s,rand_cond, rand_cond_n);
    v.prepare(rand_s);

    auto rand_start = std::chrono::high_resolution_clock::now();

    v.exec_cycle_visual(rand_s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition n" << rand_cond_n <<": " << duration.count() << "s" << std::endl;

    auto tot_inds = rand_s.get_demo_sim().get_demo_cycles().back().get_n_actives() +
            rand_s.get_demo_sim().get_demo_cycles().back().get_n_spores() +
            rand_s.get_demo_sim().get_demo_cycles().back().get_n_sporulating();
    std::cout<< "n individuals:" << tot_inds << std::endl;

    return 0;
}

int  replay_best_rand_cond (double change_freq,
                            int seed_sim,
                            int n_conditions,
                            double amplitude,
                            int seed_rand_cond,
                            int rand_cond_n,
                            int pop_max)
{
    auto rand_cond = load_random_conditions(create_name_vec_rand_cond(n_conditions, amplitude, seed_rand_cond));
    auto rand_s = load_best_ind_for_rand_cond(seed_sim,change_freq);
    rand_s.get_meta_param().get_pop_max() = pop_max;
    if(exists(create_best_random_condition_name(amplitude, change_freq, seed_sim)))
    {
        auto results = load_demographic_sim(create_best_random_condition_name(amplitude, change_freq, seed_sim));
        assert( results.get_demo_cycles()[rand_cond_n].get_env_param() == rand_cond[rand_cond_n].first );
        assert( results.get_demo_cycles()[rand_cond_n].get_ind_param() == rand_cond[rand_cond_n].second );

    }
    else {
        std::cout << "best_rand file yet not available,"
                     " cannot compare condition vector"
                     " with condition in results" << std::endl;
    }

    sim_view v;
    reproduce_rand_cond(rand_s,rand_cond, rand_cond_n);
    v.prepare(rand_s);

    auto rand_start = std::chrono::high_resolution_clock::now();
    v.exec_cycle_visual(rand_s);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition n" << rand_cond_n <<": " << duration.count() << "s" << std::endl;

    auto tot_inds = rand_s.get_demo_sim().get_demo_cycles().back().get_n_actives() +
            rand_s.get_demo_sim().get_demo_cycles().back().get_n_spores() +
            rand_s.get_demo_sim().get_demo_cycles().back().get_n_sporulating();
    std::cout<< "n individuals:" << tot_inds << std::endl;

    return 0;
}
#endif

funders_success run_evo_random_conditions_visual(const simulation& s,
                                          int number_of_sequences,
                                          int cond_per_seq,
                                          int seq_index,
                                          int pop_max,
                                          double amplitude,
                                          std::string prefix)
{
    auto random_conditions = create_rand_conditions_matrix_extreme(
                s.get_env().get_param(),
                ind_param{},
                number_of_sequences,
                cond_per_seq,
                amplitude);

    simulation rand_s = no_dem_and_fund_copy(s);

    ///For now do not allow the env and ind param
    /// to change, no matter the freq of change.
    rand_s.get_meta_param().get_change_freq() = 0;
    rand_s.get_meta_param().get_pop_max() = pop_max;

    sim_view v;
    v.prepare(rand_s);

    auto start = std::chrono::high_resolution_clock::now();

   v.exec_change(rand_s, random_conditions[seq_index]);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout<< "duration : " << duration.count() << "s" << std::endl <<
                "n_individuals: " << s.get_pop().get_pop_size() << std::endl <<
                std::endl;

    std::cout << "saving demographics" << std::endl;
    save_demographic_sim(rand_s.get_demo_sim(), create_sim_demo_name(s, "death" + prefix));
    std::cout << "saving funders" << std::endl;
    save_funders_success(rand_s.get_funders_success(), create_funders_success_name(rand_s, "death" + prefix));
    save_sim_parameters(rand_s.get_sim_param(), create_sim_par_name(rand_s, "death" + prefix));

    return rand_s.get_funders_success();
}

int run_sim_evo_rand_visual(double amplitude,
                     int change_frequency,
                     int num_of_sequences,
                     int conditions_per_seq,
                     int pop_max,
                     int seed,
                     int seq_index,
                     bool overwrite)
{

    auto s = load_sim_last_pop(seed,change_frequency);

    auto prefix = "rand_evo_extreme_a" + std::to_string(amplitude) +
            "seq_" + std::to_string(seq_index) +
            "cond_per_seq" + std::to_string(conditions_per_seq);

    if(exists(prefix + create_sim_demo_name(s)) && !overwrite)
    {
        std::cout << "The evolution in random conditions for this simulation"
                     " have already been tested!" << std::endl;
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();
    run_evo_random_conditions_visual(s,
                              num_of_sequences,
                              conditions_per_seq,
                              seq_index,
                              pop_max,
                              amplitude,
                              prefix);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout<< "random condition evo :" << duration.count() << "s" << std::endl;
    return 0;
}

void test_sim_view()//!OCLINT tests may be many
{
#ifndef IS_ON_TRAVIS

    {
        // Show the sim for one frame
        // (there will be a member function 'exec' for running the game)
        simulation s{sim_param{20}};
        sim_view v;
        v.prepare(s);
        v.show(s);
    }


    //A sim_view has a member variable that is a sf::View
    //Initialized with center on 0 and
    //with showing a 100th of the window dimensions(in pixel)
    {
        const sim_view v;

        assert((v.get_view().getCenter() == sf::Vector2f{0.f,0.f}));

        assert(v.get_view().getSize().x -
               static_cast<float>(v.get_window().getSize().x) /
               v.get_max_zoom() * v.get_scale() < 0.00001f
               &&
               v.get_view().getSize().x -
               static_cast<float>(v.get_window().getSize().x) /
               v.get_max_zoom() * v.get_scale() > -0.00001f);

        assert(v.get_view().getSize().y -
               static_cast<float>(v.get_window().getSize().y) /
               v.get_max_zoom() * v.get_scale() < 0.00001f &&
               v.get_view().getSize().y -
               static_cast<float>(v.get_window().getSize().y) /
               v.get_max_zoom() * v.get_scale() > -0.00001f);
    }


    //A sim_view is initialized with a m_zoom_step_variable,
    //that will determine how fast you will zoom in and out
    //0.1 by default,(increase/decrease size of view by 1/5th)
    //An error will be thrown if the zoom_step is more than 1!
    {
        float zoom_step = 0.314f;
        sim_view v(10,zoom_step);
        assert(v.get_zoom_step() - zoom_step < 0.0001f &&
               v.get_zoom_step() - zoom_step > -0.0001f);
        zoom_step = 5;
        try {
            sim_view v1(10, zoom_step);//!OCLINT
        }
        catch (std::string e) {
            assert(e == "zoom_step > 1... too high!\n" );
        }
    }

    //A sim_view is initialized with 2 boolean variables zoom_in and zoom_out
    //Both initialized at false
    {
        sim_view v;
        assert(v.get_zoom_in() + v.get_zoom_out() == 0);
    }

    //It is possible to zoom in or out pressing the keys + or  -
    //----------> tested graphically

    //A sim_view is initialized with a grid_view object containing a vector of vertices
    //of size == to m_sim.get_env().get_grid_size() * 4
    {
        simulation s;
        sim_view v;
        v.get_grid_view().prepare_grid(s.get_env().get_grid());
        assert(v.get_grid_view().get_grid_vert_size() == s.get_env().get_grid_size() * 4);
    }

    //When the F key is pressed the flag signalling that food should be shown on the grid will
    //turn true and the flag signalling that metabolite should not be shown will be turned to
    //false tested graphically


#endif
}
#endif

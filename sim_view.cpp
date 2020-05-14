#include "sim_view.h"
#include <SFML/Graphics.hpp>
#include <cmath>
#include <cassert>


sim_view::sim_view(sim_param start_simulation_param,
                   float start_zoom,
                   float zoom_step,
                   float pan_step,
                   float scale) :
    m_sim(start_simulation_param),//!!!Initialize this first!!!
    m_grid_view{scale},
    m_max_zoom{start_zoom},
    m_pan_step{pan_step},
    m_scale{scale},
    m_window{sf::VideoMode(1280, 720), "Sim_W_V_G_2016"},
    m_view{
        sf::Vector2f{0,0},
        sf::Vector2f(m_window.getSize().x / m_max_zoom,
                     m_window.getSize().y / m_max_zoom)
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

    m_grid_view.prepare_grid(m_sim.get_env().get_grid());
    prepare_pop();
}

sim_view::~sim_view()
{

}

void sim_view::draw_background() noexcept
{
    m_grid_view.update_grid_quads(m_sim.get_env().get_grid());
    m_window.draw(m_grid_view);
}

void sim_view::draw_food() noexcept
{
    // Get position of food
    // Position in landscape
}


void sim_view::draw_inds() noexcept
{
    update_pop();
    for(const auto& ind : m_pop_shapes)
    {
        m_window.draw(ind);
    }
}

void sim_view::exec() noexcept
{

    while(m_sim.get_cycle() != m_sim.get_meta_param().get_n_cycles())
    {
        if(exec_cycle_visual(m_sim))
            break;
        dispersal(m_sim);
        m_sim.tick_cycles();
    }
    m_window.close();
    return;
}

bool sim_view::exec_cycle_visual(simulation& s) noexcept
{
    while(s.get_timestep() != s.get_meta_param().get_cycle_duration())
    {
        bool must_quit{process_events()};
        if (must_quit)
            return true;
        if(!m_stop)
            tick(m_sim);
        show();
    }
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

void sim_view::prepare_pop() noexcept
{
    for (size_t i = 0 ; i != m_sim.get_pop().get_v_ind().size(); i++)
    {
        const auto& ind = m_sim.get_pop().get_v_ind()[i];
        // Type conversions that simplify notation
        const float r{static_cast<float>(ind.get_param().get_radius()) * m_scale};
        const float x{static_cast<float>(ind.get_x()) * m_scale};
        const float y{static_cast<float>(ind.get_y()) * m_scale};

        // Create the player sprite
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

void sim_view::show() noexcept
{
    // Start drawing the new frame, by clearing the screen
    m_window.clear();
    m_window.setView(m_view);

    draw_background();

    draw_inds();

    // Display all shapes
    m_window.display();
}

void sim_view::update_pop() noexcept
{
    if(m_pop_shapes.size() != m_sim.get_pop().get_v_ind().size())
    {
        m_pop_shapes.resize(m_sim.get_pop().get_v_ind().size(), m_pop_shapes[0]);
        for( size_t i = 0; i != m_pop_shapes.size(); i++)
        {
            float x = m_scale * static_cast<float>(m_sim.get_pop().get_v_ind()[i].get_x());
            float y = m_scale * static_cast<float>(m_sim.get_pop().get_v_ind()[i].get_y());
            m_pop_shapes[i].setPosition(x, y);
        }
    }

    for( size_t i = 0; i != m_pop_shapes.size(); i++)
    {
        if(m_sim.get_pop().get_v_ind()[i].get_phen() == phenotype::active)
        {
            m_pop_shapes[i].setFillColor(sf::Color::Blue);
            continue;
        }
        if(m_sim.get_pop().get_v_ind()[i].get_phen() == phenotype::sporulating)
        {
            m_pop_shapes[i].setFillColor(sf::Color::Magenta);
            continue;
        }
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

void test_sim_view()//!OCLINT tests may be many
{
#ifndef IS_ON_TRAVIS

    {
        // Show the game for one frame
        // (there will be a member function 'exec' for running the game)
        sim_view v(sim_param{20});
        v.show();
    }

    //A window starts showing simulation from tick 0
    {
        const sim_view v;
        assert(v.get_sim().get_timestep() == 0);
    }

    //A sim_view has a member variable that is a sf::View
    //Initialized with center on 0 and
    //with showing a 100th of the window dimensions(in pixel)
    {
        const sim_view v;
        assert((v.get_view().getCenter() == sf::Vector2f{0.f,0.f}));
        assert(v.get_view().getSize().x -
               static_cast<float>(v.get_window().getSize().x) / v.get_max_zoom() < 0.00001f &&
               v.get_view().getSize().x -
               static_cast<float>(v.get_window().getSize().x) / v.get_max_zoom() > -0.00001f);

        assert(v.get_view().getSize().y -
               static_cast<float>(v.get_window().getSize().y) / v.get_max_zoom() < 0.00001f &&
               v.get_view().getSize().y -
               static_cast<float>(v.get_window().getSize().y) / v.get_max_zoom() > -0.00001f);
    }


    //A sim_view is initialized with a m_zoom_step_variable,
    //that will determine how fast you will zoom in and out
    //0.1 by default,(increase/decrease size of view by 1/5th)
    //An error will be thrown if the zoom_step is more than 1!
    {
        float zoom_step = 0.314f;
        sim_view v(sim_param() ,10,zoom_step);
        assert(v.get_zoom_step() - zoom_step < 0.0001f &&
               v.get_zoom_step() - zoom_step > -0.0001f);
        zoom_step = 5;
        try {
            sim_view v1(sim_param(), 10, zoom_step);//!OCLINT
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
        sim_view v;
        assert(v.get_grid_view().get_grid_vert_size() == v.get_sim().get_env().get_grid_size() * 4);
    }

    //When the F key is pressed the flag signalling that food should be shown on the grid will
    //turn true and the flag signalling that metabolite should not be shown will be turned to
    //false tested graphically


#endif
}

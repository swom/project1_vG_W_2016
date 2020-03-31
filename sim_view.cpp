#include "sim_view.h"
#include <SFML/Graphics.hpp>
#include <cmath>
#include <cassert>


sim_view::sim_view(simulation start_simulation, float start_zoom, float zoom_step, double scale) :
  m_sim(start_simulation),
  m_window{sf::VideoMode(1280, 720), "Sim_W_V_G_2016"},
  m_max_zoom{start_zoom},
  m_scale{scale},
  m_view{
    sf::Vector2f{0,0},
    sf::Vector2f(m_window.getSize().x / m_max_zoom,
                 m_window.getSize().y / m_max_zoom)
    },
  m_zoom_step{zoom_step}
{
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
}

sim_view::~sim_view()
{

}

void sim_view::exec() noexcept
{
  while (m_window.isOpen())
    {
      bool must_quit{process_events()};
      if (must_quit)
        return;
      tick(m_sim);
      show();
    }
}

void sim_view::draw_background() noexcept
{
  // Draw the background
}

void sim_view::draw_food() noexcept
{
  // Get position of food
  // Position in landscape
}


void sim_view::draw_inds() noexcept
{
  for (const auto &ind : m_sim.get_pop())
    {
      // Type conversions that simplify notation
      const float r{static_cast<float>(ind.get_radius() * m_scale)};
      const float x{static_cast<float>(ind.get_x() * m_scale)};
      const float y{static_cast<float>(ind.get_y() * m_scale)};

      // Create the player sprite
      sf::CircleShape circle;
      circle.setRadius(r);
      circle.setFillColor(sf::Color::Blue);
      circle.setOutlineColor(sf::Color::Black);
      circle.setOutlineThickness(2.f);
      circle.setOrigin(r, r);
      circle.setPosition(x, y);

      // Draw the player
      m_window.draw(circle);
    }
}

bool sim_view::process_events()
{

  // User interaction
  sf::Event event;
  while (m_window.pollEvent(event))
    {
      if (event.type == sf::Event::Closed)
        {
          m_window.close();
          return true; // Sim is done
        }
      zoom_input(event);
    }
  return false; // if no events proceed with tick
}

void sim_view::show() noexcept
{
  // Start drawing the new frame, by clearing the screen
  m_window.clear();

  if(m_zoom_dir < 0 || m_zoom_dir >0)
    m_view.zoom(1 + m_zoom_step * m_zoom_dir);
  m_window.setView(m_view);

  draw_background();

  draw_inds();

  draw_food();

  // Display all shapes
  m_window.display();
}

void sim_view::set_zoom_dir(float dir) noexcept
{
  if(dir > 0)
    {
      m_zoom_dir = 1;
    }
  else if(dir < 0)
    {
      m_zoom_dir = -1;
    }
  else
    {
      m_zoom_dir = 0;
    }
}

void sim_view::zoom_input(const sf::Event& event) noexcept
{
  if(event.type == sf::Event::MouseWheelScrolled)
    {

      set_zoom_dir(event.mouseWheel.delta);
    }
  else if(event.key.code == sf::Keyboard::Add)
    {
      set_zoom_dir(1);
    }
  else if(event.key.code == sf::Keyboard::Hyphen)
    {
      set_zoom_dir(-1);
    }
  else
    {
      set_zoom_dir(0);
    }
}
void test_sim_view()
{
#ifndef NDEBUG

  {
    // Show the game for one frame
    // (there will be a member function 'exec' for running the game)
    simulation s(20);
    sim_view v(s);
    v.show();
  }

  //A window starts showing simulation from tick 0
  {
    const sim_view v;
    assert(v.get_sim().get_tick() == 0);
  }

  //A sim_view has a member variable that is a sf::View
  //Initialized with center on 0 and
  //with showing a 100th of the window dimensions(in pixel)
  {
    const sim_view v;
    assert((v.get_view().getCenter() == sf::Vector2f{0.f,0.f}));
    assert(v.get_view().getSize().x - static_cast<float>(v.get_window().getSize().x) / v.get_max_zoom() < 0.00001f &&
           v.get_view().getSize().x - static_cast<float>(v.get_window().getSize().x) / v.get_max_zoom() > -0.00001f);

    assert(v.get_view().getSize().y - static_cast<float>(v.get_window().getSize().y) / v.get_max_zoom() < 0.00001f &&
           v.get_view().getSize().y - static_cast<float>(v.get_window().getSize().y) / v.get_max_zoom() > -0.00001f);
  }
  //A sim_view is initialized with a m_zoom_step_variable, that will determine how fast you will zoom in and out
  //0.2 by default,(increase/decrease size of view by 1/5th)
  //An error will be thrown if the zoom_step is more than 1!
  {
    float zoom_step = 0.314f;
    sim_view v;
    assert(v.get_zoom_step() - zoom_step < 0.0001f &&
           v.get_zoom_step() - zoom_step > -0.0001f);
    zoom_step = 5;
    try {
      sim_view v1(simulation(),10,zoom_step);
    } catch (std::string e) {
      assert(e == "zoom_step > 1... too high!\n" );
    }
  }

  //It is possible to zoom out and in
  //By default the zoom starts from 0
  {
    sim_view v;
    auto init_view_size = v.get_view().getSize();
    v.set_zoom_dir(1);
    v.show();
    assert(v.get_view().getSize() != init_view_size);
  }
#endif
}

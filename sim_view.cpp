#include "sim_view.h"
#include <SFML/Graphics.hpp>
#include <cmath>
#include <cassert>


sim_view::sim_view(simulation start_simulation) :
  m_sim(start_simulation),
  m_window(sf::VideoMode(1280, 720), "Sim_W_V_G_2016"),
  m_view(sf::Vector2f{0,0},
         sf::Vector2f(m_window.getSize().x,m_window.getSize().y))
  {

}

sim_view::~sim_view()
{

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
    }
  return false; // if no events proceed with tick
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


void sim_view::draw_inds() noexcept //!OCLINT too long indeed, please
//! shorten
{
  for (const auto &ind : m_sim.get_pop())
    {
      // Type conversions that simplify notation
      const float r{static_cast<float>(ind.get_radius()) * 100};
      const float x{static_cast<float>(ind.get_x() * 100)};
      const float y{static_cast<float>(ind.get_y() * 100)};

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

void sim_view::show() noexcept
{
  // Start drawing the new frame, by clearing the screen
  m_window.clear();

  m_window.setView(m_view);

  draw_background();

  draw_inds();

  draw_food();

  // Display all shapes
  m_window.display();
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

#endif
}

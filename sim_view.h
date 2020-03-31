#ifndef SIM_VIEW_H
#define SIM_VIEW_H
#include "simulation.h"
#include <SFML/Graphics.hpp>


/// The game's main window
/// Displays the game class
class sim_view
{
public:

  sim_view(simulation start_simulation = simulation());
  ~sim_view();

  /// Show one frame
  void show() noexcept;

  /// Run the game until the window is closed
  void exec() noexcept;

  /// Processes events in game and ouputs false if quit
  /// is inputted
  bool process_events();

  ///Gets a const ref to m_sim
  const simulation& get_sim() const noexcept {return m_sim; }

  ///Gets constant ref to sf::RenderWindow m_window
  const sf::RenderWindow& get_window() const noexcept {return m_window; }

private:
  /// The simulation logic
  simulation m_sim;

  /// The window to draw to
  sf::RenderWindow m_window;

  ///The view of the simulation
  sf::View m_view;

  /// Draws players
  void draw_inds() noexcept;

  ///Draws the background
  void draw_background() noexcept;

  ///Draws food
  void draw_food() noexcept;

  /// Draws projectiles
  void draw_metabolite() noexcept;

};

void test_sim_view();

#endif // SIM_VIEW_H

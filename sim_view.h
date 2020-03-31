#ifndef SIM_VIEW_H
#define SIM_VIEW_H
#include "simulation.h"
#include <SFML/Graphics.hpp>


/// The game's main window
/// Displays the game class
class sim_view
{
public:

  sim_view(simulation start_simulation = simulation(), float start_zoom = 10, float zoom_step = 0.002f, double scale = 10);
  ~sim_view();

  /// Show one frame
  void show() noexcept;

  /// Run the game until the window is closed
  void exec() noexcept;

  /// Processes events in game and ouputs false if quit
  /// is inputted
  bool process_events();

  ///Gets the value of the zoom
  float get_max_zoom() const noexcept {return m_max_zoom;}

  ///Gets a const ref to m_sim
  const simulation& get_sim() const noexcept {return m_sim; }

  ///Gets const ref to m_view
  const sf::View& get_view() const noexcept {return  m_view;}

  ///Gets ref to m_view
  sf::View& get_view() noexcept {return  m_view;}

  ///Gets constant ref to sf::RenderWindow m_window
  const sf::RenderWindow& get_window() const noexcept {return m_window;}

  ///Get zoom step
  float get_zoom_step() const noexcept {return m_zoom_step;}


  ///Sets the sign of the zoom
  void set_zoom_dir(float dir) noexcept;

  ///Processes input for zooming
  void zoom_input(const sf::Event& event) noexcept;


private:
  /// The simulation logic
  simulation m_sim;

  /// The window to draw to
  sf::RenderWindow m_window;

  ///The zoom at which the simulation is shown
  /// in the beginning and the max level of zoom
  float m_max_zoom;

  ///The times the size of individuals and their coordinates are increased
  /// This is necessary to not show fractions of pixels, default 10
  double m_scale;

  ///The view of the simulation
  sf::View m_view;

  ///The direction towards which the view will zoom (+)in or (-)out
  int m_zoom_dir = 0;

  ///The step at which it is possible to zoom in or our
  float m_zoom_step;

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

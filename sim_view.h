#ifndef SIM_VIEW_H
#define SIM_VIEW_H
#include "simulation.h"
#include <SFML/Graphics.hpp>


/// The game's main window
/// Displays the game class
class sim_view
{
public:

  sim_view(simulation start_simulation = simulation(), float start_zoom = 10, float zoom_step = 0.01f, float pan_step = 2.f, double scale = 10);
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

  ///Gets the state of the flag that signals if it is requested to zoom in
  bool get_zoom_in() const noexcept {return m_zoom_in;}

  ///Gets the state of the flag that signals if it is requested to zoom out
  bool get_zoom_out() const noexcept {return m_zoom_out;}

  ///Processe when keyboard input for panning is given
  void pan_k_input_starts(const sf::Event& event) noexcept;

  ///Processe when keyboard input for panning is stopped
  void pan_k_input_ends(const sf::Event& event) noexcept;

  ///Pans camera following input from keyboard
  void k_pan() noexcept;

  ///Zooms following input from keyboard
  void k_zoom() noexcept;

  ///Processes when keyboard input for zooming is given
  void zoom_k_input_starts(const sf::Event& event) noexcept;

  ///Processes when keyboard input for zooming is stopped
  void zoom_k_input_ends(const sf::Event& event) noexcept;



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

  ///Boolean that signals if it is requested to zoom out
  bool m_zoom_out = false;

  ///Boolean that signals if it is requested to zoom out
  bool m_zoom_in = false;

  ///Boolean that signals if it is requested to pan the camera up
  bool m_pan_up = false;

  ///Boolean that signals if it is requested to pan the camera down
  bool m_pan_down = false;

  ///Boolean that signals if it is requested to pan the camera left
  bool m_pan_left = false;

  ///Boolean that signals if it is requested to pan the camera right
  bool m_pan_right = false;

  ///The step at which it is possible to move around the camera
  float m_pan_step;

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

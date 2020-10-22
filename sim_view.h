#ifndef SIM_VIEW_H
#define SIM_VIEW_H
#include "grid_view.h"
#include "simulation.h"
#include <SFML/Graphics.hpp>


/// The game's main window
/// Displays the game class
class sim_view
{
public:

    sim_view(float start_zoom = 10,
             float zoom_step = 0.01f,
             float pan_step = 2.f,
             float scale = 10);
    ~sim_view();

    /// Run the sim until the window is closed
    /// or simulation is finished
    void exec(simulation& s) noexcept;

    ///Same as exec_cycle for simulation but
    /// processes event and shows() every tick
    bool exec_cycle_visual(simulation &s) noexcept;

    ///Returns const ref to m_grid_view
    const  grid_view& get_grid_view() const noexcept {return m_grid_view;}

    ///Returns  ref to m_grid_view
    grid_view& get_grid_view() noexcept {return m_grid_view;}

    ///Gets the value of the zoom
    float get_max_zoom() const noexcept {return m_max_zoom;}

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

    ///Pans camera following input from keyboard
    void k_pan() noexcept;

    ///Zooms following input from keyboard
    void k_zoom() noexcept;

    ///Processe when keyboard input for panning is given
    void pan_k_input_starts(const sf::Event& event) noexcept;

    ///Processe when keyboard input for panning is stopped
    void pan_k_input_ends(const sf::Event& event) noexcept;

    ///Creates a vector of sf::CircleShape that will be used to draw the population
    void prepare_pop(const simulation &s) noexcept;

    /// Processes events in game and ouputs false if quit
    /// is inputted
    bool process_events();

    ///Updates the vector of shapes representing individuals if new individuals are added
    void update_pop(const simulation &s) noexcept;

    ///Updates the vector of indexes representing individuals if new individuals are added
    void update_indexes() noexcept;

    /// Show one frame
    void show(const simulation &s) noexcept;

    ///Parses input to see if food concentrations need to be shown on the grid
    void show_food_input(sf::Event& event);

    ///Parses input to see if metab concentrations need to be shown on the grid
    void show_metab_input(sf::Event& event);

    ///Parses input to see if simulation logic has to stop
    void start_stop_input(const sf::Event& event) noexcept;

    ///Parses input to see if simulation logic has to stop
    void stop_input_ends(const sf::Event& event) noexcept;

    ///Processes when keyboard input for zooming is given
    void zoom_k_input_starts(const sf::Event& event) noexcept;

    ///Processes when keyboard input for zooming is stopped
    void zoom_k_input_ends(const sf::Event& event) noexcept;



private:

    ///The class that defines the drawable object representing the environment grid
    grid_view m_grid_view;

    ///The zoom at which the simulation is shown
    /// in the beginning and the max level of zoom
    float m_max_zoom;

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

    ///Vector containing the circle shapes representing each individual in the population
    std::vector<sf::CircleShape> m_pop_shapes;

    ///Vector containing the fonts representing each individual index
    std::vector<sf::Text> m_pop_indexes;

    ///The times the size of individuals and their coordinates are increased
    /// This is necessary to not show fractions of pixels, default 10
    float m_scale;

    ///Flag that signal if the logic is active or stopped and sim_view only shows
    /// Simulation strts stopped;
    bool m_stop = false;

    /// The window to draw to
    sf::RenderWindow m_window;

    ///The view of the simulation
    sf::View m_view;

    ///The step at which it is possible to zoom in or our
    float m_zoom_step;

    /// Draws individuals
    void draw_inds(const simulation &s) noexcept;

    ///Draws the background
    void draw_background(const simulation &s) noexcept;

    ///Draws food
    void draw_food() noexcept;

    /// Draws projectiles
    void draw_metabolite() noexcept;
};

void test_sim_view();

#endif // SIM_VIEW_H

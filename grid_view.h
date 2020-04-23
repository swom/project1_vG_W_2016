#ifndef GRID_VIEW_H
#define GRID_VIEW_H
#include "SFML/Graphics.hpp"
#include "environment.h"

class grid_view : public sf::Drawable, public sf::Transformable
{
public:
  grid_view(float scale);

  ///Gets the size of the vertex_array
  int get_grid_vert_size() const noexcept {return static_cast<int>(m_grid_vertices.getVertexCount()) ;}

  ///Gets the const the vertex_array
  const sf::VertexArray& get_vertex_grid() const noexcept {return m_grid_vertices;}

  ///Returns the flag signalling if food is going to be shown on the grid
  bool is_food_shown() const noexcept {return m_food_is_shown;}

  ///Returns the flag signalling if food is going to be shown on the grid
  bool is_metab_shown() const noexcept {return m_metab_is_shown;}

  ///Prepares the VertexArray of quads and lines based on the size and food content of the environmental grid
  void prepare_grid(const std::vector<env_grid_cell> &env_grid) noexcept;

  ///Prepares the VertexArray of quads based on the size and food content of the environmental grid
  void prepare_grid_quads(const std::vector<env_grid_cell>& env_grid) noexcept;

  ///Prepares the VertexArray of lines based on the size of the environmental grid
  void prepare_grid_lines(const std::vector<env_grid_cell> &env_grid) noexcept;

  ///Sets the flag that signal that food should be shown to true
  ///and the one that ensure that metabolite should be shown to false
  void show_food_conc() noexcept {m_food_is_shown = true; m_metab_is_shown = false;}

  ///Sets the flag that signal that food should be shown to true
  ///and the one that ensure that metabolite should be shown to false
  void show_metab_conc() noexcept {m_food_is_shown = false; m_metab_is_shown = true;}

  ///Updates the color transparency of the quads based on the amount of food contained in them
  void update_grid_quads(const std::vector<env_grid_cell> &env_grid) noexcept;


private:

  ///Flag signalling that the food concentrations are being shown on the grid
  bool m_food_is_shown;

  ///Flag signalling that the metab concentrations are being shown on the grid
  bool m_metab_is_shown;

  ///Scale at which to draw things, given by sim_view at construction
  float m_scale;

  ///Vertex array to draw the squares
  sf::VertexArray m_grid_vertices;

  ///Vertex array to draw the lines between the squares
  sf::VertexArray m_grid_lines;

  virtual void draw(sf::RenderTarget& target, sf::RenderStates states) const
  {
    // apply the transform
    states.transform *= getTransform();


    // our particles don't use a texture
    states.texture = nullptr;

    // draw the vertex array of quads
    target.draw(m_grid_vertices, states);

    // draw the vertex array of lines
    target.draw(m_grid_lines, states);
  }

};

void test_grid_view() noexcept;
#endif // GRID_VIEW_H

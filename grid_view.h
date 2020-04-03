#ifndef GRID_VIEW_H
#define GRID_VIEW_H
#include "SFML/Graphics.hpp"
#include "environment.h"

class grid_view : public sf::Drawable, public sf::Transformable
{
public:
  grid_view(float scale);

  ///Gets the sie of the vertex_array
  int get_grid_vert_size() const noexcept {return static_cast<int>(m_grid_vertices.getVertexCount()) ;}

  ///Prepares the VertexArray of quads and lines based on the size and food content of the environmental grid
  void prepare_grid(const std::vector<env_grid_cell> &env_grid) noexcept;

  ///Prepares the VertexArray of quads based on the size and food content of the environmental grid
  void prepare_grid_quads(const std::vector<env_grid_cell>& env_grid) noexcept;

  ///Prepares the VertexArray of lines based on the size of the environmental grid
  void prepare_grid_lines(const std::vector<env_grid_cell> &env_grid) noexcept;

  ///Updates the color transparency of the quads based on the amount of food contained in them
  void update_grid_quads(const std::vector<env_grid_cell> &env_grid) noexcept;


private:

  ///Scale at which to draw things, given by sim_view at construction
  float m_scale;

  sf::VertexArray m_grid_vertices;

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

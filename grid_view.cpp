#include "grid_view.h"

#include <cassert>
#include <cmath>

grid_view::grid_view(float scale):
    m_scale(scale)
{

}

void grid_view::prepare_grid_lines(const std::vector<env_grid_cell> &env_grid) noexcept
{
    //Starting values
    auto grid_side = static_cast<float>((std::sqrt(env_grid.size()))) * m_scale;
    auto half_grid_side = grid_side / 2;
    m_grid_lines.setPrimitiveType(sf::Lines);
    m_grid_lines.resize(static_cast<size_t>(grid_side) * 4);
    auto x = -half_grid_side;
    auto y = - half_grid_side;

    for(size_t i = 0; i != static_cast<size_t>(grid_side); i++)
    {
        sf::Vertex* line = &m_grid_lines[i * 4];
        line[0].position = sf::Vector2f(x, - half_grid_side);
        line[1].position = sf::Vector2f(x, half_grid_side);
        line[2].position = sf::Vector2f(- half_grid_side, y);
        line[3].position = sf::Vector2f(half_grid_side, y);
        for( size_t j = 0; j != 4; j++)
        {
            line[j].color = sf::Color::Black;
        }
        x += m_scale;
        y += m_scale;
    }
}


void grid_view::prepare_grid_quads(const std::vector<env_grid_cell> &env_grid) noexcept
{

    //Starting values
    auto grid_side = static_cast<float>((std::sqrt(env_grid.size()))) * m_scale;
    auto half_grid_side = grid_side / 2;

    // resize the vertex array to fit the grid size
    m_grid_vertices.setPrimitiveType(sf::Quads);
    m_grid_vertices.resize(env_grid.size() * 4);

    auto x = -half_grid_side;
    auto y = - half_grid_side;

    for(size_t i = 0; i != env_grid.size(); i++)
    {
        // get a pointer to the current grid_cell's quad
        sf::Vertex* quad = &m_grid_vertices[i * 4];

        // define its 4 corners
        quad[0].position = sf::Vector2f(x, y);
        quad[1].position = sf::Vector2f(x + m_scale,  y);
        quad[2].position = sf::Vector2f(x + m_scale,  y + m_scale);
        quad[3].position = sf::Vector2f(x,  y + m_scale);

        //Define color
        for( size_t j = 0; j != 4; j++)
        {
            quad[j].color = sf::Color::Green;
            auto food_ratio = env_grid[i].get_food() / env_grid[i].get_max_food();
            quad[j].color.a = static_cast<sf::Uint8>(food_ratio * 255);
        }

        x += m_scale;
        if( x - half_grid_side < 0.000001f && x - half_grid_side > -0.000001F )
        {
            x = - half_grid_side;
            y += m_scale;
        }
    }
}

void grid_view::prepare_grid(const std::vector<env_grid_cell> &env_grid) noexcept
{
    prepare_grid_lines(env_grid);
    prepare_grid_quads(env_grid);
}


void grid_view::update_grid_quads(const std::vector<env_grid_cell> &env_grid) noexcept
{
    //All vertexes have the same color
    auto color = m_grid_vertices[0].color;

    for(size_t i = 0; i != env_grid.size(); i++)
    {
        // get a pointer to the current grid_cell's quad
        sf::Vertex* quad = &m_grid_vertices[i * 4];

        //Define color
        if(m_food_is_shown)
        {
            for( size_t j = 0; j != 4; j++)
            {
                if(color.g != sf::Color::Green.g)
                {quad[j].color = sf::Color::Green;}

                auto food_ratio = env_grid[i].get_food() / env_grid[i].get_max_food();
                if(food_ratio > 1)
                    food_ratio = 1;
                quad[j].color.a = static_cast<sf::Uint8>(food_ratio * 255);
            }
            continue;
        }

        for( size_t j = 0; j != 4; j++)
        {
            if(color.g == sf::Color::Green.g)
            {quad[j].color = sf::Color(255,192,203);}
            auto metab_ratio = env_grid[i].get_metab() / env_grid[i].get_max_metab();
            if(metab_ratio > 1)
                metab_ratio = 1;
            quad[j].color.a = static_cast<sf::Uint8>(metab_ratio * 255);
        }

    }
}

void test_grid_view() noexcept
{

    //When food is shown the color of the quad vertexes is green
    {
        std::vector<env_grid_cell> grid{10,10};
        grid_view gv{10};
        gv.prepare_grid(grid);
        gv.is_food_shown();
        gv.update_grid_quads(grid);
        for(size_t i = 0; i != gv.get_vertex_grid().getVertexCount(); i++)
        {
            assert(gv.get_vertex_grid()[i].color.g == sf::Color::Green.g);
        }
    }
    //When metab is shown the color of the quad vertexes is pink
    {
        std::vector<env_grid_cell> grid{10,10};
        grid_view gv{10};
        gv.prepare_grid(grid);
        gv.is_metab_shown();
        gv.update_grid_quads(grid);
        sf::Color pink(255,192,203);
        for(size_t i = 0; i != gv.get_vertex_grid().getVertexCount(); i++)
        {
            assert(gv.get_vertex_grid()[i].color.r == pink.r);
            assert(gv.get_vertex_grid()[i].color.b == pink.b);
            assert(gv.get_vertex_grid()[i].color.g == pink.g);
        }
    }

}

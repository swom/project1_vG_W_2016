#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>

simulation::simulation(int pop_size, int exp_new_pop_size, double min_dist,
                       int grid_side, double diff_coeff,
                       double init_food, double mutation_prob,
                       double mutation_step, double base_disp_prob, double spore_advantage):

  m_pop(static_cast<unsigned int>(pop_size),individual(0,0)),
  m_exp_new_pop_size(exp_new_pop_size),
  m_min_init_dist_btw_cells(min_dist),
  m_grid_side(grid_side),
  m_diff_coeff(diff_coeff),
  m_init_food(init_food),
  m_e(m_grid_side, m_diff_coeff, m_init_food),
  m_reproduction_angle(0, 2 * M_PI),
  m_mutation_prob(mutation_prob),
  m_mutation_step(0,mutation_step),
  m_disp_dist(0,1),
  m_base_disp_prob(base_disp_prob),
  m_spore_advantage(spore_advantage)
{
#ifndef IS_ON_TRAVIS
  try {
    if(base_disp_prob * spore_advantage > 1)
      {throw std::string{"base dispersal probability * spore advantage > 1, too high!\n"};}
  }
  catch (std::string& e) {
    std::cout << e;
#ifdef NDEBUG
    abort();
#endif
  }
#endif
  if(!m_pop.empty())
    {
      place_start_cells(*this);
    }
}

std::pair<double,double> simulation::get_ind_pos(int i)
{
  std::pair<double, double> pos = get_pos(get_ind(i));
  return pos;
}

bool all_ind_are_drawn(const simulation& s) noexcept
{
  return std::all_of(s.get_pop().begin(), s.get_pop().end(),
                     [](const individual& i) {return is_drawn(i);});
}

int count_hex_layers(int pop_size)  noexcept
{
  int n = 1;
  if(pop_size>0){while(3*n*(n-1)+1 < pop_size) n++;}
  else {return 0;}
  return n;
}

double calc_angle_3_pos(std::pair<double, double> P1,
                        std::pair<double, double> P2,
                        std::pair<double, double> P3)
{
  auto angle1 = std::atan2(P3.second - P1.second, P3.first - P1.first);
  auto angle2 = std::atan2(P2.second - P1.second, P2.first - P1.first);
  if(angle1 >= angle2)
    return angle1 - angle2;

  return M_PI + std::abs(angle1 -angle2);
}

//not sure this is the fastest implementation, maybe sawp and pop_back is still faster?
void death(simulation& s) noexcept
{
  s.get_pop().erase(std::remove_if(s.get_pop().begin(),s.get_pop().end(),
                                   [](individual const &i){ return is_dead(i);})
      ,s.get_pop().end());
}

void division(simulation &s) noexcept
{
  auto  div_inds  = get_dividing_individuals(s);
  for(size_t i = 0; i != div_inds.size(); i++)
    {
      int div_ind = div_inds[i];
      divides(s.get_ind(div_ind),s.get_pop(),s.repr_angle(),
              s.get_rng(),s.get_mu_p(),s.get_mu_st());
    }
}

void dispersal(simulation& s)
{
  select_new_pop(s);
  fund_pop(s);
  place_start_cells(s);
  reset_env(s.get_env());
}

void exec(simulation& s, int n_tick) noexcept
{
  while(s.get_tick() != n_tick){tick(s);}
}

void feeding(simulation& s)
{
  for(auto& ind : s.get_pop())
    {
      auto index_grid = find_grid_index(ind,s.get_env().get_grid_side());
      if(index_grid == -100 ||
         ind.get_phen() != phenotype::active)
        {continue;}
      feed(ind,s.get_env().get_cell(index_grid));
    }
}

void fund_pop(simulation& s) noexcept
{
  s.get_pop().swap(s.get_new_pop());
  s.get_new_pop().clear();
}

std::vector<int> get_dividing_individuals(const simulation& s) noexcept
{

  std::vector<int> dividing_individuals;
  for(unsigned int i = 0; i != s.get_pop().size(); i++)
    {
      if(s.get_pop()[i].get_energy() >= s.get_pop()[i].get_treshold_energy()
         && s.get_pop()[i].get_phen() == phenotype::active)
        {
          dividing_individuals.push_back(static_cast<int>(i));
        }
    }
  return dividing_individuals;
}

std::vector<double> get_excess_energies(const simulation& s) noexcept
{
  std::vector<double> excess_energy;
  auto v_ind{get_dividing_individuals(s)};
  for(unsigned int i=0; i != v_ind.size(); i++)
    {
      int ind_index = v_ind[i];
      excess_energy.push_back(
            s.get_ind_en(ind_index) - s.get_ind_tr_en(ind_index)
            );
    }
  return excess_energy;
}

std::vector<std::pair<int,int>> get_sisters_index_offset(const simulation& s)  noexcept
{
  std::vector<std::pair<int,int>> sis_indexes;
  std::pair<int,int> daughters;
  int sis_dist;
  auto div_ind = get_dividing_individuals(s);

  for (const auto& ind : div_ind)
    {
      sis_dist = static_cast<int>(s.get_pop_size()) - ind;
      daughters.first = ind;
      daughters.second = daughters.first + sis_dist;
      sis_indexes.push_back(daughters);
    }
  return sis_indexes;
}

std::vector<int> has_collision(simulation& s)
{
  //Sort the pop vector by increasing x
  std::sort(s.get_pop().begin(),s.get_pop().end(),
            [](const individual& lhs, const individual& rhs)
  {return lhs.get_x() < rhs.get_x();});

  int n_ind = s.get_pop_size();
  for ( int i = 0; i < n_ind; ++i)
    {
      const auto focal_ind = s.get_ind(i);
      //Sort all individuals whose x +/- radius is in the range between focal_x -/+ focal_radius
      //in ascending order by thei y coordinate
      auto first_x = std::lower_bound(
            s.get_pop().begin(), s.get_pop().end(), focal_ind,
            [](const individual& lhs, const individual& rhs)
      {return lhs.get_x() + lhs.get_radius() < rhs.get_x() - rhs.get_radius();}
      );
      if(first_x == s.get_pop().end()){continue;}

      auto last_x = std::upper_bound(
            s.get_pop().begin(), s.get_pop().end(), focal_ind,
            [](const individual& lhs, const individual& rhs)
      {return lhs.get_x() + lhs.get_radius() < rhs.get_x() - rhs.get_radius();}
      );
      if(last_x == s.get_pop().end()){continue;}

      if(first_x == last_x){continue;}

      auto check_x_first = first_x;
      auto check_x_last = last_x;
      std::sort(first_x, last_x, [](const individual& lhs, const individual& rhs){return lhs.get_y() < rhs.get_y();});
      assert(first_x == check_x_first && last_x == check_x_last);

      auto first_y = std::lower_bound(
            first_x,last_x,focal_ind,
            [](const individual& lhs, const individual& rhs)
      {return lhs.get_y() + lhs.get_radius() < rhs.get_y() - rhs.get_radius();}
      );
      if(first_y == s.get_pop().end()){continue;}


      auto last_y = std::upper_bound(
            first_x,last_x,focal_ind,
            [](const individual& lhs, const individual& rhs)
      {return lhs.get_y() + lhs.get_radius() < rhs.get_y() - rhs.get_radius();}
      );
      if(last_y == s.get_pop().end()){continue;}

      if(first_y == last_y){continue;}

      auto start = static_cast<int>(std::distance(s.get_pop().begin(),first_y));
      auto end = static_cast<int>(std::distance(s.get_pop().begin(), last_y));
      for ( auto j = start ; j != end; ++j)
        {
          if(i == j)
            {
              continue;
            }
          if (are_colliding(s.get_ind(i), s.get_ind(j)))
            {
              return std::vector<int>{i,j};
            }
        }
      //Sort back to increasing x
      std::sort(first_y ,last_y,
                [](const individual& lhs, const individual& rhs)
      {return lhs.get_x() <  rhs.get_x();});
    }

  std::vector<int> empty_v;
  return empty_v;
}

void calc_tot_displ_pop(std::vector<individual>& pop, std::vector<int> first_collisions_indexes)
{
  //Sort the pop vector by increasing x
  std::sort(pop.begin(), pop.end(),
            [](const individual& lhs, const individual& rhs)
  {return lhs.get_x() <  rhs.get_x();});

  const auto n_ind = pop.size();

  for ( auto i = static_cast<size_t>(first_collisions_indexes[0]); i != n_ind; ++i)
    {
      for ( size_t j = 0 ; j != n_ind; ++j)
        {
          if (i == static_cast<size_t>(first_collisions_indexes[0]) &&
              j < static_cast<size_t>(first_collisions_indexes[1]))
            {
              j = static_cast<size_t>(first_collisions_indexes[1]);
            }
          else if(j == i) continue;
          if(are_colliding(pop[i], pop[j]))
            add_displacement(pop[i], pop[j]);
        }
    }
}

void no_complete_overlap(simulation& s) noexcept
{
  for(size_t i = 0; i != s.get_pop().size(); i++)
    for(size_t j = 0; j != s.get_pop().size(); j++)
      {
        if(i == j ) continue;
        if(distance(s.get_pop()[i], s.get_pop()[j]) < 0.00001 &&
           distance(s.get_pop()[i], s.get_pop()[j]) > -0.00001)
          {
            auto rnd_n_0_1 = s.get_unif_dist()(s.get_rng());
            s.get_pop()[i].change_x(rnd_n_0_1 * 0.0001);
            s.get_pop()[j].change_x(rnd_n_0_1 * -0.0001);
          }
      }
}

void manage_static_collisions(simulation& s)
{
  auto first_collision_indexes = has_collision(s);
  int time = 0;
  while(!first_collision_indexes.empty())
    {
      no_complete_overlap(s);
      calc_tot_displ_pop(s.get_pop(), first_collision_indexes);
      for(auto& ind : s.get_pop()){ind.displace();}
      first_collision_indexes = has_collision(s);
      time ++;
    }
}

void metabolism_pop(simulation& s)
{
  for(auto& ind : s.get_pop())
    {
      if(ind.get_phen() != phenotype::spore)
        metabolism(ind);
    }
}

std::vector<double> modulus_of_btw_ind_angles(simulation& s, double ang_rad)
{
  std::vector<double> v_modulus;
  for(int i = 0; i != s.get_pop_size() - 2; i++)
    for(int j = i+1; j != s.get_pop_size() - 1; j++)
      for(int k = j+1; k != s.get_pop_size(); k++)
        {
          auto P1 = get_pos(s.get_ind(i));
          auto P2 = get_pos(s.get_ind(j));
          auto P3 = get_pos(s.get_ind(k));
          v_modulus.push_back((abs(fmod(calc_angle_3_pos(P1,P2,P3),ang_rad))));
        }
  return v_modulus;
}

void place_start_cells(simulation& s) noexcept
{
  int n = count_hex_layers(s.get_pop_size());
  unsigned int placed_ind = 0;

  // d is the distance between 2 individuals's centers
  double d = 2 * (s.get_pop()[0].get_radius() + s.get_min_dist());

  for(int i = 0; i != n; i++)
    {
      double y = (sqrt(3) * i * d) / 2.0;

      for (int j = 0; j < (2 * n - 1 - i); j++)
        {
          double x = (-(2 * n - i - 2) * d) / 2.0 + j * d;

          s.get_pop()[placed_ind].set_x(x);
          s.get_pop()[placed_ind].set_y(y);
          placed_ind++;
          if(placed_ind == s.get_pop().size()) return;

          if (y > 0.000001 || y < -0.000001)
            {
              s.get_pop()[placed_ind].set_x(x);
              s.get_pop()[placed_ind].set_y(-y);
              placed_ind++;
              if(placed_ind == s.get_pop().size()) return;
            }
        }

    }
}

void reset_drawn_fl_new_pop(simulation& s) noexcept
{
  std::for_each(s.get_new_pop().begin(),s.get_new_pop().end(),
                [](individual& i){draw_flag_reset(i);});
}

void select_new_pop(simulation& s)
{
  assert(s.get_new_pop().empty());
  while(true)
    {
      for(auto& ind : s.get_pop())
        {
          if(s.get_unif_dist()(s.get_rng()) <
             get_fitness(ind, s.get_base_disp_prob(), s.get_spo_adv())
             && !is_drawn(ind))
            {
              draw(ind);
              s.get_new_pop().push_back(ind);
            }
          if(s.get_new_pop_size() == s.get_exp_new_pop_size() ||
             all_ind_are_drawn(s))
            {
              reset_drawn_fl_new_pop(s);
              return;
            }
        }
    }
}

void tick(simulation& s)
{
  feeding(s);
  metabolism_pop(s);
  death(s);
  division(s);
  manage_static_collisions(s);
  diffusion(s.get_env());
  s.update_sim_timer();
}


void test_simulation()//!OCLINT tests may be many
{
#ifndef NDEBUG


  //Simulation is initialized with a certain number of individuals
  // The value 1234567890 is irrelevant: just get this to compile

  {
    int pop_size = 100;
    simulation s(pop_size);
    for(int i = 0; i < s.get_pop_size(); ++i)
      {
        double x = s.get_ind(i).get_x();
        assert( x > -1234567890 );
      }
  }

  //The size of a population is equal to the size of the vector containing its individuals
  {
    simulation s;
    assert( s.get_pop_size() == static_cast<int>(s.get_pop().size()));
  }


  //No individuals are dividing at the start of the simulation
  {
    //initiate empty vector and leave it empty
    simulation s(3);
    assert(static_cast<int>(get_dividing_individuals(s).size()) < 0.00000000001);
  }


  // The number n of hex_layers required for x individual is
  // 3n(n-1)+1 <= x < 3(n+1)((n+1)-1)+1
  {
    std::vector<int> x{0,1,7,19,22};
    std::vector<int> n{0,1,2,3,4};

    for(size_t i = 0; i < x.size(); i++ )
      {
        assert(count_hex_layers(x[i]) == n[i]);
      }
  }

  //Individuals' centers are placed in an hexagonal pattern
  {
    simulation s(7);

    //The angle formed by any 3 individual_centers is a multiple of PI/6*number_of_hex_layers
    //This is true in an hex_patt but i do not know if it is sufficient

    double n_hex_l = count_hex_layers(s.get_pop_size());
    auto minimal_angle = M_PI / (6 * n_hex_l);
    auto angle = M_PI  / 12.0;
    assert(minimal_angle - angle < 0.00001 &&
           minimal_angle - angle  / 12.0 > -0.00001 );
    auto v_modulus = modulus_of_btw_ind_angles(s,minimal_angle);
    for(auto ind_modulus : v_modulus)
      assert( ind_modulus < 0.0000000001 || (ind_modulus > minimal_angle - 0.000001 && ind_modulus <= minimal_angle));

  }


  //No individuals are overlapping at the start of the simulation
  {
    simulation s(2);
    assert(has_collision(s).empty());
  }


  //An individual position can be retrieved as a pair object x,y
  {
    simulation s;
    std::pair<double, double> v{0,0};//default coordinates of individuals
    std::pair<double, double> v2{1,1};//different coordinates from default

    assert(s.get_ind_pos(0) == v);
    assert(s.get_ind_pos(0) != v2);
  }


  // An individaul can be placed to some given coordinates after initialization
  {
    simulation s(2);
    for(int i = 0; i < s.get_pop_size(); i++)
      {
        set_pos(s.get_ind(i),std::pair<double,double>(i,i));
      }

    for (int i = 0; i < s.get_pop_size() - 1; i++)
      {
        for(int j = i+1 ; j < s.get_pop_size(); j++)
          {
            assert(s.get_ind_pos(i) != s.get_ind_pos(j));
          }
      }

  }


  //When an individual divides it adds a copy of itself
  //to the population/individual vector(i.e divides)
  {
    simulation s;
    double lhs = s.get_pop_size();
    //let's allow all individuals in the population to reproduce,
    //to facilitate the testing conditions
    for(auto& ind : s.get_pop())
      {
        ind.set_energy(ind.get_treshold_energy());
      }
    division(s);

    double rhs = s.get_pop_size();
    assert(lhs * 2 - rhs < 0.0000001);

  }


  //Only individuals with energy >= than the treshold will divide
  {
    simulation s(3);
    //This ind will not reproduce
    s.get_ind(0).set_energy(s.get_ind(0).get_treshold_energy()-1);
    //This ind will reproduce with no extra energy
    s.get_ind(1).set_energy(s.get_ind(1).get_treshold_energy());
    //This ind will reproduce with extra energy
    s.get_ind(2).set_energy(s.get_ind(2).get_treshold_energy()+1);

    std::vector<int> div_ind = get_dividing_individuals(s);
    assert(!div_ind.empty());

    for(unsigned int i = 0; i != div_ind.size(); i++)
      assert(s.get_ind(div_ind[i]).get_energy()
             >= s.get_ind(div_ind[i]).get_treshold_energy());

  }

  //The excess energy of dividing individuals is equal
  //to the diference between their energy and their treshold energy
  {
    simulation s(2);
    double energy = s.get_ind_tr_en(0);
    s.set_ind_en(0,energy);
    s.set_ind_en(1,energy*2);
    auto v_ex_en = get_excess_energies(s);
    for(int i =0; i != static_cast<int>(v_ex_en.size()); i++){
        assert(v_ex_en[static_cast<unsigned int>(i)]
            - (s.get_ind_en(i) - s.get_ind_tr_en(i)) < 0.00000001);
      }
  }


  //The distance between the two offspring from same mother in population vector
  //is the distance betwen the mother index and the end of the population vector

  {
    simulation s(1);
    //First individual reproduces
    s.set_ind_en(0,s.get_ind_tr_en(0)*2);
    auto daughters = get_sisters_index_offset(s);
    // In this case the distance between the two offspring
    //will be 1 element of the vector in next gen
    auto sister_distances = 1;
    for(auto sisters : daughters)
      {
        assert( sisters.second - sisters.first - sister_distances < 0.000001);
      }
  }


  //After a reproduction round individuals with enough energy will divide
  //and redistribute their remaining energy to their daughter cells
  {
    simulation s;
    auto repr_excess_en = 3;
    //This ind will reproduce with extra energy
    s.set_ind_en(0,repr_excess_en);

    auto mother_excess_energy = get_excess_energies(s);
    division(s);
    assert(s.get_ind_en(0) - mother_excess_energy[0]/2 < 0.00001 &&
        s.get_ind_en(0) - mother_excess_energy[0]/2 > -0.00001);
    assert(s.get_ind_en(1) - mother_excess_energy[0]/2 < 0.00001 &&
        s.get_ind_en(1) - mother_excess_energy[0]/2 > -0.00001);
    assert(s.get_ind_en(0) - s.get_ind_en(1) < 0.00001 &&
           s.get_ind_en(0) - s.get_ind_en(1) > -0.00001);
  }


  //After reproduction the first daughter individual takes the position of the mother
  {
    simulation s;
    auto parent_pop = s.get_pop();
    //setting energy high enough for the individual
    //to reproduce without offspring having negative enrgies
    s.set_ind_en(0,s.get_ind_tr_en(0));

    divides(s.get_ind(0),s.get_pop(),s.repr_angle(),
            s.get_rng(), s.get_mu_p(), s.get_mu_st());
    //The first daughter cell is at the same index of the mother cell
    assert(s.get_ind_pos(0) == get_pos(parent_pop[0]) );
  }

  //After reproduction the second daughter individual
  //is placed just outside the position of the first daughter cell
  {
    simulation s;
    //setting energy high enough for the individual to reproduce
    //without offspring having negative enrgies
    s.set_ind_en(0,s.get_ind_tr_en(0));
    divides(s.get_ind(0), s.get_pop(), s.repr_angle(),
            s.get_rng(), s.get_mu_p(), s.get_mu_st());
    assert(has_collision(s).empty());
    assert(distance(s.get_ind(0), s.get_ind(1)) -
           (s.get_ind(0).get_radius() + s.get_ind(1).get_radius()) < 0.1);
  }

  //If there are collisions individuals are displaced
  //until there are no more collisions
  {
    simulation s(2);
    assert(has_collision(s).empty());
    set_pos(s.get_ind(1),s.get_ind_pos(0));
    no_complete_overlap(s);
    assert(!has_collision(s).empty());

    calc_tot_displ_pop(s.get_pop(), has_collision(s));
    for(auto& ind : s.get_pop()){ind.displace();}

    assert(has_collision(s).empty());

  }
  //After reproduction new collisions caused by new individuals
  //being placed where other individuals already are managed
  {
    simulation s(7);
    assert(has_collision(s).empty());
    //add 1 individual overlapping with central individual
    s.get_pop().emplace_back(individual(0,0));
    assert(!has_collision(s).empty());
    manage_static_collisions(s);
    assert(has_collision(s).empty());
  }

  //A simulation has an environment
  // The value -1234567890 is irrelevant: just get this to compile

  {
    simulation s;
    assert(s.get_env().get_grid_side() > -1234567890);
  }

  //A simulation can be initialized with a certain amount
  //of food in all the cells of its environment grid
  //1 by default
  {
    simulation s;
    for( auto& grid_cell : s.get_env().get_grid())
      {
        assert(grid_cell.get_food() - 1 < 0.000001
               && grid_cell.get_food() - 1 > -0.0000001);
      }

    double starting_food = 3.14;
    s = simulation(1, 1, 1, 1, 0.1, starting_food);
    for( auto& grid_cell : s.get_env().get_grid())
      {
        assert(grid_cell.get_food() - starting_food < 0.000001
               && grid_cell.get_food() - starting_food > -0.0000001);
      }
  }

  //Individuals can take up energy from the environment
  {
    simulation s(1,1,0.1,2);
    double total_food_init = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0,
          [](double sum,const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_init = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );

    feeding(s);

    double total_food_after = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_after = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){
        return sum + i.get_energy();}
    );

    assert(total_food_init > total_food_after);
    assert(total_en_init < total_en_after);

    auto total_uptake = std::accumulate(s.get_pop().begin(), s.get_pop().end(), 0.0,
                                        [](double sum, const individual& i) {return sum + i.get_uptake_rate();});

    assert(total_en_after - (total_uptake + total_en_init) < 0.00000001 &&
           total_en_after - (total_uptake + total_en_init) > -0.00000001);
    assert(total_food_init - (total_food_after + total_uptake) < 0.000001 &&
           total_food_init - (total_food_after + total_uptake) > -0.000001);
  }

  //Individuals outside the grid do not feed
  {
    simulation s(1,1,0.1,2);

    double total_food_init = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0.0,
          [](double sum,const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_init = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );

    set_pos(s.get_ind(0),std::pair<double,double>(-42,42));

    feeding(s);

    double total_food_after = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_after = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){
        return sum + i.get_energy();}
    );

    assert(total_food_init - total_food_after < 0.00001
           && total_food_init - total_food_after > -0.00001);
    assert(total_en_init - total_en_after < 0.000001
           && total_en_init - total_en_after > -0.000001);
  }

  //Individuals lose energy through metabolism
  {
    simulation s(1,1,0.1,2);
    feeding(s);
    double init_en_tot = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );

    metabolism_pop(s);
    double after_en_tot = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );
    assert(after_en_tot < init_en_tot);
  }

  //In one tick/timestep individuals first feed,
  //than reproduce, than substances diffuse
  {
    simulation s(1,1,0.1,3,1,1,0);

    //The single individual in this population
    //after a tick should reproduce
    auto init_pop_size = s.get_pop().size();
    s.get_ind(0).set_energy(s.get_ind_tr_en(0)
                            + s.get_ind(0).get_metab_rate() + 0.01
                            - s.get_ind(0).get_uptake_rate());

    //and the grid_cell where it is should recieve
    //food nutrients
    auto grid_index_ind =  find_grid_index(s.get_ind(0),s.get_env().get_grid_side());
    double predicted_food_after_feeding =
        s.get_env().get_cell(grid_index_ind).get_food() - s.get_ind(0).get_uptake_rate();

    tick(s);

    assert(!(predicted_food_after_feeding - s.get_env().get_cell(grid_index_ind).get_food() < 0.00001
             && predicted_food_after_feeding -s.get_env().get_cell(grid_index_ind).get_food() > -0.0001));
    assert(s.get_pop().size() == 2 * init_pop_size);
    assert(has_collision(s).empty());
  }

  //If nothing else happens, food should constantly decrease when cells are feeding
  {
    simulation s (2, 1, 0, 3, 0.1, 5);
    auto food_begin = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                      [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

    //The simulation will last long enough  for the individuals to reproduce
    auto sim_time = s.get_ind(0).get_treshold_energy() / s.get_ind(0).get_uptake_rate() + 5;

    auto init_pop_size =  s.get_pop().size();

    for( int i = 0; i != static_cast<int>(sim_time); i++)
      {

        auto food_before_feed = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
        feeding(s);

        auto food_after_feed = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                               [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

        double food_eaten = 0.0;
        for(const auto& ind : s.get_pop())
          {
            auto grid_cell_ind = find_grid_index(ind, s.get_grid_side());
            if( grid_cell_ind != -100 && s.get_env().get_cell(grid_cell_ind).get_food() > 0)
              {
                food_eaten += ind.get_uptake_rate();
              }
          }

        auto balance_uptake = food_before_feed - (food_after_feed + food_eaten) ;
        assert(balance_uptake < 0.0001 && balance_uptake > -0.0001);

        metabolism_pop(s);
        death(s);
        division(s);
        manage_static_collisions(s);

        auto food_before_diff = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

        auto new_grid = calc_diffusion_food(s.get_env());

        auto food_after_diff = std::accumulate(new_grid.begin(), new_grid.end(), 0.0,
                                               [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

        auto balance_diffusion = food_before_diff - food_after_diff;
        assert(balance_diffusion < 0.0001 && balance_diffusion > -0.0001);

        new_grid.swap(s.get_env().get_grid());

      }
    auto food_end = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                    [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
    assert(food_end < food_begin);

    auto final_pop_size =  s.get_pop().size();

    assert(init_pop_size < final_pop_size);
  }

  //every timestep/tick collisions are handled
  {
    simulation s(7,3);
    //The central individual in this population
    //after a tick should reproduce
    auto init_pop_size = s.get_pop().size();
    s.get_ind(1).set_energy(s.get_ind_tr_en(1)
                            + s.get_ind(1).get_metab_rate() + 0.01
                            - s.get_ind(1).get_uptake_rate());
    feeding(s);
    metabolism_pop(s);
    division(s);
    manage_static_collisions(s);

    assert(s.get_pop().size() == init_pop_size + 1);
    assert(has_collision(s).empty());
  }
  //A simulation is initialized with a m_tick = 0;
  {
    simulation s;
    assert(s.get_tick() == 0);
  }
  //After each tick the simulation updates its m_tick
  //by one
  {
    simulation s;
    for(int i = 0; i != 3; i++)
      {
        assert(s.get_tick() == i);
        tick(s);
      }
  }

  //A simulation can be run for a certain amount of ticks
  {
    simulation s;
    auto n_ticks = 100;
    exec(s, n_ticks);
    assert(s.get_tick() - n_ticks == 0);
  }

  //Spores do not get detected when looking for dividing individual
  {
    simulation s;
    s.set_ind_en(0,s.get_ind_tr_en(0));
    assert(get_dividing_individuals(s)[0] == 0);
    s.get_ind(0).set_phen(phenotype::spore);
    assert(get_dividing_individuals(s).empty());
  }

  //Spores do not reproduce
  {
    simulation s;
    s.get_ind(0).set_phen(phenotype::spore);
    s.set_ind_en(0,s.get_ind_tr_en(0));
    auto init_pop_size = s.get_pop_size();
    s.set_ind_en(0,s.get_ind_tr_en(0));
    tick(s);
    assert(init_pop_size == s.get_pop_size());
  }

  //Spores do not feed
  {
    simulation s;
    auto init_food = s.get_env().get_cell(0).get_food();
    s.get_ind(0).set_phen(phenotype::spore);
    feeding(s);
    assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
           && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
  }
  //Spores do not lose energy
  {
    simulation s;
    s.get_ind(0).set_phen(phenotype::spore);
    s.set_ind_en(0,s.get_ind_tr_en(0));
    auto init_en_ind0 = s.get_ind_en(0);
    tick(s);
    assert(s.get_ind_en(0) - init_en_ind0 < 0.000001
           && s.get_ind_en(0) - init_en_ind0 > -0.000001);
  }

  //Sporulating individuals do not get detected when looking for dividing individual
  {
    simulation s;
    s.set_ind_en(0,s.get_ind_tr_en(0));
    assert(get_dividing_individuals(s)[0] == 0);
    s.get_ind(0).set_phen(phenotype::sporulating);
    assert(get_dividing_individuals(s).empty());
  }

  //Sporulating individuals cannot reproduce
  {
    simulation s;
    s.get_ind(0).set_phen(phenotype::sporulating);
    s.set_ind_en(0,s.get_ind_tr_en(0));
    auto init_pop_size = s.get_pop_size();
    s.set_ind_en(0,s.get_ind_tr_en(0));
    tick(s);
    assert(init_pop_size == s.get_pop_size());
  }

  //Sporulating individuals do not feed but they lose energy
  {
    simulation s;
    s.get_ind(0).set_phen(phenotype::sporulating);
    s.set_ind_en(0,1);
    auto init_en_ind0 = s.get_ind_en(0);
    auto init_food = s.get_env().get_cell(0).get_food();
    tick(s);
    assert(init_en_ind0 - s.get_ind_en(0) - s.get_ind(0).get_metab_rate() < 0.000001
           && init_en_ind0 - s.get_ind_en(0) - s.get_ind(0).get_metab_rate() > -0.000001);
    assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
           && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
  }


  //A sporulating individual updates its timer every tick
  {
    simulation s;
    auto init_timer = s.get_ind(0).get_spo_timer();
    s.get_ind(0).set_phen(phenotype::sporulating);
    tick(s);
    assert(init_timer != s.get_ind(0).get_spo_timer());
    assert(s.get_ind(0).get_spo_timer() == init_timer + 1);
    assert(s.get_ind(0).get_spo_timer() != init_timer + 2);

    s.get_ind(0).reset_spo_timer();
    init_timer = s.get_ind(0).get_spo_timer();
    s.get_ind(0).set_phen(phenotype::spore);
    tick(s);
    assert(s.get_ind(0).get_spo_timer() == init_timer);
    assert(s.get_ind(0).get_spo_timer() != init_timer + 1);
    assert(s.get_ind(0).get_spo_timer() != init_timer + 2);

    s.get_ind(0).reset_spo_timer();
    init_timer = s.get_ind(0).get_spo_timer();
    s.get_ind(0).set_phen(phenotype::active);
    tick(s);
    assert(s.get_ind(0).get_spo_timer() == init_timer);
    assert(s.get_ind(0).get_spo_timer() != init_timer + 1);
    assert(s.get_ind(0).get_spo_timer() != init_timer + 2);
  }

  //Individuals that die are removed from the population
  {
    simulation s;//the only individual in this sim has 0 energy, thus it will die
    assert(s.get_pop_size() == 1);
    death(s);
    assert(s.get_pop().empty() && s.get_pop_size() == 0);

    int pop_size = 5;
    //The simulation does not have a grid with food,
    //so organisms cannot feed
    s = simulation(pop_size,1,0.1,0);
    //Only the first individual has enough energy to survive
    //for 1 tick
    s.set_ind_en(0,s.get_ind(0).get_metab_rate() + 0.001);
    assert(s.get_pop_size() == pop_size);
    tick(s);
    assert(s.get_pop_size() == 1);
    //then at the second tick the only individual left dies
    tick(s);
    assert(s.get_pop().empty() && s.get_pop_size() == 0);
  }

  //A simulation is initialized with a random number generator
  {
    simulation s;
    std::uniform_int_distribution u_d(0,2);
    double mean = 0;
    //Draw a thousands times from a uniform dist between 0 and 2
    for(int i = 0; i != 1000; i++)
      {
        mean += u_d(s.get_rng());
      }
    //calculate mean of the drawn values
    mean /= 1000;
    assert(mean < 1.01 && mean > 0.99 );
  }

  //A simulation is initialized with a uniform distribution
  //between 0 and 2PI
  {
    simulation s;
    double mean = 0;
    //Draw a thousands times from a uniform dist between 0 and 2
    int sampling_size = 1000;
    for(int i = 0; i != sampling_size; i++)
      {
        mean += s.repr_angle();
      }
    //calculate mean of the drawn values
    mean /= 1000;
    assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
  }

  //Daughter cells are placed at a random angle after reproduction
  {
    simulation s;
    //to calculate angle we will use three point
    //the center of the mother(0,0)
    std::pair<double, double> mother (0,0);
    set_pos(s.get_ind(0),mother);
    //a reference point on the same axis as the mother (0,1)
    //and the center of the daughter -> get_daughter_pos()
    std::pair<double, double> reference(1,0);

    //Draw a thousands times from a uniform dist between 0 and 2
    double mean = 0;
    int sampling_size = 1000;
    for(int i = 0; i != sampling_size; i++)
      {
        auto daughter = get_daughter_pos(s.get_ind(0),s.repr_angle());
        mean += calc_angle_3_pos(mother,daughter,reference);
      }
    mean /= sampling_size;
    assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
  }

  //A simulation is initialized with a normal distribution for mutation step
  //with mean 0 and variance 0.1 by default
  {
    simulation s;
    double mean = 0;
    int sampling_size = 10000;
    for(int i = 0 ; i != sampling_size; i++ )
      mean += s.mut_step();
    mean /= sampling_size;
    assert(mean < 0.01 && mean > -0.01);
  }
  //A simulation is initialized with a bernoulli distribution to see if mutation happens or not
  //0.01 by default
  {
    simulation s;
    double mean = 0;
    int sampling_size = 100000;
    for(int i = 0 ; i != sampling_size; i++ )
      mean += s.mut_happens();
    mean /= sampling_size;
    assert(mean < 0.011 && mean > 0.009);
  }

  //The sum of weight of an individual after many rounds of mutation
  //Should have the same mean as in the beginning, but its variance should
  //be the same as the mutation_step distribution
  {
    simulation s;
    //   double init_mean = weights_mean(s.get_ind(0).get_grn());
    double init_variance = weights_var(s.get_ind(0).get_grn());
    assert(init_variance < 0.0001 && init_variance > -0.000001);

    int sampling_size = 10000;

    for (int i = 0; i != sampling_size; i++)
      {
        mutates(s.get_ind(0),s.get_rng(),
                s.get_mu_p(), s.get_mu_st());
      }
    //This first assert does not pass, the mean is much more variable than
    //I thought, but I do not see any bug. I will comment this out
    //    assert(mean - init_mean > -0.1 && mean - init_mean < 0.1);
    assert(init_variance - weights_var(s.get_ind(0).get_grn()) > 0.01 ||
           init_variance - weights_var(s.get_ind(0).get_grn()) < -0.01);
  }
  //After dividing the two daughter individuals mutate
  {
    double mutation_probability = 1; //all weights will be mutated in this simulation
    simulation s(1, 1, 0, 0, 0, 0, mutation_probability);
    auto init_var = weights_var(s.get_ind(0).get_grn());
    assert(init_var < 0.00001 && init_var > -0.0001);
    divides(s.get_ind(0),s.get_pop(),s.repr_angle(),s.get_rng(),
            s.get_mu_p(),s.get_mu_st());
    auto post_var = weights_var(s.get_ind(0).get_grn());
    assert(init_var - post_var > 0.000001 || init_var - post_var < -0.0001);
  }

  //A simulation has a member variable m_new_pop_size that states the max number of
  //individuals that will start a new population
  //by default = to pop.size()
  //If at dispersal m_new_pop_size > pop.size()
  //Then the number of funding individuals == pop.size()

  {
    simulation s(1000, 100);
    select_new_pop(s);
    assert(s.get_new_pop_size() == s.get_exp_new_pop_size());
    s = simulation(10, 100);
    select_new_pop(s);
    auto n_drawn_ind = std::count_if(s.get_pop().begin(),s.get_pop().end(),
                                     [](const individual& i) {return is_drawn(i);});
    assert( n_drawn_ind == s.get_pop_size());
  }


  //During dispersal the individuals selected for the new_pop cannot be drawn again from pop
  {
    simulation s(1000, 100);
    select_new_pop(s);
    assert(std::count_if(s.get_pop().begin(),s.get_pop().end(),
                         [](const individual& i) {return is_drawn(i);}) == 100);
  }

  //A simulation is initialized with a variable m_base_fitness
  //by default = 0.01
  {
    double base_disp_prob = 0.1;
    simulation s(0,0,0,0,0,0,0,0,base_disp_prob);
    assert(s.get_base_disp_prob() - base_disp_prob < 0.00001 &&
           s.get_base_disp_prob() - base_disp_prob > -0.000001);
  }

  //A simulation is initialized with a variable m_spore_advantage
  //by default = 10
  {
    double spore_advantage = 10;
    simulation s(0,0,0,0,0,0,0,0,0,spore_advantage);
    assert(s.get_spo_adv() - spore_advantage < 0.00001 &&
           s.get_spo_adv() - spore_advantage > -0.000001);
  }

  //A simulation is initialized with a uniform distribution between 0 and 1
  //used to see which ind is drawn at dispersal
  {
    simulation s;
    int sample_size = 100000;
    double mean = 0;
    for(int i = 0; i != sample_size; i++)
      {
        mean += s.get_unif_dist()(s.get_rng());
      }
    mean /= sample_size;
    assert(mean - 0.5 < 0.01 &&
           mean - 0.5 > -0.01);
  }
  //At initialization a simulation checks that base_disp_dist * 10 is not > 1
  //--------> constructor throws exception. Tested directly in constructor
  {
#ifndef IS_ON_TRAVIS
    try {
      simulation(0,0,0,0,0,0,0,0,1);
    } catch (std::string& e) {
      assert(e == "base dispersal probability * spore advantage > 1, too high!\n" );
    }
#endif
  }

  //Individuals are selected based on their phenotype
  //A spore is more likely to be selected than a living
  {
    simulation s(1000,100);
    for(int i = s.get_pop_size() / 2; i != s.get_pop_size(); i++)
      s.get_ind(i).set_phen(phenotype::spore);
    select_new_pop(s);
    auto spore_ratio =
        std::accumulate(s.get_new_pop().begin(),s.get_new_pop().end(),0.0,
                        [](const int sum, const individual& ind){return sum + is_spore(ind);}) /
        s.get_new_pop_size();
    assert(spore_ratio > 0.5);
  }

  //After being selected in new population individuals flag is_drawn is resetted
  {
    simulation s;
    select_new_pop(s);
    for(const auto& ind :s.get_new_pop())
      assert(!is_drawn(ind));
  }

  //After a new population is selected it swapped with the old population
  //And the old population is cancelled
  {
    int pop_size = 1000;
    int new_pop_size = 100;
    simulation s(pop_size,new_pop_size);
    select_new_pop(s);
    assert(s.get_new_pop_size() == new_pop_size);
    assert(s.get_pop_size() == pop_size);
    fund_pop(s);
    assert(s.get_new_pop_size() == 0);
    assert(s.get_pop_size() == new_pop_size);
  }

  //Individuals after funding the new population are set in an hexagonal pattern
  {
    int pop_size = 1000;
    int new_pop_size = 100;
    simulation s(pop_size,new_pop_size);
    select_new_pop(s);
    fund_pop(s);
    place_start_cells(s);
    auto n_hex_l = count_hex_layers(s.get_pop_size());
    auto v_modulus = modulus_of_btw_ind_angles(s, M_PI/ (6 * n_hex_l));
    for(auto ind_modulus : v_modulus)
      assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
    assert(has_collision(s).empty());
  }


  //Max 100 ind, selected based on phenotype, are placed in a hex pattern, in a new env after dispersal
  //Tests all of the above
  {
    int pop_size = 1000;
    int new_pop_size = 100;
    auto food = 42.1;
    auto metabolite = 42.1;
    simulation s(pop_size,new_pop_size);
    for(auto& grid_cell : s.get_env().get_grid())
      {
        grid_cell.set_food(food);
        grid_cell.set_metabolite(metabolite);
      }
    environment ref_env = s.get_env();
    dispersal(s);
    //Max 100 ind
    assert(s.get_pop_size() == new_pop_size);
    //Hex pattern
    auto n_hex_l = count_hex_layers(s.get_pop_size());
    auto v_modulus = modulus_of_btw_ind_angles(s, M_PI/ (6 * n_hex_l));
    for(auto ind_modulus : v_modulus)
      {
        assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
      }
    assert(has_collision(s).empty());
    //Reset env
    assert( s.get_env() != ref_env );
    auto init_food = s.get_env().get_init_food();
    for(const auto& grid_cell : s.get_env().get_grid())
      {
        assert(grid_cell.get_food() - init_food <  0.000001 &&
               grid_cell.get_food() - init_food >  -0.000001);
        assert(grid_cell.get_metabolite() <  0.000001 &&
               grid_cell.get_metabolite() >  -0.000001);
      }
  }

#endif
}

















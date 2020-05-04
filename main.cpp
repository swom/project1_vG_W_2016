#include <cassert>
#include <string>
#include "sim_view.h"

void test() {
  test_env_grid_cell();
  test_environment();
  test_grid_view();
  test_GRN();
  test_individual();
  test_individual_type();
  test_pop_param();
  test_population();
  test_simulation();
  test_sim_view();
  test_sim_parameters();
}


int main(int argc, char ** argv) //!OCLINT tests may be long
{

  const std::vector<std::string> args(argv, argv + argc);

#ifndef NDEBUG
  if (args.size() > 1 && args[1] == "--test")
    {
      test();
      // We've already tested, so the program is done
      return 0;
    }
  else
#else
  // In release mode, all asserts are removed from the code
  assert(1 == 2);
#endif
#ifndef LOGIC_ONLY
    if(args.size() > 1 && args[1] == "--visual")
      {
        simulation s(1,1,0.1,4,0.1,20 );
        sim_view v(s);

        int sim_time = 0;

        while (v.get_window().isOpen())
          {
            bool must_quit{v.process_events()};
            if (must_quit)
              return 0;
            if(sf::Keyboard::isKeyPressed(sf::Keyboard::A))
              {
                sim_time++;
                tick(v.get_sim());
              }
            if(sf::Keyboard::isKeyPressed(sf::Keyboard::R))
              {
               v.get_sim() = s;
               v.prepare_pop();
              }
            v.show();
          }
      }
#endif
  {
    simulation s(1,1,0.1,20,0.1,1);
    int time_now = 0;
    int time_limit = 200;
    while(time_now != time_limit)
      {
        tick(s);
        time_now++;
      }
  }

  return 0;

}


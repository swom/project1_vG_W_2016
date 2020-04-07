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
  test_simulation();
  test_sim_view();
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
#else
  // In release mode, all asserts are removed from the code
  assert(1 == 2);
#endif

#ifndef LOGIC_ONLY
  simulation s(7,0,0.1,30,0.01,20);

  if (args.size() > 1 && args[1] == "--visual")
    {
      sim_view v(s);

      while (v.get_window().isOpen())
        {
          bool must_quit{v.process_events()};
          if (must_quit)
            return 0;
          if(sf::Keyboard::isKeyPressed(sf::Keyboard::A))
            {
              tick(v.get_sim());
              v.show();
              continue;
            }
              v.show();
        }
    }
#endif

  return 0;
}


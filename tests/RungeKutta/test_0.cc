/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"

template class Sandals::RK4<2>;
template class Sandals::RK4<3>;
template class Sandals::RK4<4>;
template class Sandals::RK4<5>;
template class Sandals::RK4<6>;

int main(void) {
  std::cout << "Hi, I'm finding the sandal size for you!" << std::endl;
  return 0;
}

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

// Broyden's Ugly solver.
template class Sandals::BroydenUgly<2>;
template class Sandals::BroydenUgly<3>;
template class Sandals::BroydenUgly<4>;
template class Sandals::BroydenUgly<5>;
template class Sandals::BroydenUgly<6>;
template class Sandals::BroydenUgly<7>;
template class Sandals::BroydenUgly<8>;
template class Sandals::BroydenUgly<9>;

// Broyden's Bad solver.
template class Sandals::BroydenBad<2>;
template class Sandals::BroydenBad<3>;
template class Sandals::BroydenBad<4>;
template class Sandals::BroydenBad<5>;
template class Sandals::BroydenBad<6>;
template class Sandals::BroydenBad<7>;
template class Sandals::BroydenBad<8>;
template class Sandals::BroydenBad<9>;

// Broyden's Good solver.
template class Sandals::BroydenGood<2>;
template class Sandals::BroydenGood<3>;
template class Sandals::BroydenGood<4>;
template class Sandals::BroydenGood<5>;
template class Sandals::BroydenGood<6>;
template class Sandals::BroydenGood<7>;
template class Sandals::BroydenGood<8>;
template class Sandals::BroydenGood<9>;

// Broyden's Good solver.
template class Sandals::BroydenCombined<2>;
template class Sandals::BroydenCombined<3>;
template class Sandals::BroydenCombined<4>;
template class Sandals::BroydenCombined<5>;
template class Sandals::BroydenCombined<6>;
template class Sandals::BroydenCombined<7>;
template class Sandals::BroydenCombined<8>;
template class Sandals::BroydenCombined<9>;

// Newton solver.
template class Sandals::Newton<2>;
template class Sandals::Newton<3>;
template class Sandals::Newton<4>;
template class Sandals::Newton<5>;
template class Sandals::Newton<6>;
template class Sandals::Newton<7>;
template class Sandals::Newton<8>;
template class Sandals::Newton<9>;

int main(void) {
  std::cout << "Hi, I'm finding the sandal size for you!" << std::endl;
  return 0;
}

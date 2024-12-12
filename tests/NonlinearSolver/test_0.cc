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

using Sandals::vec2; using Sandals::mat2;
using Sandals::vec3; using Sandals::mat3;
using Sandals::vec4; using Sandals::mat4;
using Sandals::vec5; using Sandals::mat5;
using Sandals::vec6; using Sandals::mat6;
using Sandals::vec7; using Sandals::mat7;
using Sandals::vec8; using Sandals::mat8;
using Sandals::vec9; using Sandals::mat9;

// Broyden's Ugly solver.
template class Sandals::BroydenUgly<vec2, mat2>;
template class Sandals::BroydenUgly<vec3, mat3>;
template class Sandals::BroydenUgly<vec4, mat4>;
template class Sandals::BroydenUgly<vec5, mat5>;
template class Sandals::BroydenUgly<vec6, mat6>;
template class Sandals::BroydenUgly<vec7, mat7>;
template class Sandals::BroydenUgly<vec8, mat8>;
template class Sandals::BroydenUgly<vec9, mat9>;

// Broyden's Bad solver.
template class Sandals::BroydenBad<vec2, mat2>;
template class Sandals::BroydenBad<vec3, mat3>;
template class Sandals::BroydenBad<vec4, mat4>;
template class Sandals::BroydenBad<vec5, mat5>;
template class Sandals::BroydenBad<vec6, mat6>;
template class Sandals::BroydenBad<vec7, mat7>;
template class Sandals::BroydenBad<vec8, mat8>;
template class Sandals::BroydenBad<vec9, mat9>;

// Broyden's Good solver.
template class Sandals::BroydenGood<vec2, mat2>;
template class Sandals::BroydenGood<vec3, mat3>;
template class Sandals::BroydenGood<vec4, mat4>;
template class Sandals::BroydenGood<vec5, mat5>;
template class Sandals::BroydenGood<vec6, mat6>;
template class Sandals::BroydenGood<vec7, mat7>;
template class Sandals::BroydenGood<vec8, mat8>;
template class Sandals::BroydenGood<vec9, mat9>;

// Broyden's Good solver.
template class Sandals::BroydenCombined<vec2, mat2>;
template class Sandals::BroydenCombined<vec3, mat3>;
template class Sandals::BroydenCombined<vec4, mat4>;
template class Sandals::BroydenCombined<vec5, mat5>;
template class Sandals::BroydenCombined<vec6, mat6>;
template class Sandals::BroydenCombined<vec7, mat7>;
template class Sandals::BroydenCombined<vec8, mat8>;
template class Sandals::BroydenCombined<vec9, mat9>;

// Newton solver.
template class Sandals::Newton<vec2, mat2>;
template class Sandals::Newton<vec3, mat3>;
template class Sandals::Newton<vec4, mat4>;
template class Sandals::Newton<vec5, mat5>;
template class Sandals::Newton<vec6, mat6>;
template class Sandals::Newton<vec7, mat7>;
template class Sandals::Newton<vec8, mat8>;
template class Sandals::Newton<vec9, mat9>;

int main(void) {
  std::cout << "Hello World!" << std::endl;
  return 0;
}

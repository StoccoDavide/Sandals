#include "Sandals.hh"
#include "Sandals/RungeKutta/RK4.hh"

int main()
{
  Sandals::RK4<double, 1, 1> RK4_1_1;
  Sandals::RK4<double, 1, 2> RK4_1_2;
  Sandals::RK4<double, 2, 1> RK4_2_1;
  Sandals::RK4<double, 2, 2> RK4_2_2;
  return 0;
}

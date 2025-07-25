// vim:ts=2:et
//===========================================================================//
//                     "Tests/LinAlgTuplesTest.cpp":                         //
//===========================================================================//
#include "SpaceBallistics/Maths/LinAlgT.hpp"

int main()
{
  using namespace std;
  using namespace SpaceBallistics;

  auto  x = make_tuple(1.0_m, 2.0_sec,       3.0_kg);
  auto  c = 5.0_m;
  auto  y = make_tuple(1.0,   1.0_sec/1.0_m, 1.0_kg/1.0_m);

  auto  r = Add(x, c, y);
  cout << get<0>(r) << endl;
  cout << get<1>(r) << endl;
  cout << get<2>(r) << endl;
  return 0;
}

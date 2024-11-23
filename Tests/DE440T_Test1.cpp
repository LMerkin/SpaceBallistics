// vim:ts=2:et
//===========================================================================//
//                         "Tests/DE440T_Test1.cpp":                         //
//    Ecliptical (J2000.0) Co-Ords of the Sun and 9 Planes Simultaneously    //
//===========================================================================//
#include <cstdio>
#include <cstdlib>

// Asserts must be enabled for this test:
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "SpaceBallistics/PhysForces/DE440T.h"

int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  //-------------------------------------------------------------------------//
  // Test for All Bodies::                                                   //
  //-------------------------------------------------------------------------//
  if (argc != 2 && argc != 3 && argc != 4)
  {
    cerr << "PARAMS: TDB_From [TDB_To [TimeStep_days]]" << endl;
    return 1;
  }
  // Get the JD_TDB:
  TDB  from    {Time_day(atof(argv[1]))};
  TDB  to    = (argc >= 3) ?     TDB{Time_day(atof(argv[2]))} : from;
  Time step  = To_Time((argc >= 4) ? Time_day(atof(argv[3]))  : 1.0_day);

  if (from < DE440T::Bits::From || to > DE440T::Bits::To || from > to ||
      !IsPos(step))
  {
    cerr << "Invalid TDB Range or Step" << endl;
    return 1;
  }

  // Generate the Ephemerides for the 10 Major Planets (incl the Sun, with EMB
  // instead of Earth, PlChB instead of Pluto):
  for (TDB tdb = from; tdb <= to; tdb += step)
  {
    PosKVBarEcl poss[10];
    VelKVBarEcl vels[10];
    DE440T::GetPlanetsBarEclPVs(tdb, poss, vels);

    // Output:
    printf("%.2lf\n", tdb.GetJD().Magnitude());
    for (int i = 0; i < 10; ++i)
      printf("%d\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             i,
             poss[i].x().Magnitude(), poss[i].y()  .Magnitude(),
             poss[i].z().Magnitude(), LenK(poss[i]).Magnitude(),
             vels[i].x().Magnitude(), vels[i].y()  .Magnitude(),
             vels[i].z().Magnitude(), VelK(vels[i]).Magnitude());
  }
  // All Done!
  return 0;
}

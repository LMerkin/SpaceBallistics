// vim:ts=2:et
//===========================================================================//
//                      "Tests/LocationsTest.cpp":                           //
//===========================================================================//
#include "SpaceBallistics/CoOrds/GeoLocations.hpp"
#include <iostream>

int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  cout << "Plesetsk_43_3:"
       << "\n\tLong = " << Plesetsk_43_3.Longitude()
       << "\n\tPhi  = " << Plesetsk_43_3.Latitude ()
       << "\n\tX    = " << Plesetsk_43_3.PosV()[0]
       << "\n\tY    = " << Plesetsk_43_3.PosV()[1]
       << "\n\tZ    = " << Plesetsk_43_3.PosV()[2]
       << "\n\tRho  = " << Plesetsk_43_3.Rho ()
       << endl;
  cout << "Plesetsk_43_4:"
       << "\n\tLong = " << Plesetsk_43_4.Longitude()
       << "\n\tPhi  = " << Plesetsk_43_4.Latitude ()
       << "\n\tX    = " << Plesetsk_43_4.PosV()[0]
       << "\n\tY    = " << Plesetsk_43_4.PosV()[1]
       << "\n\tZ    = " << Plesetsk_43_4.PosV()[2]
       << "\n\tRho  = " << Plesetsk_43_4.Rho ()
       << endl;
  cout << "Baykonur_31_6:"
       << "\n\tLong = " << Baykonur_31_6.Longitude()
       << "\n\tPhi  = " << Baykonur_31_6.Latitude()
       << "\n\tX    = " << Baykonur_31_6.PosV()[0]
       << "\n\tY    = " << Baykonur_31_6.PosV()[1]
       << "\n\tZ    = " << Baykonur_31_6.PosV()[2]
       << "\n\tRho  = " << Baykonur_31_6.Rho ()
       << endl;
  cout << "Vostochny_1S:"
       << "\n\tLong = " << Vostochny_1S.Longitude()
       << "\n\tLat  = " << Vostochny_1S.Latitude()
       << "\n\tX    = " << Vostochny_1S.PosV()[0]
       << "\n\tY    = " << Vostochny_1S.PosV()[1]
       << "\n\tZ    = " << Vostochny_1S.PosV()[2]
       << "\n\tRho  = " << Vostochny_1S.Rho ()
       << endl;
  return 0;
}

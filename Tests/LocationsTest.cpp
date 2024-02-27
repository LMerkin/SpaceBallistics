// vim:ts=2:et
//===========================================================================//
//                      "Tests/LocationsTest.cpp":                           //
//===========================================================================//
#include "SpaceBallistics/Locations.hpp"
#include <iostream>

int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  cout << "Plesetsk_43_3:"
       << "\n\tLong = " << Plesetsk_43_3.Longitude()
       << "\n\tPhi  = " << Plesetsk_43_3.Latitude()
       << "\n\tX    = " << Plesetsk_43_3.GeoCX()
       << "\n\tY    = " << Plesetsk_43_3.GeoCY()
       << "\n\tZ    = " << Plesetsk_43_3.GeoCZ()
       << "\n\tRho  = " << Plesetsk_43_3.GeoCRho()
       << endl;
  cout << "Plesetsk_43_4:"
       << "\n\tLong = " << Plesetsk_43_4.Longitude()
       << "\n\tPhi  = " << Plesetsk_43_4.Latitude()
       << "\n\tX    = " << Plesetsk_43_4.GeoCX()
       << "\n\tY    = " << Plesetsk_43_4.GeoCY()
       << "\n\tZ    = " << Plesetsk_43_4.GeoCZ()
       << "\n\tRho  = " << Plesetsk_43_4.GeoCRho()
       << endl;
  cout << "Baykonur_31_6:"
       << "\n\tLong = " << Baykonur_31_6.Longitude()
       << "\n\tPhi  = " << Baykonur_31_6.Latitude()
       << "\n\tX    = " << Baykonur_31_6.GeoCX()
       << "\n\tY    = " << Baykonur_31_6.GeoCY()
       << "\n\tZ    = " << Baykonur_31_6.GeoCZ()
       << "\n\tRho  = " << Baykonur_31_6.GeoCRho()
       << endl;
  cout << "Vostochny_1S:"
       << "\n\tLong = " << Vostochny_1S.Longitude()
       << "\n\tLat  = " << Vostochny_1S.Latitude()
       << "\n\tX    = " << Vostochny_1S.GeoCX()
       << "\n\tY    = " << Vostochny_1S.GeoCY()
       << "\n\tZ    = " << Vostochny_1S.GeoCZ()
       << "\n\tRho  = " << Vostochny_1S.GeoCRho()
       << endl;
  return 0;
}

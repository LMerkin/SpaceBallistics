// vim:ts=2:et
//===========================================================================//
//                              "AzimuthTest.cpp":                           //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz2_LaunchPads.hpp"
#include <iostream>

int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  cout << "Pad_Plesetsk_43_3: Azimuth="
       << To_Angle_deg(Pad_Plesetsk_43_3.Azimuth()) << endl;

  cout << "Pad_Plesetsk_43_4: Azimuth="
       << To_Angle_deg(Pad_Plesetsk_43_4.Azimuth()) << endl;

  cout << "Pad_Baykonur_31_6: Azimuth="
       << To_Angle_deg(Pad_Baykonur_31_6.Azimuth()) << endl;

  cout << "Pad_Vostochny_1S : Azimuth="
       << To_Angle_deg(Pad_Vostochny_1S .Azimuth()) << endl;

  return 0;
}

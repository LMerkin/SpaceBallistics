// vim:ts=2:et
//===========================================================================//
//                       "Src/PhysForces/DE440T.cpp":                        //
//                   Functions (but not Data) for DE440T                     //
//===========================================================================//
#include "SpaceBallistics/PhysForces/DE440T.hpp"

namespace SpaceBallistics::DE440T
{
  //=========================================================================//
  // Force compilation of the following template functions:                  //
  //=========================================================================//
# ifdef  INST_BAR
# undef  INST_BAR
# endif
# define INST_BAR(EE, B) \
  template void GetPlanetBar##EE##PV<Body::B> \
    (TDB a_tdb, PosKVBar##EE<Body::B>* a_pos, VelKVBar##EE<Body::B>* a_vel);

  //-------------------------------------------------------------------------//
  // "GetPlanetBarEqPV":                                                     //
  //-------------------------------------------------------------------------//
  INST_BAR(Eq,  Sun)
  INST_BAR(Eq,  Mercury)
  INST_BAR(Eq,  Venus)
  // NB: "Earth" is specialised, no need to instantiate it explicitly...
  INST_BAR(Eq,  Mars)
  INST_BAR(Eq,  Jupiter)
  INST_BAR(Eq,  Saturn)
  INST_BAR(Eq,  Uranus)
  INST_BAR(Eq,  Neptune)
  INST_BAR(Eq,  PlChB)
  INST_BAR(Eq,  EMB)

  //-------------------------------------------------------------------------//
  // "GetPlanetBarEclPV":                                                    //
  //-------------------------------------------------------------------------//
  INST_BAR(Ecl, Sun)
  INST_BAR(Ecl, Mercury)
  INST_BAR(Ecl, Venus)
  // NB: "Earth" is specialised, no need to instantiate it explicitly...
  INST_BAR(Ecl, Mars)
  INST_BAR(Ecl, Jupiter)
  INST_BAR(Ecl, Saturn)
  INST_BAR(Ecl, Uranus)
  INST_BAR(Ecl, Neptune)
  INST_BAR(Ecl, PlChB)
  INST_BAR(Ecl, EMB)

# undef INST_BAR
}
// End namespace SpaceBallistics::DE440T

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
  //-------------------------------------------------------------------------//
  // "GetPlanetBarEqPV":                                                     //
  //-------------------------------------------------------------------------//
  template void GetPlanetBarEqPV<Body::Sun>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::Mercury>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::Venus>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  // NB: "Earth" is specialised, no need to instantiate it explicitly...

  template void GetPlanetBarEqPV<Body::Mars>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::Jupiter>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::Saturn>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::Uranus>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::Neptune>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::PlChB>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  template void GetPlanetBarEqPV<Body::EMB>
                (TDB a_tdb, PosKVBarEq* a_pos,  VelKVBarEq* a_vel);

  //-------------------------------------------------------------------------//
  // "GetPlanetBarEclPV":                                                    //
  //-------------------------------------------------------------------------//
  template void GetPlanetBarEclPV<Body::Sun>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::Mercury>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::Venus>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  // NB: "Earth" is specialised, no need to instantiate it explicitly...

  template void GetPlanetBarEclPV<Body::Mars>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::Jupiter>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::Saturn>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::Uranus>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::Neptune>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::PlChB>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);

  template void GetPlanetBarEclPV<Body::EMB>
                (TDB a_tdb, PosKVBarEcl* a_pos, VelKVBarEcl* a_vel);
}
// End namespace SpaceBallistics::DE440T

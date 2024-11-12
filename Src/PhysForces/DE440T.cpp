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
  // "GetPlanetBEqPv":                                                       //
  //-------------------------------------------------------------------------//
  template void GetPlanetBEqPV<Body::Sun>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::Mercury>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::Venus>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  // NB: "Earth" is specialised, no need to instantiate it explicitly...

  template void GetPlanetBEqPV<Body::Mars>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::Jupiter>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::Saturn>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::Uranus>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::Neptune>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::PlChB>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  template void GetPlanetBEqPV<Body::EMB>
                (TDB a_tdb, PosKVBEq* a_pos,  VelKVBEq* a_vel);

  //-------------------------------------------------------------------------//
  // "GetPlanetBEclPV":                                                      //
  //-------------------------------------------------------------------------//
  template void GetPlanetBEclPV<Body::Sun>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::Mercury>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::Venus>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  // NB: "Earth" is specialised, no need to instantiate it explicitly...

  template void GetPlanetBEclPV<Body::Mars>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::Jupiter>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::Saturn>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::Uranus>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::Neptune>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::PlChB>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);

  template void GetPlanetBEclPV<Body::EMB>
                (TDB a_tdb, PosKVBEcl* a_pos, VelKVBEcl* a_vel);
}
// End namespace SpaceBallistics::DE440T

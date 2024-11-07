// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/BodyData.hpp":                  //
//                Main Characteristics of Solar System Bodies                //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/PhysForces/DE440T.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "BodyData" Struct: Defined by Specialisations:                          //
  //=========================================================================//
  template<Body BodyName>
  struct BodyData;

  //-------------------------------------------------------------------------//
  // Earth:                                                                  //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Earth>
  {
    // Equatorial Radius (from WGS84, as well as Rp);
    // EGM2008 provides        6'378.1363_km:
    constexpr static LenK Re = 6'378.1370_km;

    // Polar      Radius:
    constexpr static LenK Rp = 6'356.752314245_km;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from EMG2008; DE440T uses a somewhat diff-
    // erent value); NOT including the Moon, but including the atmosphere:
    constexpr static GMK  K  = GMK(398600.4415);

    // XXX: The discrepancy between the EMG2008 and DE440T is rather large!
    static_assert(K.ApproxEquals(DE440T::K<Body::Earth>, 5e-8));

    // EGM2008 truncated:
    constexpr static int MaxSpherHarmDegreeAndOrder = 600;

    // Axial Rotation Angular Velocity Vector, for a given Epoch: TODO
  };

  //-------------------------------------------------------------------------//
  // Moon:                                                                   //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Moon>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 1'738.0_km;

    // Polar      Radius:
    constexpr static LenK Rp = 1'736.0_km;
    static_assert(Rp < Re);

    // Gravitational  Field  Constant (from the GRGM1200A Lunar Gravity Field
    // model of 2016; DE440T uses a slightly different value);
    constexpr static GMK  K  = GMK(4902.8001224453001);
    static_assert(K.ApproxEquals(DE440T::K<Body::Moon>, 1e-9));

    // GRGM1200A truncated:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 600;
  };

  //-------------------------------------------------------------------------//
  // Sun:                                                                    //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Sun>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 696'300.0_km;

    // Polar Radius: We disregard the (very small) flattening of the Sun:
    constexpr static LenK Rp = Re;

    // Gravitational Field Constant (from DE440T);
    constexpr static GMK  K  = DE440T::K<Body::Sun>;

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Mercury:                                                                //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Mercury>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 2'439.0007_km;

    // Polar Radius: We neglect the flattening:
    constexpr static LenK Rp = Re;

    // Gravitational Field Constant (from DE440T):
    constexpr static GMK  K  = DE440T::K<Body::Mercury>;

    // gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Venus:                                                                  //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Venus>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 6'051.8_km;

    // Polar Radius: Venus is a perfect sphere:
    constexpr static LenK Rp = Re;

    // Gravitational Field Constant (from DE440T, incl the atmosphere):
    constexpr static GMK  K  = DE440T::K<Body::Venus>;

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Mars:                                                                   //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Mars>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 3'396.2_km;

    // Polar      Radius:
    constexpr static LenK Rp = 3'376.2_km;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440T, incl the atmosphere and the
    // moons):
    constexpr static GMK  K  = DE440T::K<Body::Mars>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Jupiter:                                                                //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Jupiter>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 71'492.0_km;

    // Polar      Radius:
    constexpr static LenK Rp = 66'854.0_km;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Jupiter>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Jupiter:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Saturn:                                                                 //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Saturn>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 60'268.0_km;

    // Polar      Radius:
    constexpr static LenK Rp = 54'364.0_km;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Saturn>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Saturn:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Uranus:                                                                 //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Uranus>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 25'559.0_km;

    // Polar      Radius:
    constexpr static LenK Rp = 24'973.0_km;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Uranus>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Uranus:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Neptune:                                                                //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Neptune>
  {
    // Equatorial Radius:
    constexpr static LenK Re = 24'764.0_km;

    // Polar      Radius:
    constexpr static LenK Rp = 24'341.0_km;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Neptune>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Neptune:
    constexpr static int  MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Pluto:                                                                  //
  //-------------------------------------------------------------------------//
  // FIXME: No model for Pluto yet, as in DE440T we only got data for the PlChB
  // (Pluto-Charon System BaryCenter), not for Pluto itself...
}
// End namespace SpaceBallistics

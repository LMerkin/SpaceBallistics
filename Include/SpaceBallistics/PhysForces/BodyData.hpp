// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/BodyData.hpp":                  //
//                Main Characteristics of Solar System Bodies                //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"

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
    // EGM2008 provides       6'378'136.3_m:
    constexpr static Len Re = 6'378'137.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 6'356'752.314245_m;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from EGM2008; DE440 uses a slightly diff-
    // erent value: GM(398600.435507 * 1e9));
    // NOT including the Moon, but including the atmosphere:
    constexpr static GM  K  = GM(398600.4415 * 1e9);

    // EGM2008 truncated:
    constexpr static int MaxSpherHarmDegreeAndOrder = 600;
  };

  //-------------------------------------------------------------------------//
  // Moon:                                                                   //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Moon>
  {
    // Equatorial Radius:
    constexpr static Len Re = 1'738'000.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 1'736'000.0_m;
    static_assert(Rp < Re);

    // Gravitational  Field  Constant (from the GRGM1200A Lunar Gravity Field
    // model (2016));  DE440 provides GM(4902.800118 * 1e9):
    constexpr static GM  K  = GM(4902.8001224453001  * 1e9);

    // GRGM1200A truncated:
    constexpr static int MaxSpherHarmDegreeAndOrder = 600;
  };

  //-------------------------------------------------------------------------//
  // Sun:                                                                    //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Sun>
  {
    // Equatorial Radius:
    constexpr static Len Re = 696'300'000.0_m;

    // Polar Radius: We disregard the (very small) flattening of the Sun:
    constexpr static Len Rp = Re;

    // Gravitational Field Constant (from DE440);
    constexpr static GM  K  = GM(132712440041.279419 * 1e9);

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Mercury:                                                                //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Mercury>
  {
    // Equatorial Radius:
    constexpr static Len Re = 2'439'000.7_m;

    // Polar Radius: We neglect the flattening:
    constexpr static Len Rp = Re;

    // Gravitational Field Constant (from DE440):
    constexpr static GM  K  = GM(22031.868551 * 1e9);

    // gravitational field is assumed to be spherically-symmetric:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Venus:                                                                  //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Venus>
  {
    // Equatorial Radius:
    constexpr static Len Re = 6'051'800.0_m;

    // Polar Radius: Venus is a perfect sphere:
    constexpr static Len Rp = Re;

    // Gravitational Field Constant (from DE440, incl the atmosphere):
    constexpr static GM  K  = GM(324858.592 * 1e9);

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Mars:                                                                   //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Mars>
  {
    // Equatorial Radius:
    constexpr static Len Re = 3'396'200.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 3'376'200.0_m;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440, incl the atmosphere and the
    // moons):
    constexpr static GM  K  = GM(42828.375816 * 1e9);

    // FIXME: The gravitational field is assumed to be spherically-symmetric:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Jupiter:                                                                //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Jupiter>
  {
    // Equatorial Radius:
    constexpr static Len Re = 71'492'000.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 66'854'000.0_m;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440, for the whole system, incl the
    // atmospehere and the moons):
    constexpr static GM  K  = GM(126712764.1 * 1e9);

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is grossly-incorrect in case of Jupiter:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Saturn:                                                                 //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Saturn>
  {
    // Equatorial Radius:
    constexpr static Len Re = 60'268'000.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 54'364'000.0_m;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440, for the whole system, incl the
    // atmospehere and the moons):
    constexpr static GM  K  = GM(37940584.8418 * 1e9);

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is grossly-incorrect in case of Saturn:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Uranus:                                                                 //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Uranus>
  {
    // Equatorial Radius:
    constexpr static Len Re = 25'559'000.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 24'973'000.0_m;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440, for the whole system, incl the
    // atmospehere and the moons):
    constexpr static GM  K  = GM(5794556.4 * 1e9);

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is grossly-incorrect in case of Uranus:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Neptune:                                                                //
  //-------------------------------------------------------------------------//
  template<>
  struct BodyData<Body::Neptune>
  {
    // Equatorial Radius:
    constexpr static Len Re = 24'764'000.0_m;

    // Polar      Radius:
    constexpr static Len Rp = 24'341'000.0_m;
    static_assert(Rp < Re);

    // Gravitational Field Constant (from DE440, for the whole system, incl the
    // atmospehere and the moons):
    constexpr static GM  K  = GM(6836527.10058 * 1e9);

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is grossly-incorrect in case of Neptune:
    constexpr static int MaxSpherHarmDegreeAndOrder = 0;
  };
}
// End namespace SpaceBallistics

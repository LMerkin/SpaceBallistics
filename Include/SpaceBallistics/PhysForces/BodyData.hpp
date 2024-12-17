// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/BodyData.hpp":                  //
//               Physical Characteristics of Solar System Bodies             //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"

namespace SpaceBallistics
{
  namespace DE440T
  {
    //=======================================================================//
    // Gravitational and Mass Params of Bodies from DE440T:                  //
    //=======================================================================//
    template<Body B>
    constexpr GMK K;

    template<> constexpr inline GMK K<Body::Sun>     = GMK(132712440041.279419);
    template<> constexpr inline GMK K<Body::Mercury> = GMK(       22031.868551);
    template<> constexpr inline GMK K<Body::Venus>   = GMK(      324858.592   );
    template<> constexpr inline GMK K<Body::Earth>   = GMK(      398600.435507);
    template<> constexpr inline GMK K<Body::Moon>    = GMK(        4902.800118);
    template<> constexpr inline GMK K<Body::EMB>     = K<Body::Earth>
                                                     + K<Body::Moon>;
    template<> constexpr inline GMK K<Body::Mars>    = GMK(       42828.375816);
    template<> constexpr inline GMK K<Body::Jupiter> = GMK(   126712764.1     );
    template<> constexpr inline GMK K<Body::Saturn>  = GMK(    37940584.8418  );
    template<> constexpr inline GMK K<Body::Uranus>  = GMK(     5794556.4     );
    template<> constexpr inline GMK K<Body::Neptune> = GMK(     6836527.10058 );
    template<> constexpr inline GMK K<Body::PlChB>   = GMK(         975.5     );

    // Earth/Moon Mass Ratio:
    constexpr  inline double EMRat = 81.3005682214972154;
    static_assert(K<Body::Earth>.ApproxEquals(K<Body::Moon> * EMRat, 1e-10));
  }
  // End namespace "DE440T"

  //=========================================================================//
  // Gravitational Field Model Coeffs for a Body:                            //
  //=========================================================================//
  // The struct for Dimension-Less Spherical Harmonics Coeffs representing the
  // Gravitational Potential, with Geodesy-style normalisation (to 4*Pi):
  //
  struct SpherHCoeffs
  {
    // Data Flds (Plain Vanilla):
    int    m_l;      // 2  .. MaxDeg
    int    m_m;      // 0  .. m_l
    double m_Clm;    // Coeff at Cos
    double m_Slm;    // Coeff at Sin
  };

  //=========================================================================//
  // "BodyData" Struct: Defined by Specialisations:                          //
  //=========================================================================//
  template<Body B>
  struct   BodyData;

  //-------------------------------------------------------------------------//
  // "UNDEFINED":                                                            //
  //-------------------------------------------------------------------------//
  // Assumed to be a point of negligible size and mass:
  //
  template<>
  struct BodyData<Body::UNDEFINED>
  {
    // Equatorial and Polar Radius:
    constexpr static LenK Re = 0.0_km;
    constexpr static LenK Rp = 0.0_km;

    // The Gravitational Field Constant: 0, since Mass=0:
    constexpr static GMK  K  = GMK(0.0);

    // And obviously no harmonics in the Gravitational Fld expansion:
    constexpr static int MaxSpherHDegreeAndOrder = 0;
  };

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

    // The Gravitational Field Constant (from EMG2008; DE440T uses a somewhat
    // different value); NOT including the Moon, but including the atmosphere:
    constexpr static GMK  K  = GMK(398600.4415);

    // XXX: The discrepancy between the EMG2008 and DE440T is rather large!
    static_assert(K.ApproxEquals(DE440T::K<Body::Earth>, 5e-8));

    // EGM2008 truncated:
    constexpr static int MaxSpherHDegreeAndOrder = 600;

    // The Actual Gravitational Potential Coeffs:
    static SpherHCoeffs const GravFldModelCoeffs
      [ ((MaxSpherHDegreeAndOrder + 1) * (MaxSpherHDegreeAndOrder + 2)) / 2 ];
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

    // The Gravitational Field Constant (from the GRGM1200A Lunar Gravity Field
    // model of 2016; DE440T uses a slightly different value);
    constexpr static GMK  K  = GMK(4902.8001224453001);
    static_assert(K.ApproxEquals(DE440T::K<Body::Moon>, 1e-9));

    // GRGM1200A truncated:
    constexpr static int  MaxSpherHDegreeAndOrder = 600;

    // The Actual Gravitational Potential Coeffs:
    static SpherHCoeffs const GravFldModelCoeffs
      [ ((MaxSpherHDegreeAndOrder + 1) * (MaxSpherHDegreeAndOrder + 2)) / 2 ];
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

    // The Gravitational Field Constant (from DE440T):
    constexpr static GMK  K  = DE440T::K<Body::Sun>;

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T):
    constexpr static GMK  K  = DE440T::K<Body::Mercury>;

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T, incl the atmosphere):
    constexpr static GMK  K  = DE440T::K<Body::Venus>;

    // The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T, incl the atmosphere and
    // the moons):
    constexpr static GMK  K  = DE440T::K<Body::Mars>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Jupiter>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Jupiter:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Saturn>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Saturn:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Uranus>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Uranus:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
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

    // The Gravitational Field Constant (from DE440T, for the whole system, incl
    // the atmospehere and the moons):
    constexpr static GMK  K  = DE440T::K<Body::Neptune>;

    // FIXME: The gravitational field is assumed to be spherically-symmetric.
    // This is GROSSLY-INCORRECT in case of Neptune:
    constexpr static int  MaxSpherHDegreeAndOrder = 0;
  };

  //-------------------------------------------------------------------------//
  // Pluto:                                                                  //
  //-------------------------------------------------------------------------//
  // FIXME: No model for Pluto yet, as in DE440T we only got data for the PlChB
  // (Pluto-Charon System BaryCenter), not for Pluto itself...
}
// End namespace SpaceBallistics

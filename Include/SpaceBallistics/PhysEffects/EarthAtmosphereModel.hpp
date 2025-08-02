// vim:ts=2:et
//===========================================================================//
//          "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp":          //
//===========================================================================//
// NB: This model is the International Standard Atmosphere (ISA),  valid up to
// the altitude of ~85 km. TODO: Extensions for higher altitudes (eg for pred-
// ictions of satellite orbits decay):
//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include <tuple>

namespace SpaceBallistics::EarthAtmosphereModel
{
  //=========================================================================//
  // Consts:                                                                 //
  //=========================================================================//
  // Standard Atmospheric Pressure at Sea Level:
  constexpr inline Pressure P0       = Pressure(101325.0);

  // Standard Atmospheric Temperature at Sea Level: 15 deg C:
  constexpr inline AbsTemp  T0       = 288.15_K;

  // Standard Atmospheric Density  at Sea Level (@ T0):
  // XXX: In SU2 Standard Air model:           1.2172 -- why?
  constexpr inline Density  Rho0     = Density(1.225);

  // Specific Gas Constant for Dry Air (J/K):
  // XXX: In SU2 Standard Air model:   287.058 -- why?
  constexpr inline auto     RAir     = 287.0528 * Sqr(Vel(1.0)) / 1.0_K;

  // cP/cV Ratio for the Dry Air:
  constexpr inline double   GammaAir = 1.4;

  constexpr inline auto g0K          = To_Len_km(g0);
  constexpr inline auto RAirK        = To_Len_km(RAir);
  using     TempGrad                 = decltype(1.0_K / 1.0_km);

  //=========================================================================//
  // "LayerInfo" Struct:                                                     //
  //=========================================================================//
  // NB: Here all altitudes are GeoPotential ones above the Mean Sea Level.
  // XXX:
  // (*) Interestingly, the functions below use "g0K" rather than variable
  //     (with the altitude) acceleration of gravity:
  // (*) In CLang (major version <= 20),  C++23 "constexpr" mathematical funct-
  //     ions are not yet implemented, so some methods of this struct are decl-
  //     ared "constexpr" in GCC but not in CLang:
  //
  struct LayerInfo
  {
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    LenK     const m_baseH;
    Pressure const m_baseP;
    AbsTemp  const m_baseT;
    TempGrad const m_tempGrad;
    LenK     const m_endH;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    constexpr LayerInfo
    (
      LenK         a_base_h,
      Pressure     a_base_p,
      AbsTemp      a_base_T,
      TempGrad     a_temp_grad,
      LenK         a_end_h
    )
    : m_baseH     (a_base_h),
      m_baseP     (a_base_p),
      m_baseT     (a_base_T),
      m_tempGrad  (a_temp_grad),
      m_endH      (a_end_h)
    {
      assert(m_baseH < m_endH);
    }

    // Default Ctor is auto-generated (XXX: for some reason, we need it...):
    constexpr LayerInfo() = default;

    //-----------------------------------------------------------------------//
    // Pressure within this Layer:                                           //
    //-----------------------------------------------------------------------//
#   ifndef __clang_
    constexpr
#   endif
    Pressure P(LenK a_h) const
    {
      assert(m_baseH <= a_h && a_h <= m_endH);
      return
        (!IsZero(m_tempGrad))
        ? m_baseP *
          Exp(double(-g0K / (RAirK * m_tempGrad))   *
              Log(1.0 + double(m_tempGrad / m_baseT * (a_h - m_baseH))))

        : m_baseP *
          Exp(double(-g0K / (RAirK * m_baseT)       * (a_h - m_baseH)));
    }

    //-----------------------------------------------------------------------//
    // Temperature within this Layer:                                        //
    //-----------------------------------------------------------------------//
    constexpr AbsTemp T(LenK a_h) const
    {
      assert(m_baseH <= a_h && a_h <= m_endH);
      return m_baseT + m_tempGrad * (a_h - m_baseH);
    }

    //-----------------------------------------------------------------------//
    // Constructing the "next" ("above") Layer from this one:                //
    //-----------------------------------------------------------------------//
#   ifndef __clang__
    constexpr
#   endif
    LayerInfo MkNextLayer
    (
      TempGrad  a_next_temp_grad,
      LenK      a_next_end_h
    )
    const
    {
      assert(a_next_end_h > m_endH);
      return LayerInfo
             (m_endH, P(m_endH), T(m_endH), a_next_temp_grad, a_next_end_h);
    }
  };

  //=========================================================================//
  // All Layers (up to ~85 km):                                              //
  //=========================================================================//
# ifndef __clang__
  constexpr inline LayerInfo Layers[7]
  {
# include "SpaceBallistics/PhysEffects/EarthAtmosphereLayers.h"
  };
# else
  // In CLang (major ver <= 20), we cannot create "constexpr" Layers:
  extern LayerInfo const Layers[7];
# endif

  //=========================================================================//
  // Pressure and Temperature for Any Altitide:                              //
  //=========================================================================//
  // NB: here "a_z" is a Geometric Altitude, NOT The GeoPotential one.
  // Returns (P, Rho, T, SpeedOfSound):
# ifndef __clang__
  constexpr
# else
  inline
# endif
  std::tuple<Pressure, Density, AbsTemp, Vel> AirParams(LenK a_z)
  {
    // FIXME: For the moment, altitudes below the MSL are not allowed:
    assert(!IsNeg(a_z));

    // First, compute the GeoPotential altitude from the Geometric one:
    LenK h = BodyData<Body::Earth>::Rp * a_z /
            (BodyData<Body::Earth>::Rp + a_z);
    assert(h <= a_z);

    for (LayerInfo const& l: Layers)
      if (l.m_baseH <= h && h <= l.m_endH)
      {
        // Found the Layer which "h" belongs to:
        Pressure p   =  l.P(h);
        AbsTemp  T   =  l.T(h);
        assert(IsPos(p) &&  IsPos(T));
        Density  rho = p / (RAir * T);
        Vel      a   = SqRt(GammaAir * RAir * T);
        return std::make_tuple(p, rho, T, a);
      }
    // If we got here: We must be above all Layers:
    return std::make_tuple
           (Pressure(0.0), Density(0.0), AbsTemp(0.0), Vel(0.0));
  }
}
// End namespace SpaceBallistics::EarthAtmosphereModel

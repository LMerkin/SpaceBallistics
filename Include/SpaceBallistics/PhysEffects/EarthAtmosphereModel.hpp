// vim:ts=2:et
//===========================================================================//
//          "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp":          //
//===========================================================================//
// NB: This model is the Russian GOST 4401-81 (which coincides with the Inter-
// national Standard Atmosphere up to the altitude of 85 km).  XXX: It is not
// detailed enough to predict satellite orbits decay; for that, NRLMSISE-00 or
// JB2008 should be used (TODO):
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

  // Specific Gas Constant for Dry Air (J/K), assuming the constant Molar Mass
  // (0.0028964420 kg/mol) for the altitudes up to 94 km; above that,  the air
  // composition and molar mass are not invariant anymore:
  // XXX: In SU2 Standard Air model:   287.058 -- why?
  constexpr inline auto     RAir     = 287.05287 * Sqr(Vel(1.0)) / 1.0_K;

  // cP/cV Ratio for the Dry Air:
  constexpr inline double   GammaAir = 1.4;
  constexpr inline auto RAirK        = To_Len_km(RAir);
  using     TempGrad                 = decltype(1.0_K / 1.0_km);

  //=========================================================================//
  // "LowLayerInfo" Struct:                                                  //
  //=========================================================================//
  struct LowLayerInfo
  {
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    LenK     const m_baseH;     // GEOPOTENTIAL Altitude ("h")
    Pressure const m_baseP;
    AbsTemp  const m_baseT;
    TempGrad const m_tempGrad;
    LenK     const m_endH;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    constexpr LowLayerInfo
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
    constexpr LowLayerInfo() = default;

    //-----------------------------------------------------------------------//
    // Pressure within this Layer:                                           //
    //-----------------------------------------------------------------------//
    // Non-"constexpr" because of CLang issues:
    //
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
    // Non-"constexpr" because of CLang issues:
    //
    LowLayerInfo MkNextLayer
    (
      TempGrad  a_next_temp_grad,
      LenK      a_next_end_h
    )
    const
    {
      assert(a_next_end_h > m_endH);
      return LowLayerInfo
             (m_endH, P(m_endH), T(m_endH), a_next_temp_grad, a_next_end_h);
    }
  };

  //=========================================================================//
  // "UpperLayerInfo":                                                       //
  //=========================================================================//
  struct UpperLayerInfo
  {
    LenK     const m_baseZ;   // GEOMETRIC Altitude ("z")
    AbsTemp  const m_baseT;
    Pressure const m_baseP;
    Density  const m_baseRho;
  };

  //=========================================================================//
  // Info Tables:                                                            //
  //=========================================================================//
  // "Low"     Layers (z = 0 .. 93 km):
  constexpr inline int  NLowLayers    = 8;
  extern LowLayerInfo   const     LowLayers   [NLowLayers];

  // "Upper1"  Layers (z = 93+ .. 300   km, step = 1 km):
  constexpr inline int  NUpperLayers1 = 208;
  extern UpperLayerInfo const     UpperLayers1[NUpperLayers1];

  // "Upper2"  Layers (z = 300+ .. 500  km, step = 2 km):
  constexpr inline int  NUpperLayers2 = 101;
  extern UpperLayerInfo const     UpperLayers2[NUpperLayers2];

  // "Upper5"  Layers (z = 500+ .. 1200 km, step = 5 km):
  constexpr inline int  NUpperLayers5 = 141;
  extern UpperLayerInfo const     UpperLayers5[NUpperLayers5];

  //=========================================================================//
  // Atmospheric Conditions for Any Altitide:                                //
  //=========================================================================//
  using AtmConds = std::tuple<Pressure, Density, AbsTemp, Vel>;

  // NB: here "a_z" is a GEOMETRIC Altitude, NOT The GeoPotential one;
  // Returns (P, Rho, T, SpeedOfSound):
  //
  AtmConds GetAtmConds(LenK a_z);
}
// End namespace SpaceBallistics::EarthAtmosphereModel

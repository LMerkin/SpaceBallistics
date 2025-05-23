// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/Locations.h":                  //
//                       Locations on the Body Surfaces                      //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/AngleUtils.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include <utility>
#include <tuple>

namespace SpaceBallistics
{
  //=========================================================================//
  // "BodyLocations" Class:                                                  //
  //=========================================================================//
  // XXX: The Body is assumed to be an Ellipsoid of Rotation:
  //
  template<Body BBody>
  class Location
  {
    //-----------------------------------------------------------------------//
    // Consts: Parameters of the Ellipsoidal Body:
    //-----------------------------------------------------------------------//
  public:
    // Equatorial Radius:
    constexpr static LenK Re = BodyData<BBody>::Re;

    // Polar Radius:
    constexpr static LenK Rp = BodyData<BBody>::Rp;
    static_assert(Rp <= Re);

    // Flattening:
    constexpr static double FlatC = double(Rp / Re);
    constexpr static double Flat  = 1.0 - FlatC;

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Primary Body-detic (~Body-graphic) Co-Ords:
    Angle           m_lambda;   // Longitide (rad)
    Angle           m_phi;      // Latitude  (rad)
    Len             m_h;        // Elevation (m)

    // Derived Rectangular Co-Ords (in the "BodyCRotatingCOS"):
    PosKVRot<BBody> m_r;        // (x, y, z)
    LenK            m_rho;      // Radius-vector from Body center

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Longitude and Latitude are in Degs (fractional):
    //
    constexpr Location(Angle_deg a_lambda,  Angle_deg a_phi, Len a_h)
    : Location        (To_Angle (a_lambda), To_Angle (a_phi),    a_h)
    {}

    // Longitude and Latitude are in ('+'|' '|'-', Deg, Min, Sec):
    constexpr Location
    (
      std::tuple<char,int,int,double>  a_lambda,
      std::tuple<char,int,int,double>  a_phi,
      Len                              a_h
    )
    : Location(ToRads(a_lambda), ToRads(a_phi), a_h)
    {}

    // Longitude and Latitude are in Rads:
    constexpr Location(Angle a_lambda,  Angle a_phi,  Len a_h)
    : m_lambda(a_lambda),
      m_phi   (a_phi),
      m_h     (a_h)
    {
      // Checks:
      static_assert(IsPos(Rp) && Re >= Rp && Flat >= 0.0);

      // Derived Rectangular Co-Ords (in the "BodyCRotatingCOS"):
      // a = Re, b = Rp:
      // XXX: We do NOT apply any simplifications in the case Flat==0 here:
      auto   a2  = Sqr(Re);
      LenK   bt  = Rp * Tan(m_phi);
      auto   bt2 = Sqr(bt);
      auto   bbt = Rp * bt;
      LenK   d   = SqRt(a2 + bt2);
      // (x,z) in the cross-section through the axis and the given point:
      LenK   x   = a2  / d;
      LenK   z   = bbt / d;
      m_r[0]     = x * Cos(m_lambda);            // BodyC x
      m_r[1]     = x * Sin(m_lambda);            // BodyC y
      m_r[2]     = z;                            // BodyC z
      m_rho      = SqRt(Sqr(a2) + Sqr(bbt)) / d; // BodyC r
    }

    //-----------------------------------------------------------------------//
    // Other Ctors:                                                          //
    //-----------------------------------------------------------------------//
    // There is no Default Ctor:
    Location() = delete;

    // The CopyCtor is available:
    constexpr Location(Location const&) = default;

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    // Body-detic Co-Ords:
    constexpr Angle                  Longitude() const { return m_lambda; }
    constexpr Angle                  Latitude () const { return m_phi;    }
    constexpr Len                    Elevation() const { return m_h;      }

    // BodyCentric Rectangualar Co-Ords (in the "BodyCRotatingCOS"):
    // The object itself characterised by a "Location" is NOT a Body, hence
    // "Body::UNDEFINED" in its "PosKVRot":
    constexpr PosKVRot<BBody> const& PosKV    () const { return m_r;      }
    constexpr LenK                   RhoK     () const { return m_rho;    }

    //-----------------------------------------------------------------------//
    // Util: Azimuth(degs) computation from a Tangential Vector:             //
    //-----------------------------------------------------------------------//
    constexpr static Angle_deg GetAzimuth
    (
      Angle_deg a_from_lambda,  // From: (Longitude, Latitude)
      Angle_deg a_from_phi,     //
      Angle_deg a_to_lambda,    // To  : (Longitude, Latitude)
      Angle_deg a_to_phi
    )
    {
      // The vector should be sufficiently short, otherwise the approximations
      // used below may not be valid:
      Angle  dLambda = To_Angle(a_to_lambda - a_from_lambda);
      Angle  dPhi    = To_Angle(a_to_phi    - a_from_phi);
      assert(Abs(dLambda) < 0.001_rad && Abs(dPhi) < 0.001_rad);
      Angle  Tol(DefaultTol<double>);

      if (Abs(dPhi) > Tol)
      {
        // dPhi != 0, ie the Azimuth is not (+-Pi/2):
        Angle  avgPhi = To_Angle (a_to_phi + a_from_phi) / 2.0;
        double tgA    =
          double(dLambda /
                (dPhi * FlatC * SqRt(1.0 + Sqr(FlatC * Tan(avgPhi)))));
        Angle  A(ATan(tgA));

        // "A" is in (-Pi/2 .. +Pi/2), modify it to [0..2*Pi):
        if (IsNeg(dPhi))
        {
          assert(IsNeg(A));
          A += PI;
        }
        else
        if (IsNeg(dLambda))
        {
          assert(IsNeg(A));
          A += TWO_PI;
        }
        // We should get 0 <= A < 2*Pi, but cannot formally assert that due to
        // possible rounding errors:
        return To_Angle_deg(A);
      }
      else
      if (dLambda >  Tol)
        // Due East:
        return 90.0_deg;
      else
      if (dLambda < -Tol)
        // Due West:
        return 270.0_deg;
      else
        // Both dLambda and dPhi are ~0, so the Azimuth is undefined:
        return Angle_deg(NaN<double>);
    }
  };
}
// End namespace SpaceBallistics

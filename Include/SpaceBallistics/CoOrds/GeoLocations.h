// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/GeoLocations.h":                //
//                              Geodetic Locations                           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/Utils.hpp"
#include "SpaceBallistics/CoOrds/GeoCentricRotatingCOS.h"
#include <tuple>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Location_WGS84" Class:                                                 //
  //=========================================================================//
  // Geodetic Location in WGS84 System:
  //
  class Location_WGS84
  {
  public:
    //-----------------------------------------------------------------------//
    // Consts:                                                               //
    //-----------------------------------------------------------------------//
    // Semi-Major Axis of the Earth Ellipsoid ("Equatorial Radius"):
    constexpr static Len    a   = 6'378'137.0_m;

    // Semi-Minor Axis of the Earth Ellipsoid ("Polar Radius"):
    constexpr static Len    b   = 6'356'752.314245_m;

    // Flattening of the Earth:
    constexpr static double fc  = double(b / a);
    constexpr static double f   = 1.0 - fc;

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Primary Geodetic Co-Ords:
    Angle   m_lambda;   // Longitide (rad)
    Angle   m_phi;      // Latitude  (rad)
    Len     m_h;        // Elevation (m)

    // Derived Geocentric Rectangular Co-Ords (in the GeoCentricRotatingCOS):
    PosVGR  m_r;        // (x, y, z)
    Len     m_rho;      // Radius-vector from Earth center

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Longitude and Latitude are in Degs (fractional):
    //
    constexpr Location_WGS84(Angle_deg a_lambda,  Angle_deg    a_phi, Len a_h)
    : Location_WGS84     (To_Angle_rad(a_lambda), To_Angle_rad(a_phi),    a_h)
    {}

    // Longitude and Latitude are in ('+'|' '|'-', Deg, Min, Sec):
    constexpr Location_WGS84
    (
      std::tuple<char,int,int,double>  a_lambda,
      std::tuple<char,int,int,double>  a_phi,
      Len                              a_h
    )
    : Location_WGS84(ToRads(a_lambda), ToRads(a_phi), a_h)
    {}

    // Longitude and Latitude are in Rads:
    constexpr Location_WGS84(Angle a_lambda,  Angle a_phi,  Len a_h)
    : m_lambda(a_lambda),
      m_phi   (a_phi),
      m_h     (a_h)
    {
      // Derived Geocentric Rectangular Co-Ords:
      Area   a2  = Sqr(a);
      Len    bt  = b * Tan(double(m_phi));
      Area   bt2 = Sqr(bt);
      Area   bbt = b * bt;
      Len    d   = SqRt(a2 + bt2);
      // (x,z) in the cross-section through the axis and the given point:
      Len    x   = a2  / d;
      Len    z   = bbt / d;
      m_r[0]     = x * Cos(double(m_lambda));    // GeoCentric x
      m_r[1]     = x * Sin(double(m_lambda));    // GeoCentric y
      m_r[2]     = z;                            // GeoCentric z
      m_rho      = SqRt(Sqr(a2) + Sqr(bbt)) / d; // GeoCentric Radius-Vector
    }

    //-----------------------------------------------------------------------//
    // Other Ctors:                                                          //
    //-----------------------------------------------------------------------//
    // There is no Default Ctor:
    Location_WGS84() = delete;

    // The CopyCtor is available:
    constexpr Location_WGS84(Location_WGS84 const&) = default;

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    // Geodetic Co-Ords:
    constexpr Angle         Longitude()     const { return m_lambda; }
    constexpr Angle         Latitude ()     const { return m_phi;    }
    constexpr Len           Elevation()     const { return m_h;      }

    // GeoCentric Rectangualar Co-Ords (in the GeoCentricRotatingCOS):
    constexpr PosVGR const& GeoCentricPos() const { return m_r;      }
    constexpr Len           GeoCentricR  () const { return m_rho;    }

    //-----------------------------------------------------------------------//
    // Util: Azimuth(degs) computation from a Terrestrial Vector:            //
    //-----------------------------------------------------------------------//
    constexpr inline Angle_deg GetAzimuthDegs
    (
      Angle_deg a_from_lambda,  // From: (Longitude, Latitude)
      Angle_deg a_from_phi,     //
      Angle_deg a_to_lambda,    // To  : (Longitude, Latitude)
      Angle_deg a_to_phi
    )
    {
      // The vector should be sufficiently short, otherwise the approximations
      // used below may not be valid:
      double dLambda = double(To_Angle_rad(a_to_lambda - a_from_lambda));
      double dPhi    = double(To_Angle_rad(a_to_phi    - a_from_phi));
      assert(Abs(dLambda) < 1e-3 && Abs(dPhi) < 1e-3);

      if (Abs(dPhi) > Tol)
      {
        // dPhi != 0, ie the Azimuth is not (+-Pi/2):
        double avgPhi = double (To_Angle_rad(a_to_phi + a_from_phi)) / 2.0;
        double tgA    = dLambda / (dPhi * fc * SqRt(1 + Sqr(fc * Tan(avgPhi))));

        // XXX: std::atan() might not be "constexpr" in CLang yet; may need to
        // implement our own "constexpr" ATan function:
        Angle_deg A   = To_Angle_deg(Angle_rad(std::atan(tgA)));

        // "A" is in (-90 .. +90), modify it to [0..360):
    }
  };
}

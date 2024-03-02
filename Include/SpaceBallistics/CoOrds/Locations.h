// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/Locations.h":                  //
//                              Geodetic Locations                           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/Utils.hpp"
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
    constexpr static Len a = 6'378'137.0_m;

    // Semi-Minor Axis of the Earth Ellipsoid ("Polar Radius"):
    constexpr static Len b = 6'356'752.314245_m;

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Primary Geodetic Co-Ords:
    double m_lambda;     // Longitide (rad)
    double m_phi;        // Latitude  (rad)
    Len    m_h;          // Elevation (m)

    // Derived Geocentric Rectangular Co-Ords (in GeoCRotatingCOS):
    Len    m_r[3];       // (x, y, z)
    Len    m_rho;        // Radius-vector from Earth center

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Longitude and Latitude are in ('+'|' '|'-', Deg, Min, Sec):
    //
    constexpr Location_WGS84
    (
      std::tuple<char,int,int,double>  a_lambda,
      std::tuple<char,int,int,double>  a_phi,
      Len                              a_h
    )
    : // Primary Geodetic Co-Ords:
      m_lambda(ToRads(a_lambda)),
      m_phi   (ToRads(a_phi)),
      m_h     (a_h)
    {
      // Derived Geocentric Rectangular Co-Ords:
      Area   a2  = Sqr(a);
      Len    bt  = b * Tan(m_phi);
      Area   bt2 = Sqr(bt);
      Area   bbt = b * bt;
      Len    d   = SqRt(a2 + bt2);
      // (x,z) in the cross-section through the axis and the given point:
      Len    x   = a2  / d;
      Len    z   = bbt / d;
      m_r[0]     = x * Cos(m_lambda);            // GeoCentric x
      m_r[1]     = x * Sin(m_lambda);            // GeoCentric y
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
    constexpr double Longitude() const { return m_lambda; }
    constexpr double Latitude () const { return m_phi;    }
    constexpr Len    Elevation() const { return m_h;      }

    // GeoCentric Rectangualar Co-Ords:
    constexpr Len    GeoCX    () const { return m_r[0];   }
    constexpr Len    GeoCY    () const { return m_r[1];   }
    constexpr Len    GeoCZ    () const { return m_r[2];   }
    constexpr Len    GeoCRho  () const { return m_rho;    }
  };
}

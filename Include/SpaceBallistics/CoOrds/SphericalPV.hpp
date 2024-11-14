// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/SphericalPV.hpp":               //
//                      Spherical Positions and Velocities                   //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/TopoCentricCOSes.h"

namespace SpaceBallistics
{
  //-------------------------------------------------------------------------//
  // RATIONALE:                                                              //
  //-------------------------------------------------------------------------//
  // (*) For Rectangular CoOrds, we use separate types for the COS and for the
  //     dimensioned "Vector3D"s (Pos, Vel, Acc, Force etc) in that COS.
  // (*) However, for Spherical and Ellipsoidal co-ords,   we define types for
  //     the co-ords and their "dots" together, and those  types by themselves
  //     carry the info on the COS used.
  // (*) Such a difference between representation of Rectangular and Spherical
  //     (or Ellipsoidal) co-ords in because in the former case, we frequently
  //     need perform linear algebra operations with the corresp  "Vector3D"s;
  //     and in the latter case, the co-ords objs are typically used for stor-
  //     ing the PV data (Position and Velocity) as input or output only...
  //
  //=========================================================================//
  // "SphericalPV" Class:                                                    //
  //=========================================================================//
  // Position and Velicity in the Spherical Co-Ords corresponding to the given
  // rectangular COS.
  // In particular, for the GeoCentric  or Geo-TopoCentric Equatorial COS,  it
  // represents the classical astronomical Right Ascention, Declination, Radius-
  // Vector and their "Dots" of an Object:
  //
  template<typename COS>
  class SphericalPV
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Position:
    Angle  m_alpha;     // Right Ascention or Longitude
    Angle  m_delta;     // Declination     or Latitude
    LenK   m_rho;       // Distance (Radius-Vector)

    // Velocities:
    AngVel m_alphaDot;
    AngVel m_deltaDot;
    VelK   m_rhoDot;    // Radial Velocity

  public:
    // Default Ctor,  Copy Ctor, Assignment and Equality are auto-generated;
    // in particular, the Default Ctor initialises all flds to 0:
    SphericalPV             ()                         = default;
    SphericalPV             (SphericalPV const&)       = default;
    SphericalPV& operator=  (SphericalPV const&)       = default;
    bool         operator== (SphericalPV const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing "SphericalPV" from the Rectangular PV Vectors:
    //
    constexpr SphericalPV
    (
      PosKV<COS> const& a_pos, // Must be non-0
      VelKV<COS> const& a_vel  // May  be 0
    )
    : SphericalPV()            // Zero-out all components by default
    {
      // Position:
      m_rho      = LenK(a_pos);
      assert(IsPos(m_rho));

      double sinDelta = double(a_pos.z() / m_rho);
      double cosDelta = SqRt(1.0 - Sqr(sinDelta)); // cosDelta >= 0
      m_delta         = Angle(ASin (sinDelta));

      LenK   rhoXY    = m_rho * cosDelta;          // "rho" proj onto XY plane
      assert(Sqr(rhoXY).ApproxEquals(Sqr(a_pos.x()) + Sqr(a_pos.y())));
      double sinAlpha = double(a_pos.y() / rhoXY);
      double cosAlpha = double(a_pos.x() / rhoXY);
      m_alpha         = Angle(ATan2(a_pos.y(), a_pos.x()));
      auto   V2       = Sqr(a_vel.x()) + Sqr(a_vel.y()) + Sqr(a_vel.z());

      if (IsPos(V2))
      {
        // Radial Velocity: Projection of Velocity onto the Radius-Vector:
        m_rhoDot   = (a_pos.x() * a_vel.x()  +
                      a_pos.y() * a_vel.y()  +
                      a_pos.z() * a_vel.z()) / m_rho;

        // Decl Dot:
        m_deltaDot = 1.0_rad *
                     (cosDelta * a_vel.z() - sinDelta * sinAlpha * a_vel.y() -
                      sinDelta * cosAlpha  * a_vel.x()) / m_rho;

        // RA Dot: May be infinite for delta = +- Pi/2:
        m_alphaDot = 1.0_rad *
                     (cosAlpha * a_vel.y() - sinAlpha * a_vel.x()) / rhoXY;

        // Verification:
        DEBUG_ONLY
        (
          // Transversal Velocity Squared:
          auto   Vtr2 = V2 - Sqr(m_rhoDot);
          assert(!IsNeg(Vtr2));   // Up to rounding errors?

          assert((Sqr(m_deltaDot) + Sqr(cosDelta * m_alphaDot)).ApproxEquals
                 (Sqr(1.0_rad) * Vtr2 / Sqr(m_rho)));
        )
      }
    }

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    constexpr Angle  GetAlpha   () const { return m_alpha;    }
    constexpr Angle  GetDelta   () const { return m_delta;    }
    constexpr LenK   GetDist    () const { return m_rho;      }
    constexpr VelK   GetRadVel  () const { return m_rhoDot;   }
    constexpr AngVel GetAlphaDot() const { return m_alphaDot; }
    constexpr AngVel GetDeltaDot() const { return m_deltaDot; }

    //-----------------------------------------------------------------------//
    // Other Way Round: Pos and Vel Vectors from "SphericalPV":              //
    //-----------------------------------------------------------------------//
    // Returning both vectors together is more efficient:
    //
    constexpr std::pair<PosKV<COS>, VelKV<COS>> GetPVVectors() const
    {
      double cA = Cos(double(m_alpha));
      double sA = Sin(double(m_alpha));
      double cD = Cos(double(m_delta));
      double sD = Sin(double(m_delta));

      LenK   pos[3] { m_rho * cD * cA, m_rho * cD * sA, m_rho * sD };
      assert(Sqr(m_rho).ApproxEquals (Sqr(pos[0]) + Sqr(pos[1]) + Sqr(pos[2])));

      VelK   vel[3]
      {
        m_rhoDot * cD * cA - m_rho * sD * cA * m_deltaDot / 1.0_rad
                           - m_rho * cD * sA * m_alphaDot / 1.0_rad,
        m_rhoDot * cD * sA - m_rho * sD * sA * m_deltaDot / 1.0_rad
                           + m_rho * cD * cA * m_alphaDot / 1.0_rad,
        m_rhoDot * sD      + m_rho * cD      * m_deltaDot / 1.0_rad
      };
      DEBUG_ONLY
      (
        auto V2 = Sqr(m_rhoDot) + Sqr(m_rho * cD * m_alphaDot / 1.0_rad) +
                                  Sqr(m_rho *      m_deltaDot / 1.0_rad);
        assert(V2.ApproxEquals(Sqr(vel[0])  + Sqr(vel[1]) + Sqr(vel[2])));
      )
      return std::make_pair(PosKV<COS>(pos), VelKV<COS>(vel));
    }
  };

  //-------------------------------------------------------------------------//
  // Aliases (for Body-Centric Equatorial "SpherPV"s only):                  //
  //-------------------------------------------------------------------------//
  // XXX: Currently, they are provided for Body-Centric Equatorial COSes only:
  using HelioCentricEqSpherPV   = SphericalPV<HelioCentricEqFixCOS>;
  using HermeoCentricEqSpherPV  = SphericalPV<HermeoCentricEqFixCOS>;
  using CytheroCentricEqSpherPV = SphericalPV<CytheroCentricEqFixCOS>;
  using GeoCentricEqSpherPV     = SphericalPV<GeoCentricEqFixCOS>;
  using SelenoCentricEqSpherPV  = SphericalPV<SelenoCentricEqFixCOS>;
  using AreoCentricEqSpherPV    = SphericalPV<AreoCentricEqFixCOS>;
  using ZenoCentricEqSpherPV    = SphericalPV<ZenoCentricEqFixCOS>;
  using CronoCentricEqSpherPV   = SphericalPV<CronoCentricEqFixCOS>;
  using UranoCentricEqSpherPV   = SphericalPV<UranoCentricEqFixCOS>;
  using PoseidoCentricEqSpherPV = SphericalPV<PoseidoCentricEqFixCOS>;
  using HadeoCentricEqSpherPV   = SphericalPV<HadeoCentricEqFixCOS>;
}
// End namespace SpaceBallistics

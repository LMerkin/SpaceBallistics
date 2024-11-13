// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/SpherPV.hpp":                  //
//            Spherical and Ellipsoidal Positions and Velocities             //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/TopoCentricCOSes.h"

namespace SpaceBallistics
{
  // NB:
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
  // "SpherPV" Class:                                                        //
  //=========================================================================//
  // BodyCentric or Body-TopoCentric Spherical Co-Ords,    corresp to
  // "BodyCentric{Eq,Ecl}FixCOS" or "TopoCentric{Eq,Ecl}FixCOS", res.
  // and Body's Orbital Plane @J2000.0).
  // For the GeoCentric or Geo-TopoCentric Equatorial system, it represents the
  // classical astronomical Right Ascention, Declination, Radius-Vector and
  // their "Dots" of an Object.
  // This class provides TYPE SAFETY OF APPARENT CO-ORDS (RA and Decl), as it
  // specifies "from where they are apparent"!
  // It is normally used for Equatorial COSes, but Ecliptical ones can also be
  // used, in which case "alpha" is an Ecliptical Longitude and "delta" is an
  // Ecliptical Latitude:
  //
  template<typename COS>
  class SpherPV
  {
  private:
    //-----------------------------------------------------------------------//
    // Checks on the "COS":                                                  //
    //-----------------------------------------------------------------------//
    // The "COS" must have Fixed Axes (non-Rotating), though its Origin does
    // not need to be Fixed:
    static_assert(COS::HasFixedAxes);

    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Position:
    Angle  m_alpha;     // Right Ascention or Ecliptical Longitude
    Angle  m_delta;     // Declination     or Ecliptical Latitude
    LenK   m_rho;       // Distance (Radius-Vector)

    // Velocities:
    AngVel m_alphaDot;
    AngVel m_deltaDot;
    VelK   m_rhoDot;    // Radial Velocity

  public:
    // Default Ctor,  Copy Ctor, Assignment and Equality are auto-generated;
    // in particular, the Default Ctor initialises all flds to 0:
    SpherPV             ()                     = default;
    SpherPV             (SpherPV const&)       = default;
    SpherPV& operator=  (SpherPV const&)       = default;
    bool     operator== (SpherPV const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing "SpherPV" from the Rectangular PV Vectors:
    //
    constexpr SpherPV
    (
      PosKV<COS> const& a_pos, // Must be non-0
      VelKV<COS> const& a_vel  // May  be 0
    )
    : SpherPV()                // Zero-out all components by default
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
    constexpr Angle  GetRA     () const { return m_alpha;    } // Right Ascentn
    constexpr Angle  GetDecl   () const { return m_delta;    } // Declination
    constexpr LenK   GetDist   () const { return m_rho;      } // Distance (RV)
    constexpr VelK   GetRadVel () const { return m_rhoDot;   } // Radial Vel
    constexpr AngVel GetRADot  () const { return m_alphaDot; }
    constexpr AngVel GetDeclDot() const { return m_deltaDot; }
  };

  //-------------------------------------------------------------------------//
  // Aliases (for Body-Centric Equatorial "SpherPV"s only):                  //
  //-------------------------------------------------------------------------//
  // XXX: Currently, they are provided for Body-Centric Equatorial COSes only:
  using HelioCentricEqSpherPV   = SpherPV<HelioCentricEqFixCOS>;
  using HermeoCentricEqSpherPV  = SpherPV<HermeoCentricEqFixCOS>;
  using CytheroCentricEqSpherPV = SpherPV<CytheroCentricEqFixCOS>;
  using GeoCentricEqSpherPV     = SpherPV<GeoCentricEqFixCOS>;
  using SelenoCentricEqSpherPV  = SpherPV<SelenoCentricEqFixCOS>;
  using AreoCentricEqSpherPV    = SpherPV<AreoCentricEqFixCOS>;
  using ZenoCentricEqSpherPV    = SpherPV<ZenoCentricEqFixCOS>;
  using CronoCentricEqSpherPV   = SpherPV<CronoCentricEqFixCOS>;
  using UranoCentricEqSpherPV   = SpherPV<UranoCentricEqFixCOS>;
  using PoseidoCentricEqSpherPV = SpherPV<PoseidoCentricEqFixCOS>;
  using HadeoCentricEqSpherPV   = SpherPV<HadeoCentricEqFixCOS>;
}
// End namespace SpaceBallistics

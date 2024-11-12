// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/SpherEllipsPVs.hpp":                //
//            Spherical and Ellipsoidal Positions and Velocities             //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"

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
  // "BodyCentricSpherPV":                                                   //
  //=========================================================================//
  // Body-Centric Equatorial Spherical Co-Ords, corresp to "BodyCentricEqFixCOS"
  // (with Body's Equator and Body's Orbital Plane @J2000.0). For the GeoCentric
  // system, is used  to  represent the GeoCentric Right Ascention, Declination,
  // Radius-Vector and their "dots" for a given object:
  //
  template<Body BodyName>
  class BodyCentricSpherPV
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Position:
    Angle  m_alpha;     // Right Ascention
    Angle  m_delta;     // Declination
    LenK   m_rho;       // Body-Centric Radius-Vector

    // Velocities:
    AngVel m_alphaDot;
    AngVel m_deltaDot;
    VelK   m_rhoDot;    // Radial Velocity

  public:
    // Default Ctor,  Copy Ctor, Assignment and Equality are auto-generated
    // In particular, the Default Ctor initialises all flds to 0:
    BodyCentricSpherPV            ()                                = default;
    BodyCentricSpherPV            (BodyCentricSpherPV const&)       = default;
    BodyCentricSpherPV& operator= (BodyCentricSpherPV const&)       = default;
    bool                operator==(BodyCentricSpherPV const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing "BodyCentricSpherPV" from the Rectangular PV Vectors:
    //
    constexpr BodyCentricSpherPV
    (
      PosKV<BodyCentricEqFixCOS<BodyName>> const& a_pos, // Must be non-0
      VelKV<BodyCentricEqFixCOS<BodyName>> const& a_vel  // May  be 0
    )
    : BodyCentricSpherPV()            // Zero-out all components by default
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
  // Aliases:                                                                //
  //-------------------------------------------------------------------------//
  using HelioCentricSpherPV   = BodyCentricSpherPV<Body::Sun>;
  using HermeoCentricSpherPV  = BodyCentricSpherPV<Body::Mercury>;
  using CytheroCentricSpherPV = BodyCentricSpherPV<Body::Venus>;
  using GeoCentricSpherPV     = BodyCentricSpherPV<Body::Earth>;
  using SelenoCentricSpherPV  = BodyCentricSpherPV<Body::Moon>;
  using AreoCentricSpherPV    = BodyCentricSpherPV<Body::Mars>;
  using ZenoCentricSpherPV    = BodyCentricSpherPV<Body::Jupiter>;
  using CronoCentricSpherPV   = BodyCentricSpherPV<Body::Saturn>;
  using UranoCentricSpherPV   = BodyCentricSpherPV<Body::Uranus>;
  using PoseidoCentricSpherPV = BodyCentricSpherPV<Body::Neptune>;
  using HadeoCentricSpherPV   = BodyCentricSpherPV<Body::PlChB>;
}
// End namespace SpaceBallistics

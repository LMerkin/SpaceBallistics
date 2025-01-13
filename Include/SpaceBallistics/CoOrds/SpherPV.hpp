// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/SpherPV.hpp":                  //
//                     Spherical Positions and Velocities                    //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // RATIONALE:                                                              //
  //=========================================================================//
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
  // Position and Velicity in the Spherical Co-Ords corresponding to the given
  // rectangular COS.
  // In particular, for the GeoCentric  or Geo-TopoCentric Equatorial COS,  it
  // represents the classical astronomical Right Ascention, Declination, Radius-
  // Vector and their "Dots" of an Object:
  //
  template<typename COS, Body B = Body::UNDEFINED>
  class SpherPV
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

    // TimeStamp of the COS from "Vector3D"s used in construction of this
    // "SpherPV" (if any):
    typename Bits::TSWrapper<COS>::TS m_cosTS;

  public:
    //-----------------------------------------------------------------------//
    // Default Ctor:                                                         //
    //-----------------------------------------------------------------------//
    // Initialises all flds to 0 except the TimeStamp:
    SpherPV     ()
    : m_alpha   (),
      m_delta   (),
      m_rho     (),
      m_alphaDot(),
      m_deltaDot(),
      m_rhoDot  (),
      m_cosTS   (Bits::TSWrapper<COS>::TS::UnDef())
    {}

    // Copy Ctor and Assignment are auto-generated:
    SpherPV             (SpherPV const&)       = default;
    SpherPV& operator=  (SpherPV const&)       = default;

    // Equality must take into account the COS TS:
    bool     operator== (SpherPV const& a_right) const
    {
      assert(m_cosTS.IsUndef() || a_right.m_cosTS.IsUnDef() ||
             m_cosTS == a_right.m_cosTS);
      return
        m_alpha    == a_right.m_alpha    && m_delta    == a_right.m_delta    &&
        m_rho      == a_right.m_rho      && m_alphaDot == a_right.m_alphaDot &&
        m_deltaDot == a_right.m_deltaDot && m_rhoDot   == a_right.m_rhoDot;
    }

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing "SpherPV" from the Rectangular PV Vectors:
    //
    constexpr explicit SpherPV
    (
      PosKV<COS, B> const& a_pos,                   // Must be non-0
      VelKV<COS, B> const& a_vel = VelKV<COS, B>()  // May  be     0
    )
    : SpherPV()     // Zero-out all components by default
    {
      // Position:
      m_rho           = LenK(a_pos);
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

          // Also, COS TS of "a_vel" must be compatible with that of "a_pos":
          a_pos.CheckCOSTSs(a_vel.GetCOSTS());
        )
      }
      // Then it is sufficient to propagate the COST TS from "a_pos":
      m_cosTS = a_pos.GetCOSTS();
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
    // Other Way Round: Pos and Vel Vectors from "SpherPV":                  //
    //-----------------------------------------------------------------------//
    // Returning both vectors together for efficiency:
    //
    constexpr std::pair<PosKV<COS, B>, VelKV<COS, B>> GetPVVectors() const
    {
      double cA = Cos(m_alpha);
      double sA = Sin(m_alpha);
      double cD = Cos(m_delta);
      double sD = Sin(m_delta);

      LenK   x  = m_rho * cD * cA;
      LenK   y  = m_rho * cD * sA;
      LenK   z  = m_rho * sD;
      assert(Sqr(m_rho).ApproxEquals(Sqr(x) + Sqr(y) + Sqr(z)));

      VelK   Vx = m_rhoDot   * cD * cA
                - m_rho * sD * cA * m_deltaDot / 1.0_rad
                - m_rho * cD * sA * m_alphaDot / 1.0_rad;

      VelK   Vy = m_rhoDot   * cD * sA
                - m_rho * sD * sA * m_deltaDot / 1.0_rad
                + m_rho * cD * cA * m_alphaDot / 1.0_rad;

      VelK   Vz = m_rhoDot * sD
                + m_rho * cD      * m_deltaDot / 1.0_rad;

      DEBUG_ONLY
      (
        auto V2 = Sqr(m_rhoDot) + Sqr(m_rho * cD * m_alphaDot / 1.0_rad) +
                                  Sqr(m_rho *      m_deltaDot / 1.0_rad);
        assert(V2.ApproxEquals(Sqr(Vx) + Sqr(Vy) + Sqr(Vz)));
      )
      return std::make_pair
             (
               PosKV<COS, B>(m_cosTS,  x,  y,  z),
               VelKV<COS, B>(m_cosTS, Vx, Vy, Vz)
             );
    }
  };

  //-------------------------------------------------------------------------//
  // Aliases (for Body-Centric Equatorial "SpherPV"s only):                  //
  //-------------------------------------------------------------------------//
  // XXX: Currently, they are provided for Body-Centric Equatorial COSes  only.
  // In particular, "GeoCEqSpherPV" provides GeoCentric Right Ascentions, Decl-
  // inations and Distances:
  //
  template<Body B = Body::UNDEFINED>
  using HelioCEqSpherPV   = SpherPV<HelioCEqFixCOS,   B>;

  template<Body B = Body::UNDEFINED>
  using HermeoCEqSpherPV  = SpherPV<HermeoCEqFixCOS,  B>;

  template<Body B = Body::UNDEFINED>
  using CytheroCEqSpherPV = SpherPV<CytheroCEqFixCOS, B>;

  template<Body B = Body::UNDEFINED>
  using GeoCEqSpherPV     = SpherPV<GeoCEqFixCOS,     B>;

  template<Body B = Body::UNDEFINED>
  using SelenoCEqSpherPV  = SpherPV<SelenoCEqFixCOS,  B>;

  template<Body B = Body::UNDEFINED>
  using AreoCEqSpherPV    = SpherPV<AreoCEqFixCOS,    B>;

  template<Body B = Body::UNDEFINED>
  using ZenoCEqSpherPV    = SpherPV<ZenoCEqFixCOS,    B>;

  template<Body B = Body::UNDEFINED>
  using CronoCEqSpherPV   = SpherPV<CronoCEqFixCOS,   B>;

  template<Body B = Body::UNDEFINED>
  using UranoCEqSpherPV   = SpherPV<UranoCEqFixCOS,   B>;

  template<Body B = Body::UNDEFINED>
  using PoseidoCEqSpherPV = SpherPV<PoseidoCEqFixCOS, B>;

  template<Body B = Body::UNDEFINED>
  using HadeoCEqSpherPV   = SpherPV<HadeoCEqFixCOS,   B>;
}
// End namespace SpaceBallistics

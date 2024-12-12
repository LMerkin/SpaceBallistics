// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/EllipsPV.hpp":                 //
//                     Ellipsoidal Positions and Velocities                  //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/PhysForces/BodyData.hpp"

namespace SpaceBallistics
{
  //-------------------------------------------------------------------------//
  // RATIONALE:                                                              //
  //-------------------------------------------------------------------------//
  // Mostly useful for computation of Locations on the Ellipsoidal Bodies and
  // SpaceCraft GroundTracks on such Bodies.
  // XXX: Current Restriction: "EllipsPV" makes sense only for BodyC Rotating
  // COSes (if the corrsp Body is a Rotational Ellipsoid or Spheroid):
  //
  //=========================================================================//
  // "EllipsPV" Class:                                                       //
  //=========================================================================//
  template<Body BBody, Body B = Body::UNDEFINED>
  class EllipsPV
  {
  private:
    //-----------------------------------------------------------------------//
    // Consts:                                                               //
    //-----------------------------------------------------------------------//
    // NB: Unlike "SpherPV", "EllipsPV" requires the dimensions of the Ellipso-
    // idal Surface:
    //
    // Equatorial Radius:
    constexpr static LenK   Re    = BodyData<BBody>::Re;
    constexpr static auto   Re2   = Sqr(Re);

    // Polar Radius:
    constexpr static LenK   Rp    = BodyData<BBody>::Rp;
    constexpr static auto   Rp2   = Sqr(Rp);
    static_assert(Rp <= Re);

    // Flattening:
    constexpr static double FlatC = double(Rp / Re);
    constexpr static double Flat  = 1.0 - FlatC;

    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Position: Body-detic (Body-graphic) Co-Ords (NOT BodyCentric!):
    Angle  m_lambda;          // Longitide
    Angle  m_phi;             // Latitude  (normal to the Body Surface)
    LenK   m_h;               // Elevation (normal to the Body Surface)

    // Velocities:
    AngVel m_lambdaDot;
    AngVel m_phiDot;
    VelK   m_hDot;

  public:
    // Default Ctor,  Copy Ctor, Assignment and Equality are auto-generated;
    // in particular, the Default Ctor initialises all flds to 0:
    EllipsPV             ()                      = default;
    EllipsPV             (EllipsPV const&)       = default;
    EllipsPV& operator=  (EllipsPV const&)       = default;
    bool      operator== (EllipsPV const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing "SpherPV" from the Rectangular PV Vectors:
    //
    constexpr explicit EllipsPV
    (
      PosKVRot<BBody, B> const& a_pos,                        // Must be non-0
      VelKVRot<BBody, B> const& a_vel = VelKVRot<BBody, B>()  // May  be     0
    )
    : EllipsPV()         // Zero-out all components by default
    {
      // FIXME TODO
    }

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    constexpr Angle  GetLambda   () const { return m_lambda;    }
    constexpr Angle  GetPhi      () const { return m_phi;       }
    constexpr LenK   GetH        () const { return m_h;         }
    constexpr VelK   GetVertVel  () const { return m_hDot;      }
    constexpr AngVel GetLambdaDot() const { return m_lambdaDot; }
    constexpr AngVel GetPhiDot   () const { return m_phiDot;    }

    //-----------------------------------------------------------------------//
    // Other Way Round: Pos and Vel Vectors from "EllipsPV":                 //
    //-----------------------------------------------------------------------//
    // Returning both vectors together is more efficient:
    //
    constexpr std::pair<PosKVRot<BBody, B>, VelKVRot<BBody>, B>
    GetPVVectors() const
    {
      double   cosL  = Cos(double(m_lambda));
      double   sinL  = Sin(double(m_lambda));
      double   cosP  = Cos(double(m_phi));
      double   sinP  = Sin(double(m_phi));

      // Position:
      double   tanP  = sinP / cosP;
      double   tanP2 = Sqr(tanP);
      auto     RT2   = Rp2  * tanP2;
      auto     d2    = Re2  + RT2;
      LenK     d     = SqRt(d2);
      // (x,z) in the cross-section through the Minor (Rotational) Axis and the
      // given point:
      LenK     X     = Re2 / d  + m_h * cosP;
      LenK     Z     = RT2 / d  + m_h * sinP;
      PosKVRot<BBody, B> pos(X * cosL, X  * sinL, Z);

      // Velocity:
      VelKVRot<BBody, B> vel;   // Initially, all 0s
      if (!(IsZero(m_lambdaDot) && IsZero(m_phiDot) && IsZero(m_hDot)))
      {
        double t2d   =   tanP * (1.0 + tanP2);
        VelK   ddot  =   Rp2  * t2d  / d * m_phiDot     / 1.0_rad;
        VelK   Xdot  = - Re2  / d2   * ddot
                       + m_hDot      * cosP
                       - m_h         * sinP * m_phiDot  / 1.0_rad;

        VelK   Zdot  =   Rp2 * (2.0  * t2d / d  * m_phiDot / 1.0_rad -
                                tanP2      / d2 * ddot)
                       + m_hDot      * sinP
                       + m_h         * cosP * m_phiDot  / 1.0_rad;

        vel.x() = Xdot * cosL  - X * sinL * m_lambdaDot / 1.0_rad;
        vel.y() = Xdot * sinL  + X * cosL * m_lambdaDot / 1.0_rad;
        vel.z() = Zdot;
      }
      return std::make_pair(pos, vel);
    }
  };

  //-------------------------------------------------------------------------//
  // Aliases:                                                                //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using HelioCEllipsPV   = EllipsPV<Body::Sun,     B>;

  template<Body B = Body::UNDEFINED>
  using HermeoCEllipsPV  = EllipsPV<Body::Mercury, B>;

  template<Body B = Body::UNDEFINED>
  using CytheroCEllipsPV = EllipsPV<Body::Venus,   B>;

  template<Body B = Body::UNDEFINED>
  using GeoCEllipsPV     = EllipsPV<Body::Earth,   B>;

  template<Body B = Body::UNDEFINED>
  using SelenoCEllipsPV  = EllipsPV<Body::Moon>,   B;

  template<Body B = Body::UNDEFINED>
  using AreoCEllipsPV    = EllipsPV<Body::Mars,    B>;

  template<Body B = Body::UNDEFINED>
  using ZenoCEllipsrPV   = EllipsPV<Body::Jupiter, B>;

  template<Body B = Body::UNDEFINED>
  using CronoCEllipsPV   = EllipsPV<Body::Saturn,  B>;

  template<Body B = Body::UNDEFINED>
  using UranoCEllipsPV   = EllipsPV<Body::Uranus,  B>;

  template<Body B = Body::UNDEFINED>
  using PoseidoCEllipsPV = EllipsPV<Body::Neptune, B>;

  // XXX: No Ellipsoid data for Pluto, hence no "EllipsPV"...
}
// End namespace SpaceBallistics

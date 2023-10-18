// vim:ts=2:et
//===========================================================================//
//                            "ConstrElement.cpp":                           //
//                Geometrical Objects as Construction Elements               //
//===========================================================================//
#include "SpaceBallistics/ConstrElement.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // Default Ctor:                                                           //
  //=========================================================================//
  ConstrElement::ConstrElement()
  : m_mass        (0.0),
    m_surfArea    (0.0),
    m_vol         (0.0),
    m_CoM         {0.0_m, 0.0_m, 0.0_m},
    m_MoIY        (0.0),
    m_massIsFinal (false)
  {}

  //=========================================================================//
  // Copy Ctor:                                                              //
  //=========================================================================//
  ConstrElement::ConstrElement(ConstrElement const& a_right)
  : m_mass        (a_right.m_mass),
    m_surfArea    (a_right.m_surfArea),
    m_vol         (a_right.m_vol),
    m_CoM         {a_right.m_CoM[0], a_right.m_CoM[1], a_right.m_CoM[2]},
    m_MoIY        (a_right.m_MoIY),
    m_massIsFinal (a_right.m_massIsFinal)
  {}

  //=========================================================================//
  // Non-Default Ctor:                                                       //
  //=========================================================================//
  // For internal use only. NB: The "a_mass" param may be final or provisional,
  // and the "a_mass_is_valid" is set accordingly:
  ConstrElement::ConstrElement
  (
    Mass a_mass,  Area a_surf_area,  Vol a_vol,
    Len  a_xc,    Len  a_yc,         Len a_zc,
    MoI  a_moi_y, bool a_mass_is_final
  )
  : m_mass        (a_mass),
    m_surfArea    (a_surf_area),
    m_vol         (a_vol),
    m_CoM         {a_xc, a_yc, a_zc},
    m_MoIY        (a_moi_y),
    m_massIsFinal (a_mass_is_final)
  {
    assert(!(IsNeg(m_mass) || IsNeg(m_surfArea) || IsNeg(m_vol) ||
             IsNeg(m_MoIY)));
  }

  //=========================================================================//
  // "SetMass":                                                              //
  //=========================================================================//
  // XXX: A "setter" for the Mass (invoked when we know the final mass, set it
  // and pro-rate the mass-dependent params:
  //
  void ConstrElement::SetMass(Mass a_new_mass)
  {
    if (UNLIKELY(m_massIsFinal))
      // Cannot adjust the Mass once it has been finalised:
      throw std::logic_error("ConstrElement::SetMass: Already set");

    if (UNLIKELY(!(IsPos(m_mass) && IsPos(a_new_mass))))
      throw std::logic_error
            ("ConstrElement::SetMass: Curr or New Mass is <= 0");

    // If OK:
    double scaleCoeff = double(a_new_mass / m_mass);
    assert(scaleCoeff > 0.0);

    // Set the Mass and adjust the MoI:
    m_mass        = a_new_mass;
    m_MoIY       *= scaleCoeff;
    m_massIsFinal = true;
  }

  //=========================================================================//
  // "ShellTrCone" Generator:                                                //
  //=========================================================================//
  // See the detailed comments in the Header file:
  //
  ConstrElement ConstrElement::ShellTrCone
  (
    Len     a_x0,
    Len     a_y0,
    Len     a_z0,
    double  a_alpha, // Dimension-less, in (-Pi/2 .. Pi/2)
    Len     a_d0,    // Base diameter at x0
    Len     a_d1,    // Base diameter at the other section (X = x1 > x0)
    Len     a_h,     // Must be > 0
    Mass    a_mass
  )
  {
    if (UNLIKELY
       (!(IsPos(a_h)     && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
           !(IsZero(a_d0) && IsZero(a_d1)) && std::fabs(a_alpha) < M_PI_2)))
      throw std::invalid_argument
            ("ConstrElement::ShellTrCone: Param(s) out of range");

    // At least one of "z0", "y0" must be 0:
    bool inXY    = IsZero(a_z0);
    bool inXZ    = IsZero(a_y0);
    if (UNLIKELY(!(inXY || inXZ)))
      throw std::invalid_argument
            ("ConstrElement::ShellTrCone: Y0 or Z0 must be 0 (or both)");
    // The mass may or may not be given:
    bool hasMass = IsPos(a_mass);

    double cosA  = std::cos(a_alpha);
    double sinA  = std::sin(a_alpha);
    double tanA  = sinA / cosA;

    auto   r     = a_d0 / 2.0;
    auto   R     = a_d1 / 2.0;
    auto   R2    = Sqr(R);
    auto   r2    = Sqr(r);
    auto   s     = R + r;
    auto   th2   = 2.0 * Sqr(a_h);
    // Moment of Inertia wrt OY per Surface Density. Must be strictly > 0.
    // In the XY case, the result is invariant of the sign of "alpha";
    // in the XZ case, it is invariant iunder the transform
    //                 z0 <-> (-z0), alpha <-> (-alpha)   :
    auto   L4    =
      (M_PI/4.0) * SqRt(Sqr(R - r) + Sqr(a_h)) *
      (inXY
       ? (cosA * ((16.0/3.0) * (R + r/2.0) * a_h * a_x0 -
                  cosA * (R2 * s  + (r2 - th2) * R + (r2 - th2/3.0) * r)) +
          2.0 * s * (R2 + r2 + 2.0 * Sqr(a_x0)))
       : (a_h *  ((16.0/3.0) * (R + r/2.0) * (a_x0 * cosA + a_z0 * sinA)  +
                  2.0 * a_h  * (R + r/3.0))                               +
          s * (R2 + r2 + 4.0 * (Sqr(a_x0) + Sqr(a_z0))))
      );
    assert(IsPos(L4));

    // Side Surface Area and Notional Volume:
    auto S       = SideSurfArea_TrCone(a_d0, a_d1, a_h);
    auto V       = Volume_TrCone      (a_d0, a_d1, a_h);

    // If the mass is given, calculate the Surface Density; otherwise, assume
    // it is 1.0. Then set the Mass and MoI:
    auto surfDens = hasMass ? (a_mass / S) : SurfDens(1.0);
    auto M        = hasMass ?  a_mass      : (S * surfDens);
    auto MoIY     = L4 * surfDens;

    // Co-ords of the CoM:
    // The X-coord is such that the above L4 is minimal. Other co-ords are
    // computed from the rotation axis geometry:
    auto deltaXC = cosA * (2.0 * R + r) / (3.0 * (R + r)) * a_h;
    auto xC      = a_x0 + deltaXC;
    auto yC      = inXY ? (a_y0 + deltaXC * tanA) : 0.0_m;
    auto zC      = inXZ ? (a_z0 + deltaXC * tanA) : 0.0_m;

    return ConstrElement(M, S, V, xC, yC, zC, MoIY, hasMass);
  }

  //=========================================================================//
  // "ShellSpherSegm" Generator:                                             //
  //=========================================================================//
  // See the detailed comments in the Header:
  //
  ConstrElement ConstrElement::ShellSpherSegm
  (
    bool    a_is_pos_facing,
    Len     a_x0,           // X-coord of the base center (NOT of the pole!)
    Len     a_y0,
    Len     a_z0,
    double  a_alpha,        // Dimension-less, in (-Pi/2 .. Pi/2)
    Len     a_d,            // Base diameter
    Len     a_h,            // Height: 0 < a_h <= a_d/2.0
    Mass    a_mass
  )
  {
    // NB: We must always have 0 < h <= R:
    if (UNLIKELY(!(IsPos(a_d) && IsPos(a_h) && a_h <= a_d/2.0 &&
                 std::fabs(a_alpha) < M_PI_2)))
      throw std::invalid_argument
            ("ConstrElement::ShellSpeherSegm: Param(s) out of range");
    // The mass may or may not be given:
    bool hasMass = IsPos(a_mass);

    // At least one of "z0", "y0" must be 0:
    bool inXY    = IsZero(a_z0);
    bool inXZ    = IsZero(a_y0);
    if (UNLIKELY(!(inXY || inXZ)))
      throw std::invalid_argument
            ("ConstrElement::ShellSpherSegm: Y0 or Z0 must be 0 (or both)");

    double cosA  = std::cos(a_alpha);
    double sinA  = std::sin(a_alpha);
    double tanA  = sinA / cosA;

    // X- and Z-coords of the Pole; the Y-coord is not used:
    double sgn   = a_is_pos_facing ? 1.0 : -1.0;
    auto   xP    = a_x0 +  sgn  * cosA * a_h;
    auto   zP    = inXZ ? (a_z0 + sgn  * sinA * a_h) : 0.0_m;
    auto   r     = a_d / 2.0;                    // Base   radius
    auto   R     = (Sqr(r) / a_h + a_h) / 2.0;   // Sphere radius
    // Therefore, R >= r >= h, but we cannot formally assert this due to possi-
    // ble rounding effects...

    // NB: Formally, "L4" is invariant under the transform
    // isPosFacing -> !isPosFacing,   alpha -> alpha + Pi,
    // as one can expect; however, in reality, we always have |alpha| < Pi/2:
    auto   L4    =
      M_PI * R * a_h *
      (inXY
       ? (2.0 * (Sqr(xP)   + a_h * (R - sgn  * xP   * cosA)) -
          a_h * ((2.0/3.0) * a_h + (R - a_h) * cosA * cosA))
       : (2.0 * (Sqr(xP)   + Sqr(zP))                        +
          a_h * (a_h/3.0   - sgn * 2.0 * (xP * cosA + zP * sinA) + R))
      );
    assert(IsPos(L4));

    // Side Surface Area and Notional Volume:
    auto S       = SideSurfArea_SpherSegm(a_d, a_h);
    auto V       = Volume_SpherSegm      (a_d, a_h);

    // If the mass is given, calculate the Surface Density; otherwise, assume
    // it is 1.0. Then set the Mass and MoI:
    auto surfDens = hasMass ? (a_mass / S) : SurfDens(1.0);
    auto M        = hasMass ?  a_mass      : (S * surfDens);
    auto MoIY     = L4 * surfDens;

    // Co-ords of the CoM:
    // The X-coord is such that the above L4 is minimal. Other co-ords are
    // computed from the rotation axis geometry:
    // PosFacing (sgn = 1): xC < xP;
    // NegFacing (sgn =-1): xC > xP:
    auto deltaXC = -sgn * a_h * cosA / 2.0;
    auto xC      = a_x0 + deltaXC;
    auto yC      = inXY ? (xP + deltaXC * tanA) : 0.0_m;
    auto zC      = inXZ ? (zP + deltaXC * tanA) : 0.0_m;
    assert(( a_is_pos_facing && xC < xP) ||
           (!a_is_pos_facing && xC > xP));

    return ConstrElement(M, S, V, xC, yC, zC, MoIY, hasMass);
  }

  //=========================================================================//
  // "PointMass":                                                            //
  //=========================================================================//
  ConstrElement ConstrElement::PointMass
    (Len a_x0, Len a_y0, Len a_z0, Mass a_mass)
  {
    if (UNLIKELY(!IsPos(a_mass)))
      throw std::invalid_argument
            ("ConstrElement::PosintMass: Mass must be > 0");

    auto MoIY =  a_mass * (Sqr(a_x0) + Sqr(a_y0)); // Distance^2 to OY axis
    // The mass is considered to be final:
    return ConstrElement
           (a_mass, Area(0.0), Vol(0.0), a_x0, a_y0, a_z0, MoIY, true);
  }

  //=========================================================================//
  // Addition:                                                               //
  //=========================================================================//
  // The summands must both have final or non-final masses. Only one summand is
  // allowed to have 0 mass:
  //
  ConstrElement& ConstrElement::operator+= (ConstrElement const& a_right)
  {
    assert(!(IsNeg(m_mass)         || IsNeg(m_surfArea)         ||
             IsNeg(m_MoIY)         ||
             IsNeg(a_right.m_mass) || IsNeg(a_right.m_surfArea) ||
             IsNeg(a_right.m_MoIY)));

    if (UNLIKELY(m_massIsFinal != a_right.m_massIsFinal))
      throw std::invalid_argument
            ("ConstrElement::Add: MassIsFinal flags do not match");

    if (UNLIKELY(IsZero(m_mass) && IsZero(a_right.m_mass)))
      throw std::invalid_argument
            ("ConstrElement::Add: Masses must not be both 0");

    // Masses, SurfAreas, Vols and MoIs are directly-additive:
    auto m0     = m_mass;
    m_mass     += a_right.m_mass;
    m_surfArea += a_right.m_surfArea;
    m_vol      += a_right.m_vol;
    m_MoIY     += a_right.m_MoIY;
    // For the CoM, do the weighted average (but the total mass must be non-0):
    assert(IsPos(m_mass));
    double mu0 = double(m0             / m_mass);
    double mu1 = double(a_right.m_mass / m_mass);
    m_CoM[0]   = mu0 * m_CoM[0]  + mu1 * a_right.m_CoM[0];
    m_CoM[1]   = mu0 * m_CoM[1]  + mu1 * a_right.m_CoM[1];
    m_CoM[2]   = mu0 * m_CoM[2]  + mu1 * a_right.m_CoM[2];
    return *this;
  }
}
// End namespace SpaceBallistics

// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/ME/MechElement.hpp":                 //
//                   Geometric Objects as "Mechanical Elements"              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/MechElement.h"
#include <cmath>
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // "MechElement::Init":                                                   //
  //=========================================================================//
  // For use by Derived Classes. It is similar to a Non-Default Ctor, but is
  // intended to be invoked AFTER the Derived flds have been initialised, so
  // implemented as a Static Method, not a Ctor.
  // XXX: This method is made "static"; otherwise, it could not be called from
  // methods of derived classes on OTHER objs, even though it is "protected":
  // XXX: The params are just arrays, not ECOS-parameterised vectors/tensors:
  // since this function is internal, that is OK:
  //
  template<LVSC LVSCKind>
  constexpr void MechElement<LVSCKind>::Init
  (
    MechElement*    a_me,
    TT              a_ecos_ts,
    Len  const      a_com     [3],
    Vel  const      a_com_dots[3],
    Mass            a_mass,
    MassRate        a_mass_dot,
    Vol             a_encl_vol,
    MoI     const   a_mois    [3],
    MoIRate const   a_moi_dots[3],
    bool            a_is_final
  )
  {
    // NB: Even if "a_mass" is not final, it must be positive. In the context
    // under consideration, MassDot and MoIDots cannot be positive, since the
    // MechElement's mass is either constant or being exhausted:
    //
    assert(a_me != nullptr      &&
           IsPos(a_mass)        && !IsNeg(a_mois    [0]) &&
          !IsNeg(a_mois    [1]) && !IsNeg(a_mois    [2]) &&
          !IsPos(a_mass_dot)    && !IsPos(a_moi_dots[0]) &&
          !IsPos(a_moi_dots[1]) && !IsPos(a_moi_dots[1]) &&
          !IsNeg(a_encl_vol));

    a_me->m_CoM      .Init(a_ecos_ts, a_com);
    a_me->m_CoMDots  .Init(a_ecos_ts, a_com_dots);
    a_me->m_mass    = a_mass;
    a_me->m_massDot = a_mass_dot;
    a_me->m_enclVol = a_encl_vol;
    a_me->m_MoIs     .Init(a_ecos_ts, a_mois);
    a_me->m_MoIDots  .Init(a_ecos_ts, a_moi_dots);
    a_me->m_isFinal = a_is_final;

    // NB: HERE all "Vector3D"s are constructed, with the same "a_ecos_ts"!
  }

  //=========================================================================//
  // "MechElement" Default Ctor:                                             //
  //=========================================================================//
  // NB:
  // (*) All numeric fields except the TimeStamp are initialised to 0s, not to
  //     NaNs, so the "empty" obj (constructed by the Default Ctor) can be used
  //     as an initial value for summation ("+", "+=", etc);
  // (*) "IsFinal" is set to "true", because the empty "MechElement" is typic-
  //     ally used as the base for "+" which requires Final operands:
  // (*) Again, the same "TT::UnDef()" is installed in all "Vector3D"s:
  //
  template<LVSC LVSCKind>
  constexpr MechElement<LVSCKind>::MechElement()
  : m_CoM    (TT::UnDef(),         0.0_m,         0.0_m,        0.0_m),
    m_CoMDots(TT::UnDef(),     Vel(0.0),      Vel(0.0),     Vel(0.0)),
    m_mass   (0.0),
    m_massDot(0.0),
    m_enclVol(0.0),
    m_MoIs   (TT::UnDef(), MoI    (0.0), MoI     (0.0), MoI    (0.0)),
    m_MoIDots(TT::UnDef(), MoIRate(0.0), MoIRate (0.0), MoIRate(0.0)),
    m_isFinal(true)
  {}

  //=========================================================================//
  // "MechElement" Non-Default Ctor:                                         //
  //=========================================================================//
  template<LVSC LVSCKind>
  constexpr MechElement<LVSCKind>::MechElement
  (
    TT             a_ecos_ts,
    Len      const a_com     [3],
    Vel      const a_com_dots[3],
    Mass           a_mass,
    MassRate       a_mass_dot,
    Vol            a_encl_vol,
    MoI      const a_mois    [3],
    MoIRate  const a_moi_dots[3],
    bool           a_is_final
  )
  {
    Init(this,       a_ecos_ts, a_com,      a_com_dots,  a_mass,  a_mass_dot,
         a_encl_vol, a_mois,    a_moi_dots, a_is_final);
  }

  //=========================================================================//
  // "MechElement::GetMassScale":                                            //
  //=========================================================================//
  // Given a list of "MechElement"s  (all of them must have Non-Final Masses)
  // and the real Total Mass of them,  returns  a dimension-less factor which
  // could be applied (in "ProRateMass") to each "MechElement" to set its cor-
  // rect Final Mass. It assumes that the densities  (surface  or vol) of all
  // "MechElement"s in the list are the same, ie their relative masses remain
  // unchanged after scaling:
  //
  template<LVSC LVSCKind>
  constexpr double MechElement<LVSCKind>::GetMassScale
  (
    std::initializer_list<MechElement const*> a_mes,
    Mass                                      a_total_mass
  )
  {
    assert(IsPos(a_total_mass));

    // All "a_mes" must NOT be Final yet. Calculate their nominal total mass:
    Mass nomTotal = 0.0_kg;
    for (MechElement const* me: a_mes)
    {
      assert(me != nullptr && !(me->m_isFinal));
      nomTotal  += me->m_mass;
    }
    // Determine the Scale Factor (exists iff nomTotal > 0):
    assert(IsPos(nomTotal));
    double scale = double(a_total_mass / nomTotal);

    assert(scale > 0.0);
    return scale;
  }

  //=========================================================================//
  // "MechElement::ProRateMass":                                             //
  //=========================================================================//
  // Returns a new object (of any type derived from "MechElement") which diff-
  // ers from the original one  by the Mass, MoIs, MassDot and MoIDots being
  // multiplied by the ScaleFactor.
  // There are 2 distinct use cases of this function:
  // (1) "a_der" is a Non-Final ME; "a_der" was obtained by "GetMassScale"
  //     from a list of such Non-Final MEs (typically 2D Surfaces)  with a
  //     known total mass and a 
  // (2) "a_der" is a formally Final ME which  is a 3D Propellant Bulk; its
  //     mass is scaled because the Propellant is being replaced by Pressu-
  //     risation Gas; in this case,
  //
  template<LVSC LVSCKind>
  template<bool IsPressnGas, typename Derived>
  constexpr Derived MechElement<LVSCKind>::ProRateMass
  (
    Derived const& a_der,
    double         a_scale,
    DensRate       a_dens_dot
  )
  {
    static_assert(std::is_base_of_v<MechElement, Derived>);
    assert(a_scale > 0.0);
    assert(IsPressnGas || IsZero(a_dens_dot));

    // Create a copy of "a_der" using the Copy Ctor of Derived ( which must
    // indeed be derived from "MechElement"). In particulat, the ECOS TS is
    // copied:
    Derived copy(a_der);

    // Cannot adjust the Mass once it has been finalised, or if any Rates are
    // non-0, unless this is the "PressnGas" case:
    assert(IsPressnGas  ||
          (!copy.m_isFinal           && IsZero(copy.m_massDot)    &&
           IsZero(copy.m_MoIDots[0]) && IsZero(copy.m_MoIDots[1]) &&
           IsZero(copy.m_MoIDots[2])));

    // Set the Mass and adjust the MoIs; but the CoM is unchanged:
    copy.m_mass    *= a_scale;
    copy.m_MoIs[0] *= a_scale;
    copy.m_MoIs[1] *= a_scale;
    copy.m_MoIs[2] *= a_scale;
    copy.m_isFinal  = true;

    if constexpr(IsPressnGas)
    {
      // Adjust the Mass and MoI Dots using "a_dens_dot".
      // We use the fact that the Mass and MoIs are linear in Density, and
      // the ME must have a non-0 volume:
      assert(IsPos(copy.m_mass) && IsPos(copy.m_enclVol));
      Density    dens  = copy.m_mass  /  copy.m_enclVol;

      copy.m_massDot    *= a_scale;
      copy.m_massDot    += copy.m_enclVol * a_dens_dot;

      copy.m_MoIDots[0] *= a_scale;
      copy.m_MoIDots[0] += copy.m_MoIs[0] * a_dens_dot / dens;

      copy.m_MoIDots[1] *= a_scale;
      copy.m_MoIDots[1] += copy.m_MoIs[1] * a_dens_dot / dens;

      copy.m_MoIDots[2] *= a_scale;
      copy.m_MoIDots[2] += copy.m_MoIs[2] * a_dens_dot / dens;
    }
    return copy;
  }

  //=========================================================================//
  // "MechElement::GetECOSTS":                                               //
  //=========================================================================//
  template<LVSC  LVSCKind>
  constexpr TT MechElement<LVSCKind>::GetECOSTS() const
  {
    // ECOS TSs of all "Vector3D"s must be identical. Check that, and return
    // the common TS:
    TT ts0 = m_CoM    .GetCOSTS();
    DEBUG_ONLY
    (
    TT ts1 = m_CoMDots.GetCOSTS();
    TT ts2 = m_MoIs   .GetCOSTS();
    TT ts3 = m_MoIDots.GetCOSTS();
    )
    assert
      ((ts0.IsUnDef() &&  ts1.IsUnDef() &&  ts2.IsUnDef() && ts3.IsUnDef()) ||
       ((ts0 == ts1)  && (ts1 == ts2)   && (ts2 == ts3)));
    return ts0;
  }

  //=========================================================================//
  // "MechElement::UnifyECOSTSs":                                            //
  //=========================================================================//
  template<LVSC  LVSCKind>
  constexpr void MechElement<LVSCKind>::UnifyECOSTSs
    (MechElement const& a_right)
  {
    // Get the Common TimeStamps of the LHS and the RHS "Vector3D"s:
    TT tsL = GetECOSTS();
    TT tsR = a_right.GetECOSTS();

    // They must be compatible:
    DEBUG_ONLY(bool hasUnDef =  tsL.IsUnDef() || tsR.IsUnDef();)
    assert(hasUnDef || tsL == tsR);

    // Update all TimeStamps of the LHS, if the UnifiedTS is different from the
    // existing "tsL". This can happen only if "tsL" was NOT "UnDef", but "tsR"
    // is "UnDef", so the Unified TS becomes "UnDef" as well, and we need to in-
    // stall in in the LHS "Vector3D"s:
    //
    if (!tsL.IsUnDef() && tsR.IsUnDef())
    {
      m_CoM    .GetCOSTS() = TT::UnDef();
      m_CoMDots.GetCOSTS() = TT::UnDef();
      m_MoIs   .GetCOSTS() = TT::UnDef();
      m_MoIDots.GetCOSTS() = TT::UnDef();
    }
  }

  //=========================================================================//
  // Addition / Subtraction of "MechElements"s:                              //
  //=========================================================================//
  // The operands must both have FINAL masses. Only one summand is allowed to
  // have zero mass, otherwise we cannot compute the CoM.  WE ASSUME that the
  // summands DO NOT INTERSECT in space, and for subtraction, the 2nd operand
  // entirely belongs to the 1st one:
  //-------------------------------------------------------------------------//
  // "+="                                                                    //
  //-------------------------------------------------------------------------//
  template<LVSC  LVSCKind>
  constexpr MechElement<LVSCKind>& MechElement<LVSCKind>::operator+=
    (MechElement const& a_right)
  {
    assert(m_isFinal && a_right.m_isFinal &&
           !(IsZero(m_mass) && IsZero(a_right.m_mass)));

    // Check  that the ECOS TS of the LHS and the RHS Vectors are consistent,
    // and install the unified TS:
    UnifyECOSTSs(a_right);

    // Memoise the orig Mass and MassDot -- required later:
    Mass     m0    = m_mass;
    MassRate m0dot = m_massDot;

    // Masses, SurfAreas, Vols and MoIs are directly-additive:
    m_mass        += a_right.m_mass;
    m_massDot     += a_right.m_massDot;
    m_enclVol     += a_right.m_enclVol;

    // For the CoM and CoMDots, do the weighted avg (but the total mass must be
    // non-0):
    assert(IsPos(m_mass));
    double mu0  = double(m0 / m_mass);
    assert(0.0 <= mu0 && mu0 <= 1.0);
    double mu1  = 1.0 -  mu0;

    // XXX: BEWARE that "mu0" and "mu1" may themselves be functions of time,
    // which must be taken into account in "CoMDots" computation; obviously,
    // mu1dot = - mu0dot:
    auto mu0dot = (m0dot - mu0 * m_massDot) / m_mass;

    for (int i = 0; i < 3; ++i)
    {
      m_MoIs   [i] += a_right.m_MoIs   [i];
      m_MoIDots[i] += a_right.m_MoIDots[i];

      // BEWARE: Modify "CoMDots" first, because they required the OLD vals
      // of "CoM":
      m_CoMDots[i]  = mu0 * m_CoMDots  [i] + mu1 * a_right.m_CoMDots[i] +
                      mu0dot   * (m_CoM[i] -       a_right.m_CoM    [i]);
      m_CoM    [i]  = mu0 * m_CoM      [i] + mu1 * a_right.m_CoM    [i];
    }
    return *this;
  }

  //-------------------------------------------------------------------------//
  // "+":                                                                    //
  //-------------------------------------------------------------------------//
  template<LVSC  LVSCKind>
  constexpr MechElement<LVSCKind> MechElement<LVSCKind>::operator+
    (MechElement const& a_right)
  const
  {
    MechElement res = *this;
    res += a_right;
    return res;
  }

  //-------------------------------------------------------------------------//
  // "-=":                                                                   //
  //-------------------------------------------------------------------------//
  template<LVSC  LVSCKind>
  constexpr MechElement<LVSCKind>& MechElement<LVSCKind>::operator-=
    (MechElement const& a_right)
  {
    assert(m_isFinal && a_right.m_isFinal &&
           !(IsZero(m_mass) && IsZero(a_right.m_mass)));

    // Check  that the ECOS TS of the LHS and the RHS Vectors are consistent,
    // and install the unified TS:
    UnifyECOSTSs(a_right);

    // Memoise the orig Mass and MassDot -- required later:
    Mass     m0    = m_mass;
    MassRate m0dot = m_massDot;

    // Masses, SurfAreas, Vols and MoIs are directly-additive:
    m_mass        -= a_right.m_mass;
    m_massDot     -= a_right.m_massDot;
    m_enclVol     -= a_right.m_enclVol;

    // For the CoM and CoMDots, do the weighted avg (but the total mass and
    // volume must be > 0):
    assert(IsPos(m_mass) && IsPos(m_enclVol));
    double mu0    = double(m0             / m_mass);
    double mu1    = double(a_right.m_mass / m_mass);

    // XXX: BEWARE that "mu0" and "mu1" may themselves be functions of time,
    // which must be taken into account in "CoMDots" computation:
    //
    auto mu0dot = (m0dot             - mu0 * m_massDot) / m_mass;
    auto mu1dot = (a_right.m_massDot - mu1 * m_massDot) / m_mass;

    for (int i = 0; i < 3; ++i)
    {
      m_MoIs   [i] -= a_right.m_MoIs    [i];
      m_MoIDots[i] -= a_right.m_MoIDots [i];

      // Again, modify "CoMDots" first, because they require the OLD vals of
      // "CoM":
      m_CoMDots[i]  = mu0    * m_CoMDots[i] - mu1    * a_right.m_CoMDots[i] +
                      mu0dot * m_CoM    [i] - mu1dot * a_right.m_CoM    [i];
      m_CoM    [i]  = mu0    * m_CoM    [i] - mu1    * a_right.m_CoM    [i];
    }
    return *this;
  }

  //-------------------------------------------------------------------------//
  // "-":                                                                    //
  //-------------------------------------------------------------------------//
  template<LVSC  LVSCKind>
  constexpr MechElement<LVSCKind> MechElement<LVSCKind>::operator-
    (MechElement const& a_right)
  const
  {
    MechElement res = *this;
    res -= a_right;
    return res;
  }

  //=========================================================================//
  // "PointMass" Non-Default Ctor:                                           //
  //=========================================================================//
  template<LVSC LVSCKind>
  constexpr PointMass<LVSCKind>::PointMass
  (
    TT   a_ecos_ts,
    Len  a_x0,
    Len  a_y0,
    Len  a_z0,
    Mass a_mass
  )
  : MechElement<LVSCKind>()     // Slightly sub-optimal...
  {
    Len  pt[3] { a_x0, a_y0, a_z0 };

    // Compute the MoIs:
    Len2 x2 = Sqr(a_x0);
    Len2 y2 = Sqr(a_y0);
    Len2 z2 = Sqr(a_z0);
    assert(IsPos(a_mass));

    MoI  mois[3]
    {
      a_mass * (y2 + z2), // Distance^2 to the OX axis
      a_mass * (x2 + z2), // Distance^2 to the OY axis
      a_mass * (x2 + y2), // Distance^2 to the OZ axis
    };

    // Now the actual Base Class initialisation. NB: The point-mass is consi-
    // dered to be Final, all Rates are 0:
    constexpr MoIRate Rates0[3] {MoIRate(0.0), MoIRate(0.0), MoIRate (0.0)};
    constexpr Vel     Vel0  [3] {Vel    (0.0), Vel    (0.0), Vel     (0.0)};

    MechElement<LVSCKind>::Init
      (this,   a_ecos_ts, pt, Vel0, a_mass, MassRate(0.0), Vol(0.0), mois,
       Rates0, true);
  }

  //=========================================================================//
  // "RotationShell::Init":                                                  //
  //=========================================================================//
  // Similar to "MechElement::Init": Initialiser similar to Non-Default Ctor
  // but used in a slightly different way:
  //
  template<LVSC LVSCKind,       typename Derived>
  constexpr void RotationShell<LVSCKind, Derived>::Init
  (
    // Params for the Base Class ("MechElement"):
    TT      a_ecos_ts,
    Area    a_side_surf_area,
    Vol     a_encl_vol,
    Mass    a_mass,     // If 0, then auto-calculated with SurfDens=1

    // Params for "RotationShell" itself. Providing both Cos and Sin of Alpha
    // and Psi is for optimisation only,  because these angles are typically
    // the same for several "ME"s. Also, (x0, y0, z0) are not independent of
    // Psi but provided for convenience:
    double  a_cosA,     // cos(alpha)
    double  a_sinA,     // sin(alpha)
    double  a_cosP,     // cos(psi)
    double  a_sinP,     // sin(psi)
    Len     a_x0,       // The Up or Low end of the Xi axis
    Len     a_y0,       //
    Len     a_z0,       //
    bool    a_0is_up,   // So is (x0,y0,z0) the Up or Low end?
    Len     a_h,        // Over-all length  (along the rotation axis)

    // "Empty" MoI Coeffs (wrt the Low end of the rotation axis):
    Len4    a_je0,
    Len4    a_je1,
    Len3    a_ke,

    // Propellant Density (0 if no Propellant is held in this Body):
    Density a_rho
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    // 0 <= alpha < Pi/2, so:
    assert(a_cosA > 0.0 && a_sinA >= 0.0);

    // If a_y0==a_z0==0,   ie there is a point on the Xi rotation axis segment
    // which also lies on the OX axis, then we assume both axes should coinci-
    // de:
    assert(!(IsZero(a_y0)  && IsZero(a_z0)) ||
            (a_cosA == 1.0 && a_sinA == 0.0));

    double      cosA2 = Sqr(a_cosA);
    double      sinA2 = Sqr(a_sinA);
    double      cosP2 = Sqr(a_cosP);
    double      sinP2 = Sqr(a_sinP);
    assert(ApproxEqual(cosA2 + sinA2, 1.0) &&
           ApproxEqual(cosP2 + sinP2, 1.0));

    // Cos(Psi), Sin(Psi) may be of any sign, but they are actually determined
    // by "a_y0", "a_z0"; they are passed explicitly only for the sake of opt-
    // imisation and accuracy:
    DEBUG_ONLY
    (
      Len2   y02  = Sqr(a_y0);
      Len2   z02  = Sqr(a_z0);
      Len2   u02  = y02 + z02;
      assert((u02 * cosP2).ApproxEquals(y02) &&
             (u02 * sinP2).ApproxEquals(z02));
    )
    assert(IsPos(a_h) && IsPos(a_je0) && IsPos(a_je1) && IsPos(a_ke) &&
          !IsNeg(a_rho));
    // NB:
    // (*) In general, "a_ke" may be of any sign, but because it is a moment
    //     wrt the lower Xi end,  in this case it should be positive as well;
    // (*) a_rho==0 is OK if there is no Propellant in this Shell...

    //-----------------------------------------------------------------------//
    // Initialise the "RotationShell" flds FIRST:                            //
    //-----------------------------------------------------------------------//
    // Directing Cosines of the Xi axis (Low -> Up):
    m_h        =  a_h;
    m_xi[0]    =  a_cosA;
    m_xi[1]    = -a_sinA * a_cosP;
    m_xi[2]    = -a_sinA * a_sinP;

    // The Up (Larger-X) and Low (Smaller-X) ends of the rotation axis:
    if (a_0is_up)
    {
      m_up [0] = a_x0;
      m_up [1] = a_y0;
      m_up [2] = a_z0;
      m_low[0] = a_x0 - m_h * m_xi[0];
      m_low[1] = a_y0 - m_h * m_xi[1];
      m_low[2] = a_z0 - m_h * m_xi[2];
    }
    else
    {
      m_low[0] = a_x0;
      m_low[1] = a_y0;
      m_low[2] = a_z0;
      m_up [0] = a_x0 + m_h * m_xi[0];
      m_up [1] = a_y0 + m_h * m_xi[1];
      m_up [2] = a_z0 + m_h * m_xi[2];
    }
    assert((Sqr(m_up[0] - m_low[0]) + Sqr(m_up[1] - m_low[1]) +
            Sqr(m_up[2] - m_low[2])).ApproxEquals(Sqr(m_h)));
    assert(m_low[0] < m_up[0]);

    // Side Surace Area and Nominal Volume Enclosed:
    m_sideSurfArea = a_side_surf_area;
    assert(IsPos(m_sideSurfArea) && IsPos(a_encl_vol));

    // Propellant Params:
    m_rho          = a_rho;
    m_propMassCap  = m_rho * a_encl_vol;

    // IMPORTANT!
    // Coeffs for (Jx, Jy, Jz) wrt (J0, J1, K, SurfOrVol), where J0, J1, K
    // are computed relative to the Low end-point of the rotation axis (Xi).
    // This is just a matter of choice for 2D Shells, but becomes important
    // for 3D Propellant Volumes: as Propellants are spent, the Low end of
    // of the Propellant Volume remains unchanged, whereas the Up one moves
    // with the decreasing level of Propellant:
    //
    Len  xL  = m_low[0];
    Len  yL  = m_low[1];
    Len  zL  = m_low[2];
    Len2 xL2 = Sqr(xL);
    Len2 yL2 = Sqr(yL);
    Len2 zL2 = Sqr(zL);

    m_Jx0    =  sinA2;
    m_Jx1    =  1.0 + cosA2;
    m_JxK    = -2.0 * a_sinA  * (a_cosP * yL + a_sinP * zL);
    m_JxSV   =  yL2 + zL2;

    m_Jy0    =  1.0 - sinA2   * cosP2;
    m_Jy1    =  1.0 + sinA2   * cosP2;
    m_JyK    =  2.0 * (a_cosA * xL  - a_sinA * a_sinP * zL);
    m_JySV   =  xL2 + zL2;

    m_Jz0    =  1.0 - sinA2   * sinP2;
    m_Jz1    =  1.0 + sinA2   * sinP2;
    m_JzK    =  2.0 * (a_cosA * xL  - a_sinA * a_cosP * yL);
    m_JzSV   =  xL2 + yL2;

    //-----------------------------------------------------------------------//
    // Now the Base Class ("MechElement"):                                   //
    //-----------------------------------------------------------------------//
    // The mass may or may not be given. If it is given,  calculate the Surf-
    // ace Density; otherwise, assume the SurfDens to be 1.0.   Then set the
    // Mass and MoIs:
    bool     isFinal  = IsPos(a_mass);
    SurfDens surfDens =
      isFinal ? (a_mass /  a_side_surf_area) : SurfDens(1.0);
    Mass emptyMass    =
      isFinal ?  a_mass : (a_side_surf_area  * surfDens);

    // "Empty" CoM and MoIs (2D) for the Base Class:
    Len       emptyCoM    [3];
    Vel       emptyCoMDots[3];
    MoI       emptyMoIs   [3];
    MoIRate   emptyMoIDots[3];
    constexpr Len4Rate JDot0 = Len4Rate(0.0);
    constexpr Len3Rate KDot0 = Len3Rate(0.0);
    constexpr Len2Rate SDot0 = Len2Rate(0.0);
    constexpr MassRate MDot0 = MassRate(0.0);
    MoIsCoM
    (
      a_je0, a_je1,    a_ke,     JDot0, JDot0, KDot0,     a_side_surf_area,
      SDot0, surfDens, emptyCoM, emptyCoMDots, emptyMoIs, emptyMoIDots
    );
    // CoMDots and MoIDots must all be 0s, of course:
    assert(IsZero(emptyCoMDots[0]) && IsZero(emptyCoMDots[1]) &&
           IsZero(emptyCoMDots[2]) && IsZero(emptyMoIDots[0]) &&
           IsZero(emptyMoIDots[1]) && IsZero(emptyMoIDots[2]));

    // Finally, invoke the Base Class Initialiser:
    MechElement<LVSCKind>::Init
       (this,       a_ecos_ts, emptyCoM,     emptyCoMDots, emptyMass, MDot0,
        a_encl_vol, emptyMoIs, emptyMoIDots, isFinal);
  }

  //=========================================================================//
  // "RotationShell::MoIsCoM":                                               //
  // Computation of XYZ MoIs and CoM from "intrinsic" MoIs:                  //
  //=========================================================================//
  // The same function is applicable to both MoI per SurfDens (Len4) and MoI
  // per (Volume) Density (Len5), hence:
  // J=Len4, K=Len3, SV=Len2, OR
  // J=Len5, K=Len4, SV=Len3:
  //
  template<LVSC LVSCKind,       typename Derived>
  template<typename J>
  constexpr void RotationShell<LVSCKind, Derived>::MoIsCoM
  (
    J                          a_j0,      // Len4 or Len5, wrt Low axis end
    J                          a_j1,      // ditto
    decltype(a_j0  /1.0_m)     a_k,       // Len3 or Len4
    decltype(a_j0  /1.0_sec)   a_j0_dot,  // Time derivatives of J0, J1, K
    decltype(a_j0  /1.0_sec)   a_j1_dot,  //
    decltype(a_k   /1.0_sec)   a_k_dot,   //
    decltype(a_j0  /Len2(1.0)) a_sv,      // Len2 or Len3 (ie SurfArea or Vol)
    decltype(a_sv  /1.0_sec)   a_sv_dot,  // SurfAreaDot: 0; VolDot: <= 0
    decltype(1.0_kg/a_sv)      a_dens,    // SurfDens or Density
    Len                        a_com     [3],    // Output
    Vel                        a_com_dots[3],    //
    MoI                        a_mois    [3],    //
    MoIRate                    a_moi_dots[3]     //
  )
  const
  {
    // "Intrinsic" MoI components must be > 0 for any finite body size,  as
    // well as "a_k" (because it's a moment relative to the Low (Smallest-X)
    // point). On the contrary, their time derivatives must be <= 0:
    assert( IsPos(a_j0) && IsPos(a_j1) && IsPos(a_sv) && IsPos(a_k) &&
            IsPos(a_sv) &&
          !(IsPos(a_j0_dot) || IsPos(a_j1_dot) || IsPos(a_k_dot)) );

    // If it's about a Surface,  all Rates must be 0;
    // otherwise, it's a Volume, and all Rates must be <= 0:
    assert
      ((std::is_same_v<J, Len4> &&  IsZero(a_sv_dot)  &&
        IsZero(a_j0_dot)        &&  IsZero(a_j1_dot)  && IsZero(a_k_dot)) ||
       (std::is_same_v<J, Len5> && !(IsPos(a_sv_dot)  || IsPos (a_j0_dot) ||
        IsPos(a_j1_dot)         ||   IsPos (a_k_dot))));

    // MoIs:
    J Jx = m_Jx0 * a_j0 + m_Jx1 * a_j1 + m_JxK * a_k + m_JxSV * a_sv;
    J Jy = m_Jy0 * a_j0 + m_Jy1 * a_j1 + m_JyK * a_k + m_JySV * a_sv;
    J Jz = m_Jz0 * a_j0 + m_Jz1 * a_j1 + m_JzK * a_k + m_JzSV * a_sv;

    // MoI Rates: Linear transforms similar to above:
    using JRate  = decltype(a_j0 / 1.0_sec);
    JRate JxDot  =
      m_Jx0   * a_j0_dot + m_Jx1 * a_j1_dot + m_JxK * a_k_dot +
      m_JxSV  * a_sv_dot;
    JRate JyDot   =
      m_Jy0   * a_j0_dot + m_Jy1 * a_j1_dot + m_JyK * a_k_dot +
      m_JySV  * a_sv_dot;
    JRate JzDot   =
      m_Jz0   * a_j0_dot + m_Jz1 * a_j1_dot + m_JzK * a_k_dot +
      m_JzSV  * a_sv_dot;

    a_mois[0]     = a_dens * Jx;
    a_mois[1]     = a_dens * Jy;
    a_mois[2]     = a_dens * Jz;
    assert(!(IsNeg(a_mois[0]) || IsNeg(a_mois[1]) || IsNeg(a_mois[2])));

    a_moi_dots[0] = a_dens * JxDot;
    a_moi_dots[1] = a_dens * JyDot;
    a_moi_dots[2] = a_dens * JzDot;

    assert((std::is_same_v<J, Len4> &&  IsZero(a_moi_dots[0])  &&
            IsZero(a_moi_dots[1])   &&  IsZero(a_moi_dots[2])) ||
           (std::is_same_v<J, Len5> && !(IsPos(a_moi_dots[0])  ||
            IsPos(a_moi_dots[1])    ||   IsPos(a_moi_dots[2]))));

    // Now the CoM: "xiCoM" is its co-ord along the rotation axis, relative
    // to the Low (Smallest-X) axis end, hence positive:
    Len  xiCoM    = a_k / a_sv;
    assert(IsPos(xiCoM));
    a_com[0]      = m_low[0] + m_xi[0] * xiCoM;
    a_com[1]      = m_low[1] + m_xi[1] * xiCoM;
    a_com[2]      = m_low[2] + m_xi[2] * xiCoM;

    // Similar for CoMDots:
    Vel xiCoMDot  = a_k_dot / a_sv - a_k * a_sv_dot / Sqr(a_sv);
    a_com_dots[0] = m_xi[0]  * xiCoMDot;
    a_com_dots[1] = m_xi[1]  * xiCoMDot;
    a_com_dots[2] = m_xi[2]  * xiCoMDot;
  }

  //=========================================================================//
  // "RotationShell::GetPropBulkME":                                         //
  //=========================================================================//
  // Constructs a "MechElement" which holds the CoM and MoI of the Propellant
  // Bulk filling this RotationShell, with the current propellant mass given
  // by "a_prop_mass". The resulting "MechElement" does NOT include the Shell!
  // Returns a ficticious "MechElement" obj  acting as a container for the
  // results computed, and suitable for the "+" operation. May also return
  // the curr Propellant Level via the 2nd arg (if non-NULL), primarily for
  // debugging purposes.
  // Although declared "constexpr", this function will mostly be called at
  // Run-Time. For the info, it also returns the PropLevel if the output ptr
  // is non-NULL:
  //
  template<LVSC LVSCKind, typename Derived>
  constexpr MechElement<LVSCKind>
  RotationShell<LVSCKind, Derived>::GetPropBulkME
  (
    TT       a_ecos_ts,        // Usually DEFINED at Run-Time
    Mass     a_prop_mass,      // >= 0; +oo means "FullPropellantLoad"
    MassRate a_prop_mass_dot,  // Must be <= 0 in general
    Len*     a_prop_level      // May be NULL, then not computed
  )
  const
  {
    // The Propellant Level, relative to the Low (Smallest-X) Base,  and its
    // time derivative. XXX: We assume that the propellant surface is always
    // orthogonal to the rotation axis, due to the tank pressurisation:
    Len l      (NaN<double>);
    Vel lDot   (NaN<double>);
    Vol propVol(NaN<double>);

    VolRate propVolDot = a_prop_mass_dot / GetPropDens();
    assert(!IsPos(propVolDot));

    if (!IsFinite(a_prop_mass) && IsZero(a_prop_mass_dot))
    {
      // This is a special "static" case, Full Propellant Load:
      propVol = ME::GetEnclVol();
      l       = m_h;
      lDot    = Vel(0.0);
    }
    else
    {
      // "a_prop_mass" is explicitly specified. Check its limits. For the 2nd
      // inequaluty, allow some floating-point tolerance:
      assert(!IsNeg(a_prop_mass) &&
             a_prop_mass <= GetPropMassCap() * TolFact);

      // The Propellant Volume: Make sure it is within the limits to avoid
      // rounding errors:
      propVol = a_prop_mass  /  GetPropDens();
      assert(!IsNeg(propVol) && propVol <= ME::GetEnclVol() * TolFact);
      propVol = std::min(propVol, ME::GetEnclVol()); // Enforce it for safety

      // Solve for (l, lDot) using the "Derived" geometry:
      auto lld =
        static_cast<Derived const*>(this)->Derived::PropLevelOfVol
          (propVol, propVolDot);
      l    = lld.first;
      lDot = lld.second;

      // Mathematically, we must have 0 <= l <= h,   where l=0 corresponds to
      // the empty Element, and l=h to the full one. We enforce this interval
      // explicitly to prevent rounding errors:
      assert(!IsNeg(l) && l <= m_h * TolFact);
      l        = std::min(l,   m_h);
    }
    assert(IsFinite(propVol) && !IsNeg(propVol));

    // Memoise the "l" if required:
    if (a_prop_level != nullptr)
       *a_prop_level  = l;

    // If an infinite (placeholder) "a_prop_mass" was provided, restore its
    // actual value:
    if (!IsFinite(a_prop_mass))
      a_prop_mass = GetPropMassCap();

    // Get the MoIs and the CoM of the Popellant. NB: the "Intrinsic" MoI
    // components of the Propellant are computed by "Derived":
    Len5     JP0,    JP1;
    Len4     KP;
    Len5Rate JP0Dot, JP1Dot;
    Len4Rate KPDot;

    static_cast<Derived const*>(this)->Derived::PropMoIComps
      (l, lDot, &JP0, &JP1, &KP, &JP0Dot, &JP1Dot, &KPDot);

    Len     com    [3];
    Vel     comDots[3];
    MoI     mois   [3];
    MoIRate moiDots[3];

    MoIsCoM
      (JP0, JP1,  KP, JP0Dot, JP1Dot,  KPDot, propVol, propVolDot,
       GetPropDens(), com,    comDots, mois,  moiDots);

    // Construct the result (XXX: this is slightly sub-optimal, as involves
    // copying of "com", "mois" and their "Dots").
    // NB:
    // (*) We must  install the correct Mass of the Propellant in the "res",
    //     so that  it can then participate in "operator+";
    // (*) SurfArea is not computed and remains 0;, we do not compute it;
    //     this is OK because by SurfArea we mean that of the Shell (with
    //     some non-0 SurfDens), whereas in this case, the Shell is not
    //     included in the result;
    // (*) "a_ecos_ts" from the CallER is installed in the results;
    // (*) the result is Final:
    //
    return MechElement<LVSCKind>
      (a_ecos_ts, com, comDots, a_prop_mass, a_prop_mass_dot, propVol, mois,
       moiDots,   true);
  }
}
// End namespace SpaceBallistics

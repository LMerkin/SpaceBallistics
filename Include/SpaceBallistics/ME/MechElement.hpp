// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/ME/MechElement.hpp":                 //
//                   Geometric Objects as "Mechanical Elements"              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/LVSC/LVSC.h"
#include "SpaceBallistics/CoOrds/EmbeddedCOS.h"
#include <boost/container/static_vector.hpp>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <initializer_list>

namespace SpaceBallistics
{
  //=========================================================================//
  // "MechElement": Base for Constant-or-Variable Mass Mechanical Elements:  //
  //=========================================================================//
  // General characteristics of CONSTANT-or-VARIABLE-MASS MECHANICAL ELEMENTS:
  // Currently, 2D Shells (Truncated Cones and Spherical Segments), 3D Propel-
  // lant Bulks (of the same shape as Shells)  and Point Masses are supported,
  // via the corresp Derived Classes of "MechElement".  For each them, we can
  // compute the Masses, Centers of Masses (CoM) and Moments of Inertia (MoI)
  // wrt  OX, OY, OZ axes,   as well as time derivatives ("Dots") of CoMs and
  // MoIs.
  // NB:
  // It is assumed that the co-ords system OXYZ is the Embedded CS  of the obj-
  // ect given by the "LVSKind" template param. Then OX is the object's princi-
  // pal axis of symmetry, with the positive direction pointing UPWARDS (ie
  // from the Tail to the Nose, normally making a SHARP angle with the velocity
  // vector). This is important when it comes  to the CoM and MoIs of the cont-
  // ained Propellant: We assume that the Propellant is concentrated at the
  // bottom (Lower, Smaller-X) part of the "MechElement", due to the gravity and
  // tanks pressurisation (FIXME what about the SCs in the Free Space???):
  //
  template<LVSC LVSCKind>
  class MechElement
  {
  public:
    //=======================================================================//
    // Types and Consts:                                                     //
    //=======================================================================//
    // The Embedded CoOrds System for this "LVSCKind":
    using ECOS      = EmbeddedCOS<LVSCKind>;

    // Position Vector in the Embedded COS:
    using PosVE     = PosV    <ECOS>;

    // Velocity in the Embedded COS (XXX: here it is only used for the motion
    // of the CoM, so not very important):
    using VelVE     = VelV    <ECOS>;

    // Force in the Embedded COS:
    using ForceVE   = ForceV  <ECOS>;

    // The MoI Tensor (XXX: currently the main diagonal only), also in the ECOS:
    using MoITE     = MoIT    <ECOS>;

    // The Time Rate of the MoI Tensor:
    using MoIRateTE = MoIRateT<ECOS>;

    // "UnKnownMass": Param to be used for "MechElement"s when their mass is
    // not yet known (all real Masses are of course strictly positive):
    //
    constexpr static Mass UnKnownMass = 0.0_kg;

  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Over-All Dynamical Params of this "MechElement":
    //
    PosVE       m_CoM;         // (X,Y,Z) co-ords of the Center of Masses
    VelVE       m_CoMDots;     // Velocity of the moving CoM
    Mass        m_mass;        // Mass
    MassRate    m_massDot;     // d(Mass)/dt, typically <= 0
    MoITE       m_MoIs;        // Moments of Inertia wrt the OX, OY and OZ axes
    MoIRateTE   m_MoIDots;     // d(MoIs)/dt
    bool        m_isFinal;     // False => Mass, MoIs & their Dots not valid
                               //   yet (but the CoM is valid!)
  protected:
    //=======================================================================//
    // "Init":                                                               //
    //=======================================================================//
    // For use by Derived Classes. It is similar to a Non-Default Ctor, but is
    // intended to be invoked AFTER the Derived flds have been initialised, so
    // implemented as a Static Method, not a Ctor.
    // XXX: This method is made "static"; otherwise, it could not be called from
    // methods of derived classes on OTHER objs, even though it is "protected":
    // XXX: The params are just arrays, not ECOS-parameterised vectors/tensors:
    // since this function is internal, that is OK:
    //
    constexpr static void Init
    (
      MechElement*    a_me,
      Len  const      a_com     [3],
      Vel  const      a_com_dots[3],
      Mass            a_mass,
      MassRate        a_mass_dot,
      MoI     const   a_mois    [3],
      MoIRate const   a_moi_dots[3],
      bool            a_is_final
    )
    {
      // NB: Even if "a_mass" is not final, it must be positive. In the context
      // under consideration, MassDot and MoIDots cannot be positive, since the
      // ContrElement's mass is either constant or being exhausted:
      //
      assert(a_me != nullptr      &&
             IsPos(a_mass)        && !IsNeg(a_mois    [0]) &&
            !IsNeg(a_mois    [1]) && !IsNeg(a_mois    [2]) &&
            !IsPos(a_mass_dot)    && !IsPos(a_moi_dots[0]) &&
            !IsPos(a_moi_dots[1]) && !IsPos(a_moi_dots[1]));

      a_me->m_CoM[0]     = a_com     [0];
      a_me->m_CoM[1]     = a_com     [1];
      a_me->m_CoM[2]     = a_com     [2];

      a_me->m_CoMDots[0] = a_com_dots[0];
      a_me->m_CoMDots[1] = a_com_dots[1];
      a_me->m_CoMDots[2] = a_com_dots[2];

      a_me->m_mass       = a_mass;
      a_me->m_massDot    = a_mass_dot;

      a_me->m_MoIs[0]    = a_mois    [0];
      a_me->m_MoIs[1]    = a_mois    [1];
      a_me->m_MoIs[2]    = a_mois    [2];

      a_me->m_MoIDots[0] = a_moi_dots[0];
      a_me->m_MoIDots[1] = a_moi_dots[1];
      a_me->m_MoIDots[2] = a_moi_dots[2];

      a_me->m_isFinal    = a_is_final;
    }

  public:
    //=======================================================================//
    // Default Ctor:                                                         //
    //=======================================================================//
    // NB:
    // (*) All numeric fields are initialised to 0s, not to NaNs, so the "empty"
    //     obj (constructed by the Default Ctor) can be used as an initial value
    //     for summation ("+", "+=", etc);
    // (*) "IsFinal" is set to "true", because the empty "MechElement" is typic-
    //     ally used as the base for "+" which requires Final operands:
    //
    constexpr MechElement()
    : m_CoM    {{0.0_m,  0.0_m,        0.0_m}},
      m_CoMDots{{Vel    (0.0), Vel    (0.0), Vel    (0.0)}},
      m_mass   (0.0),
      m_massDot(0.0),
      m_MoIs   {{MoI    (0.0), MoI    (0.0), MoI    (0.0)}},
      m_MoIDots{{MoIRate(0.0), MoIRate(0.0), MoIRate(0.0)}},
      m_isFinal(true)
    {}

    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr MechElement
    (
      Len      const a_com     [3],
      Vel      const a_com_dots[3],
      Mass           a_mass,
      MassRate const a_mass_dot,
      MoI            a_mois    [3],
      MoIRate  const a_moi_dots[3],
      bool           a_is_final
    )
    {
      Init(this, a_com, a_com_dots, a_mass, a_mass_dot, a_mois, a_moi_dots,
           a_is_final);
    }
    // Copy Ctor, Assignment and Equality are auto-generated...

    //=======================================================================//
    // "GetMassScale":                                                       //
    //=======================================================================//
    // Given a list of "MechElement"s  (all of them must have Non-Final Masses)
    // and the real Total Mass of them,  returns  a dimension-less factor which
    // could be applied (in "ProRateMass") to each "MechElement" to set its cor-
    // rect Final Mass. It assumes that the densities  (surface  or vol) of all
    // "MechElement"s in the list are the same, ie their relative masses remain
    // unchanged after scaling:
    //
    constexpr static double GetMassScale
      (std::initializer_list<MechElement const*> a_mes, Mass a_total_mass)
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

    //=======================================================================//
    // "ProRateMass":                                                        //
    //=======================================================================//
    // Returns a new object (of any type derived from "MechElement") which diff-
    // ers from the original one  by the Mass and MoIs being multiplied by the
    // ScaleFactor.
    // XXX:
    // (*) it is assumed that any other flds of "Derived" are UNAFFECTED by this
    //     scaling;
    // (*) we also assume that this operation is only applicable to Fixed-Mass
    //     "MechElement"s, so all Rates must be 0, and are not affected by the
    //     scaling:
    //
    template<typename Derived>
    constexpr static  Derived ProRateMass(Derived const& a_der, double a_scale)
    {
      static_assert(std::is_base_of_v<MechElement, Derived>);
      assert(a_scale > 0.0);

      // Create a copy of "a_der" using the Copy Ctor of Derived (which must
      // indeed be derived from "MechElement"):
      Derived copy(a_der);

      // Cannot adjust the Mass once it has been finalised, or if any Rates are
      // non-0:
      assert(!copy.m_isFinal           && IsZero(copy.m_massDot)    &&
             IsZero(copy.m_MoIDots[0]) && IsZero(copy.m_MoIDots[1]) &&
             IsZero(copy.m_MoIDots[2]));

      // Set the Mass and adjust the MoIs; but the CoM is unchanged:
      copy.m_mass    *= a_scale;
      copy.m_MoIs[0] *= a_scale;
      copy.m_MoIs[1] *= a_scale;
      copy.m_MoIs[2] *= a_scale;
      copy.m_isFinal  = true;
      return copy;
    }

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    // Geometrical params are always available, even if the Mass is not final
    // yet:
    constexpr PosVE const& GetCoM() const { return m_CoM; }

    // But the Mass and MoIs  may or may not be available.
    // If not, assert failure will be signaled (NB: in C++ >= 14, "assert" is
    // allowed in "constexpr" functions, whereas throwing exceptions is not):
    //
    constexpr bool IsFinal() const { return m_isFinal; }

    constexpr Mass GetMass() const
    {
      assert(m_isFinal);
      return m_mass;
    }

    constexpr MassRate         GetMassDot() const
    {
      assert(m_isFinal);
      return m_massDot;
    }

    constexpr MoITE     const& GetMoIs()    const
    {
      assert(m_isFinal);
      return m_MoIs;
    }

    constexpr MoIRateTE const& GetMoIDots() const
    {
      assert(m_isFinal);
      return m_MoIDots;
    }

    //=======================================================================//
    // Addition / Subtraction:                                               //
    //=======================================================================//
    // The summands must both have FINAL masses. Only one summand is allowed to
    // have zero mass, otherwise we cannot compute the CoM.  WE ASSUME that the
    // summands DO NOT INTERSECT in space:
    //
    constexpr MechElement& operator+= (MechElement const& a_right)
    {
      assert(m_isFinal && a_right.m_isFinal &&
             !(IsZero(m_mass) && IsZero(a_right.m_mass)));

      // Masses, SurfAreas, Vols and MoIs are directly-additive:
      Mass      m0  = m_mass;
      m_mass       += a_right.m_mass;
      m_massDot    += a_right.m_massDot;

      m_MoIs   [0] += a_right.m_MoIs   [0];
      m_MoIs   [1] += a_right.m_MoIs   [1];
      m_MoIs   [2] += a_right.m_MoIs   [2];

      m_MoIDots[0] += a_right.m_MoIDots[0];
      m_MoIDots[1] += a_right.m_MoIDots[1];
      m_MoIDots[2] += a_right.m_MoIDots[2];

      // For the CoM, do the weighted avg (but the total mass must be non-0):
      assert(IsPos(m_mass));
      double mu0  = double(m0             / m_mass);
      double mu1  = double(a_right.m_mass / m_mass);
      m_CoM[0]    = mu0 * m_CoM[0]  + mu1 * a_right.m_CoM[0];
      m_CoM[1]    = mu0 * m_CoM[1]  + mu1 * a_right.m_CoM[1];
      m_CoM[2]    = mu0 * m_CoM[2]  + mu1 * a_right.m_CoM[2];
      return *this;
    }

    constexpr MechElement operator+ (MechElement const& a_right) const
    {
      MechElement res = *this;
      res += a_right;
      return res;
    }
  };

  //=========================================================================//
  // "PointMass" Class:                                                      //
  //=========================================================================//
  // A positive mass concentrated in the (x0, y0, z0) point:
  //
  template<LVSC LVSCKind>
  class PointMass final: public MechElement<LVSCKind>
  {
  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    constexpr PointMass(Len a_x0, Len a_y0, Len a_z0, Mass a_mass)
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
        (this, pt, Vel0, a_mass, MassRate(0.0), mois, Rates0, true);
    }
  };

  //=========================================================================//
  // "RotationShell" Class:                                                  //
  //=========================================================================//
  // SubClass of "MechElement", SuperClass of "TrCone", "SpherSegment" and
  // others. Provides common functionality of  Rotation Surfaces:
  //
  template<LVSC LVSCKind, typename Derived>
  class RotationShell: public MechElement<LVSCKind>
  {
  protected:
    //=======================================================================//
    // Types and Consts:                                                     //
    //=======================================================================//
    using ME = MechElement<LVSCKind>;

  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Orientation of the rotation axis: "alpha" (|alpha| < Pi/2)  is the angle
    // between the positive direction of the OX axis and the rotation axis
    // (from the Low (Smaller-X) to the Up (Larger-X) end points). The latter is
    // assumed to lie in the OXY plane or in the OXZ plane (or both, in which
    // case it coincides with OX);  MoIs are computed wrt (xL, yL, zL)  which is
    // the Low (Smaller-X) rotation axis end, because this point is invariant
    // under the change of the Propellant level. We must have yL=0 (the rotation
    // axis is in OXY) or zL=0 (in OXZ), or both (then alpha=0 as well):
    //
    bool        m_inXY;         // If false, then inXZ holds (both may be true)
    bool        m_inXZ;         //
    double      m_cosA;         // cos(alpha)
    double      m_sinA;         // sin(alpha)
    ME::PosVE   m_up;           // Upper (Smaller-X) axis end
    Len         m_h;            // Over-all body length along the rotation axis
    ME::PosVE   m_low;          // Lower (Larger-X)  axis end: MoI ORIGIN
    Len         m_yzL;          // Low[1] or Low[2]

    // Geometric Properties:
    Area        m_sideSurfArea; // W/o the Bases
    Vol         m_enclVol;      // Nominal Volume enclosed (with imag. Bases)

    // Propellant-related flds (ie it is assumed that this rotation body may
    // contain Propellant):
    //
    Density     m_rho;          // Propellant Density (0 if no propellant)
    Mass        m_propMassCap;  // Propellant Mass Capacity

    // Coeffs for translation of "instrinsic" MoI params (J0, J1, K) into XYZ
    // MoI params (Jx, Jin, Jort and ultimately to Jx, Jy, Jz):
    // For Jx:
    double      m_Jx0;
    double      m_Jx1;
    Len         m_JxK;
    Len2        m_JxSV;
    // For Jin:
    double      m_Jin0;
    double      m_Jin1;
    Len         m_JinK;
    Len2        m_JinSV;
    // For Jort (Jort0 = Jort1 = 1):
    Len         m_JortK;
    Len2        m_JortSV;

  protected:
    // Default and Copy Ctors, Assignment and Equality are auto-generated, and
    // they are Protected:
    constexpr RotationShell()                                       = default;
    constexpr RotationShell            (RotationShell const&)       = default;
    constexpr RotationShell& operator= (RotationShell const&)       = default;
    constexpr bool           operator==(RotationShell const&) const = default;
    constexpr bool           operator!=(RotationShell const&) const = default;

    //=======================================================================//
    // "Init":                                                               //
    //=======================================================================//
    // Similar to "MechElement::Init": Initialiser similar to Non-Default Ctor
    // but used in a slightly different way:
    //
    constexpr void Init
    (
      // Params for the Base Class ("MechElement"):
      Area        a_side_surf_area,
      Vol         a_encl_vol,
      Mass        a_mass,     // If 0, then auto-calculated with SurfDens=1

      // Params for "RotationShell" itself:
      double      a_alpha,    // Rotation axis orientation
      Len         a_x0,       // The Up or Low rotation axis end
      Len         a_y0,       //
      Len         a_z0,       //
      bool        a_0is_up,   // Is (x0,y0,z0) the Up or Low end?
      Len         a_h,        // Over-all body length  (along the rotation axis)

      // "Empty" MoI Coeffs (wrt the Low end of the rotation axis):
      Len4        a_je0,
      Len4        a_je1,
      Len3        a_ke,

      // Propellant Density (may be 0 if no Propellant is held in this Body):
      Density     a_rho
    )
    {
      //---------------------------------------------------------------------//
      // Initialise the "RotationShell" flds FIRST:                          //
      //---------------------------------------------------------------------//
      // Rotation axis orientation angle, and the over-all length:
      m_inXY = IsZero(a_z0);
      m_inXZ = IsZero(a_y0);
      assert(  m_inXY || m_inXZ);
      assert(!(m_inXY && m_inXZ) || a_alpha == 0.0);

      m_cosA  = Cos(a_alpha);
      m_sinA  = Sin(a_alpha);
      m_h     = a_h;
      assert(m_cosA > 0.0 && IsPos(m_h));

      // The Up (Larger-X) and Low (Smaller-X) ends of the rotation axis:
      if (a_0is_up)
      {
        m_up   [0] = a_x0;
        m_up   [1] = a_y0;
        m_up   [2] = a_z0;
        Len dyz    =          m_sinA * m_h;
        m_low  [0] = a_x0   - m_cosA * m_h;
        m_low  [1] = m_inXY ? (a_y0  - dyz) : 0.0_m;
        m_low  [2] = m_inXZ ? (a_z0  - dyz) : 0.0_m;
      }
      else
      {
        m_low  [0] = a_x0;
        m_low  [1] = a_y0;
        m_low  [2] = a_z0;
        Len dyz    =          m_sinA * m_h;
        m_up   [0] = a_x0   + m_cosA * m_h;
        m_up   [1] = m_inXY ? (a_y0  + dyz) : 0.0_m;
        m_up   [2] = m_inXZ ? (a_z0  + dyz) : 0.0_m;
      }
      m_yzL = m_inXY ? m_low[1] : m_low[2];

      // Side Surace Area and Nominal Volume Enclosed:
      m_sideSurfArea = a_side_surf_area;
      assert(IsPos(m_sideSurfArea));

      m_enclVol      = a_encl_vol;
      assert(IsPos(m_enclVol));

      // Propellant Params:
      m_rho          = a_rho;
      m_propMassCap  = m_rho * a_encl_vol;

      // IMPORTANT!
      // Coeffs for (Jx, Jy, Jz) wrt (J0, J1, K, SurfOrVol), where J0, J1, K
      // are computed relative to the Low end-point of the rotation axis.
      // This is just a matter of choice for 2D Shells, but becomes important
      // for 3D Propellant Volumes: as Propellants are spent, the Low end of
      // of the Propellant Volume remains unchanged, whereas the Up one moves
      // with the decreasing level of Propellant:
      //
      // Jx =  m_sinA^2   * J0 + (1.0 + m_cosA^2)    * J1 +
      //       yzL * (yzL * SurfOrVol + 2.0 * m_sinA * K) :
      m_Jx0    = Sqr(m_sinA);
      m_Jx1    = 1.0 + Sqr(m_cosA);
      m_JxK    = 2.0 * m_sinA * m_yzL;
      m_JxSV   = Sqr(m_yzL);

      // Jin = m_cosA^2   * J0 + (1.0 + m_sinA^2)    * J1 +
      //       xL  * (xL  * SurfOrVol + 2.0 * m_cosA * K) :
      m_Jin0   = Sqr(m_cosA);
      m_Jin1   = 1.0 + Sqr(m_sinA);
      m_JinK   = 2.0 * m_cosA * m_low[0];
      m_JinSV  = Sqr(m_low[0]);

      // Jort =  J0  + J1    + (Sqr(xL)  + Sqr(yzL)) * SurfOrVol +
      //         2.0 * (cosA * xL + sinA * yzL) * K :
      m_JortK  = 2.0 * (m_cosA  * m_low[0] + m_sinA * m_yzL);
      m_JortSV = Sqr(m_low[0])  + Sqr(m_yzL);

      //---------------------------------------------------------------------//
      // Checks:                                                             //
      //---------------------------------------------------------------------//
      // The rotation axis lies in the OXY plane, or OXZ plane, or in both:
      assert(m_inXY || m_inXZ);

      // If the rotation axis lies in both OXY and OXZ, then it must coincide
      // with OX, and therefore, we must have alpha=0:
      assert(!(m_inXY && m_inXZ) || (m_cosA == 1.0 && m_sinA == 0.0));

      // In general, |alpha| < Pi/2:
      assert(m_cosA > 0.0);

      // Propellant Density must be non-negative (0 is OK if no propellant):
      assert(!IsNeg(m_rho));

      //---------------------------------------------------------------------//
      // Now the Base Class ("MechElement"):                               //
      //---------------------------------------------------------------------//
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
      // Check: "m_inXY", "m_inXZ" derived from "a_{xyz}0" must be consistent
      // with "emptyCoM":
      assert(!m_inXY || IsZero(emptyCoM[2]));
      assert(!m_inXZ || IsZero(emptyCoM[1]));

      // CoMDots and MoIDots must all be 0s, of course:
      assert(IsZero(emptyCoMDots[0]) && IsZero(emptyCoMDots[1]) &&
             IsZero(emptyCoMDots[2]) && IsZero(emptyMoIDots[0]) &&
             IsZero(emptyMoIDots[1]) && IsZero(emptyMoIDots[2]));

      // Finally, invoke the Base Class Initialiser:
      MechElement<LVSCKind>::Init
        (this, emptyCoM, emptyCoMDots, emptyMass, MDot0, emptyMoIs,
         emptyMoIDots,   isFinal);
    }

    //=======================================================================//
    // "MoIsCoM": Computation of XYZ MoIs and CoM from "intrinsic" MoIs:     //
    //=======================================================================//
    // The same functions are applicable to both MoI per SurfDens (Len4) and
    // MoI per (Volume) Density (Len5), hence:
    // J=Len4, K=Len3, SV=Len2, OR
    // J=Len5, K=Len4, SV=Len3:
    //
    template<typename J>
    constexpr void MoIsCoM
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
      J Jx   = m_Jx0  * a_j0 + m_Jx1  * a_j1 + m_JxK   * a_k + m_JxSV   * a_sv;
      J Jin  = m_Jin0 * a_j0 + m_Jin1 * a_j1 + m_JinK  * a_k + m_JinSV  * a_sv;
      J Jort =          a_j0 +          a_j1 + m_JortK * a_k + m_JortSV * a_sv;
      J Jy   = m_inXY ? Jin  : Jort;
      J Jz   = m_inXZ ? Jin  : Jort;

      // MoI Rates: Linear transforms similar to above:
      using JRate   = decltype(a_j0 / 1.0_sec);
      JRate JxDot   =
        m_Jx0   * a_j0_dot   + m_Jx1  * a_j1_dot + m_JxK    * a_k_dot +
        m_JxSV  * a_sv_dot;
      JRate JinDot  =
        m_Jin0  * a_j0_dot   + m_Jin1 * a_j1_dot + m_JinK   * a_k_dot +
        m_JinSV * a_sv_dot;
      JRate JortDot =
        a_j0_dot + a_j1_dot  + m_JortK * a_k_dot + m_JortSV * a_sv_dot;

      JRate JyDot   = m_inXY ? JinDot : JortDot;
      JRate JzDot   = m_inXZ ? JinDot : JortDot;

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

      // Now the CoM: "xiC" is its co-ord along the rotation axis, relative to
      // the Low (Smallest-X) axis end, hence positive:
      Len  xiC      = a_k / a_sv;
      assert(IsPos(xiC));
      a_com[0]      = m_low[0]     + m_cosA * xiC;
      Len  yzC      = m_yzL        + m_sinA * xiC;
      a_com[1]      = m_inXY ? yzC : 0.0_m;
      a_com[2]      = m_inXZ ? yzC : 0.0_m;

      // Similar for CoMDots:
      Vel xiCDot    = a_k_dot / a_sv - a_k * a_sv_dot / Sqr(a_sv);
      a_com_dots[0] = m_cosA * xiCDot;
      Vel yzCDot    = m_sinA * xiCDot;
      a_com_dots[1] = m_inXY ? yzCDot : Vel(0.0);
      a_com_dots[2] = m_inXZ ? yzCDot : Vel(0.0);
    }

  public:
    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    constexpr Area     GetSideSurfArea()  const { return m_sideSurfArea; }
    constexpr Vol      GetEnclVol()       const { return m_enclVol;      }
    constexpr Len      GetHeight ()       const { return m_h;            }

    // Upper and Lower Points:
    constexpr typename ME::PosVE const& GetUp () const { return m_up;    }
    constexpr typename ME::PosVE const& GetLow() const { return m_low;   }

    // The Maximum Propellant Mass (Capacity):
    constexpr Mass     GetPropMassCap()   const { return m_propMassCap;  }

    // The Propellant Density:
    constexpr Density  GetPropDens()      const { return m_rho;          }

    //=======================================================================//
    // "GetPropBulkME":                                                      //
    //=======================================================================//
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
    constexpr MechElement<LVSCKind> GetPropBulkME
    (
      Mass     a_prop_mass,
      MassRate a_prop_mass_dot = MassRate(0.0),    // Must be <= 0 in general
      Len*     a_prop_level    = nullptr
    )
    const
    {
      // Check the limits of the Propellant Mass. For the 2nd inequaluty, allow
      // some floating-point tolerance:
      assert(!IsNeg(a_prop_mass) && a_prop_mass <= GetPropMassCap() * TolFact);

      // The Propellant Volume: Make sure it is within the limits to avoid
      // rounding errors:
      Vol     propVol    = a_prop_mass     / GetPropDens();
      VolRate propVolDot = a_prop_mass_dot / GetPropDens();
      assert (!IsNeg(propVol) && propVol <= GetEnclVol() * TolFact &&
              !IsPos(propVolDot));
      propVol = std::min(propVol, GetEnclVol()); // Enforce it for safety

      // The Propellant Level, relative to the Low (Smallest-X) Base.
      // XXX: We assume that the propellant surface is always orthogonal to the
      // rotation axis, due to the tank pressurisation:
      auto [l, lDot] =
        static_cast<Derived const*>(this)->Derived::PropLevelOfVol
          (propVol, propVolDot);

      // Mathematically, we must have 0 <= l <= h, where l=0 corresponds to the
      // empty Element, and l=h to the full one. We enforce this interval expli-
      // citly to prevent rounding errors:
      assert(!IsNeg(l) && l <= m_h * TolFact);
      l        = std::min(l,   m_h);

      if (a_prop_level != nullptr)
        * a_prop_level  = l;

      // Get the MoIs and the CoM of the Popellant. NB: Needless to say, use the
      // current "propVol" here, NOT the maximum "GetEnclVol()". NB: "Intrinsic"
      // MoI components of the Propellant are computed by "Derived":
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
      // (*) the result is Final:
      //
      return MechElement<LVSCKind>
        (com, comDots, a_prop_mass, a_prop_mass_dot, mois, moiDots, true);
    }
  };
}
// End namespace SpaceBallistics

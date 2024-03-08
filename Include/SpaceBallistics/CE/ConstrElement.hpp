// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CE/ConstrElement.hpp":                //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include <boost/container/static_vector.hpp>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <initializer_list>

namespace SpaceBallistics
{
  //=========================================================================//
  // "ConstrElement": Base Class for Construction Elements:                  //
  //=========================================================================//
  // General characteristics of standard construction elements:
  // Currently, 2D Shells (Truncated Cones and Spherical Segments), 3D Volumes
  // (of the same shape as Shells) and Point Masses are supported (via the cor-
  // resp Derived Classes of "ConstrElement").  For each them,  we can compute
  // Centers of Masses (CoM) and Moments of Inertia (MoI) wrt OX, OY, OZ axes.
  // NB:
  // It is assumed that the co-ords system OXYZ is such that  OX  is  the LV's
  // principal axis of symmetry, with the positive direction pointing UPWARDS
  // (ie from the Tail to the Nose). This is important when it comes to the CoM
  // and MoIs of the contained Propellant:    we assume that the Propellant is
  // concentrated at the bottom (Lower, Smaller-X) part of the "ConstrElement":
  //
  class ConstrElement
  {
  public:
    //=======================================================================//
    // Consts and Types:                                                     //
    //=======================================================================//
    // "UnKnownMass": Param to be used for "ConstElement"s when their mass is
    // not yet known (all real Masses are of course strictly positive):
    //
    constexpr static Mass UnKnownMass = 0.0_kg;

    // The following types might be useful for derived classes:
    using Len2 = decltype(Sqr     (1.0_m));    // Same as "Area"
    using Len3 = decltype(Cube    (1.0_m));    // Same as "Vol"
    using Len4 = decltype(Sqr(Sqr (1.0_m)));
    using Len5 = decltype(Sqr(Sqr (1.0_m)) * 1.0_m);
    using Len6 = decltype(Sqr(Cube(1.0_m)));

  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Over-All Dynamical Params of this "ConstrElement":
    //
    Len         m_CoM [3];  // (X,Y,Z) co-ords of the Center of Masses
    Mass        m_mass;     // Mass
    MoI         m_MoIs[3];  // Moments of Inertia wrt the OX, OY and OZ axes
    bool        m_isFinal;  // If not set, "m_mass", "m_MoIs" are not valid yet
                            //   (but the CoM still is!)
  protected:
    //=======================================================================//
    // "Init":                                                               //
    //=======================================================================//
    // For use by Derived Classes. It is similar to a Non-Default Ctor, but is
    // intended to be invoked AFTER the Derived flds have been initialised, so
    // implemented as a Static Method, not a Ctor.
    // XXX: It is Static becaise without that, it could not be called from me-
    // thods of derived classes on "non-this" objs:
    //
    constexpr static void Init
    (
      ConstrElement*  a_ce,
      Len  const      a_com [3],
      Mass            a_mass,
      MoI  const      a_mois[3],
      bool            a_is_final
    )
    {
      assert(a_ce != nullptr);
      a_ce->m_CoM [0]  = a_com [0];
      a_ce->m_CoM [1]  = a_com [1];
      a_ce->m_CoM [2]  = a_com [2];
      a_ce->m_mass     = a_mass;
      a_ce->m_MoIs[0]  = a_mois[0];
      a_ce->m_MoIs[1]  = a_mois[1];
      a_ce->m_MoIs[2]  = a_mois[2];
      a_ce->m_isFinal  = a_is_final;

      // NB: Even if "a_mass" is not final, it must be positive:
      assert( IsPos(a_ce->m_mass)    && !IsNeg(a_ce->m_MoIs[0]) &&
             !IsNeg(a_ce->m_MoIs[1]) && !IsNeg(a_ce->m_MoIs[2]));
    }

  public:
    //=======================================================================//
    // Default Ctor:                                                         //
    //=======================================================================//
    // NB:
    // (*) All numeric fields are initialised to 0s, not to NaNs, so the "empty"
    //     obj (constructed by the Default Ctor) can be used as an initial value
    //     for summation ("+", "+=", etc);
    // (*) IsFinal is set to "true", because the empty "ConstrElement" is typi-
    //     cally used as the base for "+" which requires Final components:
    //
    constexpr ConstrElement()
    : m_CoM       {0.0_m,    0.0_m,    0.0_m},
      m_mass      (0.0),
      m_MoIs      {MoI(0.0), MoI(0.0), MoI(0.0)},
      m_isFinal   (true)     // !!!
    {}

    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr ConstrElement
    (
      Len         a_com [3],
      Mass        a_mass,
      MoI         a_mois[3],
      bool        a_is_final
    )
    { Init(this, a_com, a_mass, a_mois, a_is_final); }

    // Copy Ctor, Assignment and Equality are auto-generated...

    //=======================================================================//
    // "GetMassScale":                                                       //
    //=======================================================================//
    // Given a list of "ConstrElement"s  (all of them must have Non-Final Mass-
    // es) and the real Total Mass of them,  returns  a  dimension-less  factor
    // which could be applied (in "ProRateMass") to each "ConstrELement" to set
    // its correct Final Mass. It assumes that the densities  (surface  or vol)
    // of all "ConstrElement"s in the list are the same, ie their relative mas-
    // ses remain unchanged after scaling:
    //
    constexpr static double GetMassScale
      (std::initializer_list<ConstrElement const*> a_ces, Mass a_total_mass)
    {
      assert(IsPos(a_total_mass));

      // All "a_ces" must NOT be Final yet. Calculate their nominal total mass:
      Mass nomTotal = 0.0_kg;
      for (ConstrElement const* ce: a_ces)
      {
        assert(ce != nullptr && !(ce->m_isFinal));
        nomTotal  += ce->m_mass;
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
    // Returns a new object (of any type derived from "ConstrElement") which
    // differes from the original one  by the Mass and MoIs being multiplied
    // by the ScaleFactor. XXX: It is assumed that any other flds of "Derived"
    // are UNAFFECTED by this scaling:
    //
    template<typename Derived>
    constexpr static  Derived ProRateMass(Derived const& a_der, double a_scale)
    {
      assert(a_scale > 0.0);

      // Create a copy of "a_der" using the Copy Ctor of Derived (which must
      // indeed be derived from "ConstrElement"):
      Derived copy(a_der);

      // Cannot adjust the Mass once it has been finalised:
      assert(!copy.m_isFinal);

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
    using     Point = decltype(m_CoM);
    using     MoIs  = decltype(m_MoIs);

    constexpr Point const& GetCoM() const { return m_CoM; }

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

    constexpr MoIs const& GetMoIs() const
    {
      assert(m_isFinal);
      return m_MoIs;
    }

    //=======================================================================//
    // Addition / Subtraction:                                               //
    //=======================================================================//
    // The summands must both have FINAL masses. Only one summand is allowed to
    // have zero mass, otherwise we cannot compute the CoM.  WE ASSUME that the
    // summands DO NOT INTERSECT in space:
    //
    constexpr ConstrElement& operator+= (ConstrElement const& a_right)
    {
      assert(m_isFinal && a_right.m_isFinal &&
             !(IsZero(m_mass) && IsZero(a_right.m_mass)));

      // Masses, SurfAreas, Vols and MoIs are directly-additive:
      Mass    m0  = m_mass;
      m_mass     += a_right.m_mass;
      m_MoIs[0]  += a_right.m_MoIs[0];
      m_MoIs[1]  += a_right.m_MoIs[1];
      m_MoIs[2]  += a_right.m_MoIs[2];

      // For the CoM, do the weighted avg (but the total mass must be non-0):
      assert(IsPos(m_mass));
      double mu0  = double(m0             / m_mass);
      double mu1  = double(a_right.m_mass / m_mass);
      m_CoM[0]    = mu0 * m_CoM[0]  + mu1 * a_right.m_CoM[0];
      m_CoM[1]    = mu0 * m_CoM[1]  + mu1 * a_right.m_CoM[1];
      m_CoM[2]    = mu0 * m_CoM[2]  + mu1 * a_right.m_CoM[2];
      return *this;
    }

    constexpr ConstrElement operator+ (ConstrElement const& a_right) const
    {
      ConstrElement res = *this;
      res += a_right;
      return res;
    }

    // Subtraction is also possible, but ONLY IF the resulting SurfArea, Vol
    // and Mass remain valid. USE THIS OP WITH CARE: It is intended for removal
    // of particular components from "ConstrElement" systems:
    //
    constexpr ConstrElement& operator-= (ConstrElement const& a_right)
    {
      assert(m_isFinal && a_right.m_isFinal &&
             !(IsZero(m_mass) && IsZero(a_right.m_mass)));

      // Masses, SurfAreas, Vols and MoIs are directly-additive/subtractable:
      Mass m0     = m_mass;
      m_mass     -= a_right.m_mass;
      m_MoIs[0]  -= a_right.m_MoIs[0];
      m_MoIs[1]  -= a_right.m_MoIs[1];
      m_MoIs[2]  -= a_right.m_MoIs[2];
      // CHECKS:
      assert(IsPos(m_mass)    && IsPos(m_MoIs[0]) &&
             IsPos(m_MoIs[1]) && IsPos(m_MoIs[2]));

      // For the CoM, do the weighted avg:
      double mu0  = double(m0             / m_mass);
      double mu1  = double(a_right.m_mass / m_mass);
      m_CoM[0]    = mu0 * m_CoM[0]  - mu1 * a_right.m_CoM[0];
      m_CoM[1]    = mu0 * m_CoM[1]  - mu1 * a_right.m_CoM[1];
      m_CoM[2]    = mu0 * m_CoM[2]  - mu1 * a_right.m_CoM[2];
      return *this;
    }

    constexpr ConstrElement operator- (ConstrElement const& a_right) const
    {
      ConstrElement res = *this;
      res -= a_right;
      return res;
    }
  };

  //=========================================================================//
  // "PointMass" Class:                                                      //
  //=========================================================================//
  // A positive mass concentrated in the (x0, y0, z0) point:
  //
  class PointMass final: public ConstrElement
  {
  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    constexpr PointMass(Len a_x0, Len a_y0, Len a_z0, Mass a_mass)
    : ConstrElement()  // Slightly sub-optimal...
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
      // dered to be Final:
      ConstrElement::Init(this, pt, a_mass, mois, true);
    }
  };

  //=========================================================================//
  // "RotationBody" Class:                                                   //
  //=========================================================================//
  // SubClass of "ConstrElement", SuperClass of "TrCone", "SpherSegment" and
  // others. Provides common functionality of  Rotation Surfaces and Bodies:
  //
  class RotationBody: public ConstrElement
  {
  protected:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // Computation Tolerances:
    constexpr static double Tol     = 100.0 * Eps<double>;
    constexpr static double TolFact = 1.0   + Tol;

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
    Len         m_up [3];       // Upper (Smaller-X) axis end
    Len         m_h;            // Over-all body length along the rotation axis
    Len         m_low[3];       // Lower (Larger-X)  axis end: MoI ORIGIN
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

    // Coeffs for computation of "intrinsic" MoI params (J0, J1, K) as polynom-
    // ial functions of the propellant level.  XXX: Interestingly, the degrees
    // of such polynomials are the same for the concrete hapes currently imple-
    // mented ("TrCone", "SpherSegm", ...), though extra coeffs may potentially
    // be required for other shapes:
    double      m_JP05;     // coeff(JP0, l^5)
    Len         m_JP04;     // coeff(JP0, l^4)
    Len2        m_JP03;     // coeff(JP0, l^3)
    double      m_JP15;     // coeff(JP1, l^5)
    Len         m_JP14;     // coeff(JP1, l^4)
    Len2        m_JP13;     // coeff(JP1, l^3)
    Len3        m_JP12;     // coeff(JP1, l^2)
    Len4        m_JP11;     // coeff(JP1, l)
    double      m_KP4;      // coeff(KP,  l^4)
    Len         m_KP3;      // coeff(KP,  l^3)
    Len2        m_KP2;      // coeff(KP,  l^2)

    // Emulating virtual methods:
    // The following function ptr is provided by derived classes:
    typedef    Len (*LevelOfVol) (Vol, RotationBody const*);
    LevelOfVol m_LoV;

  protected:
    // Default and Copy Ctors, Assignment and Equality are auto-generated, and
    // they are Protected:
    constexpr RotationBody()                                      = default;
    constexpr RotationBody            (RotationBody const&)       = default;
    constexpr RotationBody& operator= (RotationBody const&)       = default;
    constexpr bool          operator==(RotationBody const&) const = default;
    constexpr bool          operator!=(RotationBody const&) const = default;

    //=======================================================================//
    // "Init":                                                               //
    //=======================================================================//
    // Similar to "ConstrElement::Init": Initialiser similar to Non-Default Ctor
    // but used in a slightly different way:
    //
    constexpr void Init
    (
      // Params for the Base Class ("ConstrElement"):
      Area        a_side_surf_area,
      Vol         a_encl_vol,
      Mass        a_mass,     // If 0, then auto-calculated with SurfDens=1

      // Params for "RotationBody" itself:
      double      a_alpha,    // Rotation axis orientation
      Len         a_x0,       // The Up or Low rotation axis end
      Len         a_y0,       //
      Len         a_z0,       //
      bool        a_0is_up,   // Is (x0,y0,z0) the Up or Low end?
      Len         a_h,        // Over-all body length  (along the rotation axis)

      // "Empty" MoI Coeffs (wrt the LOW end of the rotation axis):
      Len4        a_je0,
      Len4        a_je1,
      Len3        a_ke,

      // Propellant Density (may be 0):
      Density     a_rho,

      // Propellant Vol -> Propellant Level:
      LevelOfVol  a_lov,

      // Propellant MoI Coeffs as functions of Propellant Level (also wrt the
      // LOW end of the rotation axis):
      double      a_jp05,
      Len         a_jp04,
      Len2        a_jp03,

      double      a_jp15,
      Len         a_jp14,
      Len2        a_jp13,
      Len3        a_jp12,
      Len4        a_jp11,

      double      a_kp4,
      Len         a_kp3,
      Len2        a_kp2
    )
    {
      //---------------------------------------------------------------------//
      // Initialise the "RotationBody" flds FIRST:                           //
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
      m_LoV          = a_lov;

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

      // Coeffs for computing (J0, J1, K) as polynomial functions of the Propel-
      // lant Level:
      m_JP05 = a_jp05;
      m_JP04 = a_jp04;
      m_JP03 = a_jp03;

      m_JP15 = a_jp15;
      m_JP14 = a_jp14;
      m_JP13 = a_jp13;
      m_JP12 = a_jp12;
      m_JP11 = a_jp11;

      m_KP4  = a_kp4;
      m_KP3  = a_kp3;
      m_KP2  = a_kp2;

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

      // The LevelOfVolume function ptr must be non-NULL:
      assert(m_LoV != nullptr);

      //---------------------------------------------------------------------//
      // Now the Base Class ("ConstrElement"):                               //
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
      Len emptyCoM [3];
      MoI emptyMoIs[3];
      MoIsCoM
        (a_je0, a_je1, a_ke, a_side_surf_area, surfDens, emptyCoM, emptyMoIs);

      // Check: "m_inXY", "m_inXZ" derived from "a_{xyz}0" must be consistent
      // with "emptyCoM":
      assert(!m_inXY || IsZero(emptyCoM[2]));
      assert(!m_inXZ || IsZero(emptyCoM[1]));

      // Finally, invoke the Base Class Initialiser:
      ConstrElement::Init(this, emptyCoM, emptyMass, emptyMoIs, isFinal);
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
      decltype(J(1.0)/1.0_m)     a_k,       // Len3 or Len4
      decltype(J(1.0)/Len2(1.0)) a_sv,      // Len2 or Len3 (ie SurfArea or Vol)
      decltype(1.0_kg/a_sv)      a_dens,    // SurfDens or Density
      Len                        a_com [3], // Output
      MoI                        a_mois[3]  // ditto
    )
    const
    {
      // "Intrinsic" MoI components must be > 0 for any finite body size,  as
      // well as "a_k" (because it's a moment relative to the Low (Smallest-X)
      // point):
      assert(IsPos(a_j0) && IsPos(a_j1) && IsPos(a_sv) && IsPos(a_k) &&
             IsPos(a_sv));
      // MoIs:
      J Jx   = m_Jx0  * a_j0 + m_Jx1  * a_j1 + m_JxK   * a_k + m_JxSV   * a_sv;
      J Jin  = m_Jin0 * a_j0 + m_Jin1 * a_j1 + m_JinK  * a_k + m_JinSV  * a_sv;
      J Jort =          a_j0 +          a_j1 + m_JortK * a_k + m_JortSV * a_sv;
      J Jy   = m_inXY ? Jin  : Jort;
      J Jz   = m_inXZ ? Jin  : Jort;

      a_mois[0]  = a_dens * Jx;
      a_mois[1]  = a_dens * Jy;
      a_mois[2]  = a_dens * Jz;
      assert(!(IsNeg(a_mois[0]) || IsNeg(a_mois[1]) || IsNeg(a_mois[2])));

      // Now the CoM: "xiC" is its co-ord along the rotation axis, relative to
      // the Low (Smallest-X) axis end, hence positive:
      Len  xiC   = a_k / a_sv;
      assert(IsPos(xiC));
      a_com[0]   = m_low[0]     + m_cosA * xiC;
      Len  yzC   = m_yzL        + m_sinA * xiC;
      a_com[1]   = m_inXY ? yzC : 0.0_m;
      a_com[2]   = m_inXZ ? yzC : 0.0_m;
    }

  public:
    //=======================================================================//
    // Elementary Accessors:                                                 //
    //=======================================================================//
    constexpr Area GetSideSurfArea() const { return m_sideSurfArea; }
    constexpr Vol  GetEnclVol()      const { return m_enclVol;      }
    constexpr Len  GetHeight ()      const { return m_h;            }

    //=======================================================================//
    // "GetPropCE":                                                          //
    //=======================================================================//
    // Constructs a "ConstrElement" which holds the CoM and MoI of the Propel-
    // lant filling this RotationBody, with the current propellant mass given
    // by "a_prop_mass".   The resulting "ConstrElement" does NOT include the
    // Shell!
    // Uses the "LevelOfVol" func installed by the derived class, in the NON-
    // virtual way.
    // Returns a ficticious "ConstrElement" obj acting as a container for the
    // results computed, and suitable for the "+" operation. May also return
    // the curr Propellant Level via the 2nd arg (if non-NULL), primarily for
    // debugging purposes.
    // Although declared "constexpr", this function will mostly be called at
    // Run-Time. For the info, it also returns the PropLevel if the output ptr
    // is non-NULL:
    //
    constexpr ConstrElement   GetPropCE
      (Mass a_prop_mass, Len* a_prop_level = nullptr) const
    {
      // Check the limits of the Propellant Mass. For the 2nd inequaluty, allow
      // some floating-point tolerance:
      assert(!IsNeg(a_prop_mass) && a_prop_mass <= GetPropMassCap() * TolFact);

      // The Propellant Volume: Make sure it is within the limits to avoid
      // rounding errors:
      Vol propVol = a_prop_mass / GetPropDens();
      assert    (!IsNeg (propVol) && propVol <= GetEnclVol() * TolFact);
      propVol = std::min(propVol, GetEnclVol());

      // The Propellant Level, relative to the Low (Smallest-X) Base.
      // XXX: We assume that the propellant surface is always orthogonal to the
      // rotation axis, due to the tank pressurisation:
      Len l = m_LoV(propVol, this);

      // Mathematically, we must have 0 <= l <= h, where l=0 corresponds to the
      // empty Element, and l=h to the full one. We enforce this interval expli-
      // citly to prevent rounding errors:
      assert(!IsNeg(l) && l <= m_h * TolFact);
      l        = std::min(l,   m_h);

      if (a_prop_level != nullptr)
        * a_prop_level  = l;

      // "Intrinsic" MoIs components of the Priopellant:
      Len2 l2  = Sqr (l);
      Len3 l3  = l2 * l;
      Len5 JP0 = (   (m_JP05 * l + m_JP04) * l + m_JP03) * l3;
      Len5 JP1 = (((((m_JP15 * l + m_JP14) * l + m_JP13) * l) + m_JP12)
                             * l + m_JP11) * l;
      Len4 KP  = (   (m_KP4  * l + m_KP3 ) * l + m_KP2 ) * l2;
      assert(!(IsNeg(JP0) || IsNeg(JP1) || IsNeg(KP)));

      // Get the MoIs and the CoM of the Popellant (returned straight to the
      // CallER). NB: Needless to say, use the current "propVol" here, NOT
      // the maximum "GetEnclVol()":
      Len com [3];
      MoI mois[3];
      MoIsCoM(JP0, JP1, KP, propVol, GetPropDens(), com, mois);

      // Construct the result (XXX: this is slightly sub-optimal, as involves
      // copying of "com" and "mois"). NB:
      // (*) We must  install the correct Mass of the Propellant in the "res",
      //     so that  it can then participate in "operator+";
      // (*) SurfArea is not computed and remains 0;, we do not compute it;
      //     this is OK because by SurfArea we mean that of the Shell (with
      //     some non-0 SurfDens), whereas in this case, the Shell is not
      //     included in the result;
      // (*) the result is Final:
      //
      return ConstrElement(com, a_prop_mass, mois, true);
    }

    //=======================================================================//
    // Accessors (for use by Derived Classes):                               //
    //=======================================================================//
    constexpr Point const& GetLow()          const { return m_low; }
    constexpr Point const& GetUp ()          const { return m_up;  }

    // The Maximum Propellant Mass (Capacity):
    constexpr Mass         GetPropMassCap()  const { return m_propMassCap; }

    // The Propellant Density:
    constexpr Density      GetPropDens()     const { return m_rho; }
  };
}
// End namespace SpaceBallistics

// vim:ts=2:et
//===========================================================================//
//                      "SpaceBallistics/ME/MechElement.h":                  //
//                   Geometric Objects as "Mechanical Elements"              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/LVSC/LVSC.h"
#include "SpaceBallistics/CoOrds/EmbeddedCOS.h"
#include <boost/container/static_vector.hpp>
#include <cassert>
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
  // tanks pressurisation:
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
    using MoIVE     = MoIV    <ECOS>;

    // The Time Rate of the MoI Tensor:
    using MoIRateVE = MoIRateV<ECOS>;

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
    Vol         m_enclVol;     // Enclosed Volume
    MoIVE       m_MoIs;        // Moments of Inertia wrt the OX, OY and OZ axes
    MoIRateVE   m_MoIDots;     // d(MoIs)/dt
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
      Vol             a_encl_vol,
      MoI     const   a_mois    [3],
      MoIRate const   a_moi_dots[3],
      bool            a_is_final
    );

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
    constexpr MechElement();

    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr MechElement
    (
      Len      const a_com     [3],
      Vel      const a_com_dots[3],
      Mass           a_mass,
      MassRate       a_mass_dot,
      Vol            a_encl_vol,
      MoI      const a_mois    [3],
      MoIRate  const a_moi_dots[3],
      bool           a_is_final
    );

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
      (std::initializer_list<MechElement const*> a_mes, Mass a_total_mass);

    //=======================================================================//
    // "ProRateMass":                                                        //
    //=======================================================================//
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
    template<bool IsPressnGas = false, typename Derived>
    constexpr static  Derived ProRateMass
    (
      Derived const& a_der,
      double         a_scale,
      DensRate       a_dens_dot = DensRate(0.0)
    );

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

    constexpr Vol              GetEnclVol() const
      { return m_enclVol; }    // For both Final and Non-Final MEs

    constexpr MoIVE     const& GetMoIs()    const
    {
      assert(m_isFinal);
      return m_MoIs;
    }

    constexpr MoIRateVE const& GetMoIDots() const
    {
      assert(m_isFinal);
      return m_MoIDots;
    }

    constexpr VelVE     const& GetCoMDots() const
    {
      assert(m_isFinal);
      return m_CoMDots;
    }

    //=======================================================================//
    // Addition / Subtraction:                                               //
    //=======================================================================//
    // The operands must both have FINAL masses. Only one summand is allowed to
    // have zero mass, otherwise we cannot compute the CoM.  WE ASSUME that the
    // summands DO NOT INTERSECT in space, and for subtraction, the 2nd operand
    // entirely belongs to the 1st one:
    //
    constexpr MechElement& operator+= (MechElement const& a_right);

    constexpr MechElement  operator+  (MechElement const& a_right) const;

    constexpr MechElement& operator-= (MechElement const& a_right);

    constexpr MechElement  operator-  (MechElement const& a_right) const;
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
    constexpr PointMass(Len a_x0, Len a_y0, Len a_z0, Mass a_mass);
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
    // Orientation of the rotation axis:
    // (*) "alpha" (0 <= alpha < Pi/2)  is the angle between the positive direc-
    //     tion of the OX axis and the positive direction rotation axis (Xi, ie
    //     from the Low (Smaller-X) to the Up (Larger-X) end point); the  angle
    //     is in the counter-clock-wise direction as usual;
    // (*) "psi" (in 0..2*Pi or -Pi..+Pi) is the polar angle in the OYZ plane of
    //     the intersection line between the OYZ plane and the plane made by the
    //     OX axis and the rotation axis (Xi), viewed from the positive directi-
    //     on of the OX axis.
    // MoIs are computed wrt (xL, yL, zL)  which is the Low (Smaller-X) rotation
    // axis end, because for Propellant Bulks, this point is INVARIANT under the
    // change of the Propellant level.
    //
    ME::PosVE   m_up;           // Up  (Larger-X)  axis end
    Len         m_h;            // Over-all body length along the rotation axis
    ME::PosVE   m_low;          // Low (Smaller-X) axis end:  MoI ORIGIN
    double      m_xi[3];        // Directing Cosines of the Xi axis, Low->Up

    // Geometric Properties:
    Area        m_sideSurfArea; // W/o the Bases

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
    // For Jy:
    double      m_Jy0;
    double      m_Jy1;
    Len         m_JyK;
    Len2        m_JySV;
    // For Jz:
    double      m_Jz0;
    double      m_Jz1;
    Len         m_JzK;
    Len2        m_JzSV;

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
    );

    //=======================================================================//
    // "MoIsCoM": Computation of XYZ MoIs and CoM from "intrinsic" MoIs:     //
    //=======================================================================//
    // The same function is applicable to both MoI per SurfDens (Len4) and MoI
    // per (Volume) Density (Len5), hence:
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
    const;

  public:
    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    constexpr Area     GetSideSurfArea()  const { return m_sideSurfArea; }
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
      Mass     a_prop_mass     = Mass(Inf<double>), // >= 0; +oo => FullLoad
      MassRate a_prop_mass_dot = MassRate(0.0),     // Must be <= 0 in gen
      Len*     a_prop_level    = nullptr            // Not computed if NULL
    )
    const;
  };
}
// End namespace SpaceBallistics

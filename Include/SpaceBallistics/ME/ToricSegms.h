// vim:ts=2:et
//===========================================================================//
//                      "SpaceBallistics/ME/ToricSegms.h":                   //
//                  Geometric Objects as "Mechanical Elements"               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/MechElement.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "ToricSegm" Class:                                                      //
  //=========================================================================//
  // Describes the Upper or Lower part of a Toric Shell; "r" and "h"  are the
  // Radius and Height of the Segment (similar to those of SpherSegm, but un-
  // like the "SpherSegm", the "ToricSegm" is obtained by rotation of a circu-
  // lar arc around an axis which is parallel to the arc axis, but not coinci-
  // ding with the latter). "R" is the Minor Radius of the Torus (derived from
  // "r" and "h" similar to "SpherSegm"),  and "Q" is the Major Radius of the
  // Torus, that is, the distance between the aforementioned arc axis and the
  // rotation axis; we may assume r <= R <= Q.  Up-/Low-Facing orientation is
  // similar to that of a "SpherSegm":
  //
  template<LVSC LVSCKind>
  class ToricSegm final: public RotationShell<LVSCKind, ToricSegm<LVSCKind>>
  {
  private:
    using ME = MechElement  <LVSCKind>;
    using RS = RotationShell<LVSCKind, ToricSegm<LVSCKind>>;

    //=======================================================================//
    // Data Flds: Toric Segment's Geometry:                                  //
    //=======================================================================//
    Len     m_Q;        // Major Radius of the Toric Segm (NOT Torus!)
    Len     m_R;        // Minor Radius of the Torus      (not just Segm!)
    bool    m_facingUp; // Orientation
    double  m_q2;       // Memoised Coeffs
    double  m_q2L;      //
    Len     m_twoH;     //
    Len2    m_h2;       //
    Vol     m_NV;       //
    Len4    m_pQR3;     //
    Len5    m_pQR4;     //
    Len5    m_JPL0;     // MoI components for Full-Propellant in Low-Facing Sgm
    Len5    m_JPL1;     //
    Len4    m_KPL;      //

    // Default Ctor is deleted:
    ToricSegm() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr ToricSegm
    (
      bool    a_facing_up,
      Len     a_xb,          // Over-All   Segment Base Center
      Len     a_yb,          //
      Len     a_zb,          //
      double  a_cosA,        // See "RotationShell"...
      double  a_sinA,        //
      double  a_cosP,        //
      double  a_sinP,        //
      Len     a_d,           // Cross-Section Base Diameter
      Len     a_h,           // Cross-Section Height
      Len     a_D,           // Over-All   Segment Diameter
      Density a_rho,         // Propellant Density (may be 0 if no Propelt)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Non-Default Ctors, Simple Cases:                                      //
    //=======================================================================//
    // Rotation axis coinciding with OX, so alpha=0 and we can then assume psi=0
    // as well:
    //
    constexpr ToricSegm
    (
      bool    a_facing_up,
      Len     a_xb,       // Over-All Segment Base Center (X-CoOrd)
      Len     a_d,        // Cross-Section Base Diameter
      Len     a_h,        // Cross-Section Height
      Len     a_D,        // Over-All Segment Diameter
      Density a_rho,      // 0 may be OK (if holds no Propellant)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    // As above, but with d/2 = h, ie the Cross-Section is a HemiSpehere:
    //
    constexpr ToricSegm
    (
      bool    a_facing_up,
      Len     a_xb,       // Over-All Segment Base Center (X-CoOrd)
      Len     a_d,        // Cross-Section Base Diameter
      Len     a_D,        // Over-All Segment Diameter
      Density a_rho,      // 0 may be OK (if holds no Propellant)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr std::pair<Len,Vel> PropLevelOfVol(Vol a_v, VolRate a_v_dot) const;

  private:
    //-----------------------------------------------------------------------//
    // For the Low-Facing "ToricSegm":                                       //
    //-----------------------------------------------------------------------//
    constexpr std::pair<Len,Vel> PropLevelOfVolLow
      (Vol a_v, VolRate a_v_dot) const;

    //-----------------------------------------------------------------------//
    // For the Up-Facing "ToricSegm":                                        //
    //-----------------------------------------------------------------------//
    constexpr std::pair<Len,Vel> PropLevelOfVolUp
      (Vol a_v, VolRate a_v_dot) const;

  public:
    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    constexpr void PropMoIComps
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot
    )
    const;

    //-----------------------------------------------------------------------//
    // For the Low-Facing "ToricSegm":                                       //
    //-----------------------------------------------------------------------//
    // XXX: This method is also called from the "ToricSegm" Ctor, when the par-
    // ent classes are not yet initialised, so it MUST NOT call any methods of
    // "RotationShell":
    //
    constexpr void PropMoICompsLow
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot,
      Vol*      a_vp,       // May be NULL
      VolRate*  a_vp_dot    // Ditto
    )
    const;

    //-----------------------------------------------------------------------//
    // For the Up-Facing "ToricSegm":                                        //
    //-----------------------------------------------------------------------//
    // Here there is no optional return values "VP" and "VPDot":
    //
    constexpr void PropMoICompsUp
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot
    )
    const;
  };

  //=========================================================================//
  // "DoubleCylinder" Class:                                                 //
  //=========================================================================//
  // Provides a "cylindrical torus" (a body obtained by rotating a rectangle
  // around an outside axis parallel to the cylinder's main axis):
  //
  template<LVSC LVSCKind>
  class DoubleCylinder:
    public RotationShell<LVSCKind, DoubleCylinder<LVSCKind>>
  {
  private:
    using ME = MechElement  <LVSCKind>;
    using RS = RotationShell<LVSCKind, DoubleCylinder<LVSCKind>>;

    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Memoised Coeffs for Propellant MoI Components:
    Area    m_S;       // Base Area
    Len2    m_cJP0;
    Len2    m_dJP0;
    Len4    m_cJP1;
    Len2    m_cKP;
    Len2    m_dKP;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // The General Case:                                                     //
    //-----------------------------------------------------------------------//
    constexpr DoubleCylinder
    (
      Len     a_xu,          // Upper  Base Center
      Len     a_yu,          //
      Len     a_zu,          //
      double  a_cosA,        // As in "RotationShell"...
      double  a_sinA,        //
      double  a_cosP,        //
      double  a_sinP,        //
      Len     a_D,           // Outer Diameter
      Len     a_d,           // Inner Diameter
      Len     a_h,           // Height
      Density a_rho,         // Propellant Density (may be 0 if no Propelt)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //-----------------------------------------------------------------------//
    // Simple Case: Rotation Axis coincides with OX:                         //
    //-----------------------------------------------------------------------//
    // Then alpha=0 and we can assume psi=0 as well:
    //
    constexpr DoubleCylinder
    (
      Len     a_xu,          // Upper  Base Center
      Len     a_D,           // Outer Diameter
      Len     a_d,           // Inner Diameter
      Len     a_h,           // Height
      Density a_rho,         // Propellant Density (may be 0 if no Propelt)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr std::pair<Len,Vel> PropLevelOfVol(Vol a_v, VolRate a_v_dot) const;

    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    constexpr void PropMoIComps
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot
    )
    const;
  };
}
// End namespace SpaceBallistics

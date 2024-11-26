// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/ME/TrConeSpherSegm.h":                //
//                  Geometric Objects as "Mechanical Elements"               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/MechElement.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TrCone" Class: Truncated Cone Shell (w/o and with Propellant):         //
  //=========================================================================//
  // Object  : Truncated (in general) conical shell (side surface only).
  // Geometry: Upper Base Diameter "du" (for Larger-X) and Lower Base Diameter
  //           "dl" (for Smaller-X), height "h"; may have  du<=>dl;
  //           either "du" or "dl" may be be 0 (ie full cone), but  not both.
  // Location: (xu,yu,zu) -- co-ords of the Upper (Larger-X) Base center; XXX:
  //           either "zu" or "yu" must be 0. The cone axis is lying in the XY
  //           or XZ plane, resp., at the angle  "alpha" to the positive direc-
  //           tion of the OX axis; we assume that |alpha| < Pi/2;  obviously,
  //           the "du" base diameter corresponds to the base center (xu,yu)
  //           (if zu=0, ie XY plane) or (xu,zu) (if yu=0, ie XZ plane); the
  //           "dl" diameter corresponds to the other (Smaller-X) end of the
  //           cone axis segment;
  //           when alpha=0, the axis of the Spherical Segment coincides with
  //           OX, which is the most common case.
  // Mass:     If the "mass" param is set, its value is used  for the element
  //           mass which is in this case final. If this param is 0, the mass
  //           is set under the assumtion of SutfDensity=1, and it is NOT fi-
  //           nal yet (needs to be set later via "ProRateMass").
  // Return value: The resulting "MechElement" object:
  //
  template<LVSC LVSCKind>
  class TrCone final: public RotationShell<LVSCKind, TrCone<LVSCKind>>
  {
  private:
    using ME = MechElement  <LVSCKind>;
    using RS = RotationShell<LVSCKind, TrCone<LVSCKind>>;

    //=======================================================================//
    // Data Flds: Truncated Cone's Geometry:                                 //
    //=======================================================================//
    // Memoise pre-computed coeffs of the Cardano formula in "LevelOfVol":
    Len       m_deltaR;
    Len2      m_clD;
    Len3      m_clVol;
    Len2      m_rh;
    Len6      m_rh3;
    Area      m_baseArea; // Meaningful for a Cylinder only
    Len2      m_PiR2;

    // Coeffs for computation of "intrinsic" MoI params (J0, J1, K) as polynom-
    // ial functions (of degress 5 for JP0, JP1, and of degree 4 for KP) of the
    // Propellant Level, used in "PropMoIComps":
    double    m_JP05;     // coeff(JP0, l^5)
    Len       m_JP04;     // coeff(JP0, l^4)
    Len2      m_JP03;     // coeff(JP0, l^3)
    double    m_JP15;     // coeff(JP1, l^5)
    Len       m_JP14;     // coeff(JP1, l^4)
    Len2      m_JP13;     // coeff(JP1, l^3)
    Len3      m_JP12;     // coeff(JP1, l^2)
    Len4      m_JP11;     // coeff(JP1, l)
    double    m_KP4;      // coeff(KP,  l^4)
    Len       m_KP3;      // coeff(KP,  l^3)
    Len2      m_KP2;      // coeff(KP,  l^2)

    // And similar for JP0, JP1 and K Rates ("Dots"):
    double    m_JP0D4;
    Len       m_JP0D3;
    Len2      m_JP0D2;
    double    m_JP1D4;
    Len       m_JP1D3;
    Len2      m_JP1D2;
    Len3      m_JP1D1;
    Len4      m_JP1D0;
    double    m_KPD3;
    Len       m_KPD2;
    Len2      m_KPD1;

    // Default Ctor is deleted:
    TrCone() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    // See "RotationShell" for the meaning of "{cos|sin}{A|P}":
    //
    constexpr TrCone
    (
      Len     a_xu,        // Upper Base Center (Upper=Larger-X axis end)
      Len     a_yu,        //
      Len     a_zu,        //
      double  a_cosA,      // cos(alpha), >  0
      double  a_sinA,      // sin(alpha), >= 0
      double  a_cosP,      // cos(psi),   any sign
      double  a_sinP,      // sin(psi),   any sign
      Len     a_du,        // Upper (Larger-X)  Base Diameter
      Len     a_dl,        // Lower (Smaller-X) Base Diameter
      Len     a_h,         // Height
      Density a_rho,       // Propellant Density (0 if no Propellant)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Non-Default Ctors, Simple Cases (TrCone/Cylinder with OX axis):       //
    //=======================================================================//
    // "TrCone" with the rotation axis coinsiding with Ox: yu=zu=alpha=0, psi
    // is irrelevant and is also set to 0:
    //
    constexpr TrCone
    (
      Len     a_xu,       // Upper (Larger-X)  Base Center
      Len     a_du,       // Upper (Larger-X)  Base Diameter
      Len     a_dl,       // Lower (Smaller-X) Base Diameter
      Len     a_h,        // Height
      Density a_rho,      // 0 may be OK
      Mass    a_empty_mass = ME::UnKnownMass
    );

    // As above, but with dU=dL, ie a Cylinder:
    constexpr TrCone
    (
      Len     a_xu,       // Upper (Larger-X)  Base Center
      Len     a_d,        // Diameter of both  Bases
      Len     a_h,        // Height
      Density a_rho,      // 0 may be OK
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    // Returns (Level, LevelDot):
    //
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

  //=========================================================================//
  // "SpherSegm" Class:                                                      //
  //=========================================================================//
  // Similar to "TrCone", but for a Spherical Segment of the base diameter "d"
  // and the height (from the base plane to the pole) "h", where we assume
  // h <= d/2 (the equality corresponds to the case of a HemiSphere).
  // NB: this is a Shperical Segment,  NOT a Spherical Slice,  ie it always
  // contains a pole. The Spherical Slice would be a more general case and a
  // more close analogy of the Truncated Cone, but it would (arguably) have
  // almost no real applications.
  // Furthermore, "alpha" is the angle between the segment's axis and the po-
  // sitive direction of the X axis; we assume that |alpha| < Pi/2;
  // (xB, yB, zB) are the co-ords of the base center (NOT of the pole!), and
  // either "zB" or "yB" must be 0;  in the former case, the Segment's rotat-
  // ion axis is assumed to lie in the XY plane, in the latter -- in XZ.
  // The spherical segment may be facing "Up" (ie towards the positive direc-
  // tion of the OX axis), or otherise, as given by the "facing_up" flag.
  // When alpha=0, the axis of the Spherical Segment coincides with OX, which
  // is the most common case:
  //
  template<LVSC LVSCKind>
  class SpherSegm final: public RotationShell<LVSCKind, SpherSegm<LVSCKind>>
  {
  private:
    using ME = MechElement  <LVSCKind>;
    using RS = RotationShell<LVSCKind, SpherSegm<LVSCKind>>;

    //=======================================================================//
    // Data Flds: Spherical Segment's Geometry:                              //
    //=======================================================================//
    // Over-all sizes and Orientation:
    Len     m_R;          // Sphere  Radius
    Len3    m_R3;         // Sphere  Radius^3
    Len3    m_cR3;        // Pi*R^3 / 3 = HemiSphereVol / 2
    bool    m_facingUp;   // Segment Ortientation

    // Coeffs for computation of "intrinsic" MoI params (J0, J1, K) as polynom-
    // ial functions (of degress 5 for JP0, JP1, and of degree 4 for KP) of the
    // Propellant Level -- similar to "TrCone":
    //
    double  m_JP05;       // coeff(JP0, l^5)
    Len     m_JP04;       // coeff(JP0, l^4)
    Len2    m_JP03;       // coeff(JP0, l^3)
    double  m_JP15;       // coeff(JP1, l^5)
    Len     m_JP14;       // coeff(JP1, l^4)
    Len2    m_JP13;       // coeff(JP1, l^3)
    Len3    m_JP12;       // coeff(JP1, l^2)
    Len4    m_JP11;       // coeff(JP1, l)
    double  m_KP4;        // coeff(KP,  l^4)
    Len     m_KP3;        // coeff(KP,  l^3)
    Len2    m_KP2;        // coeff(KP,  l^2)

    // And similar for JP0, JP1 and K Rates ("Dots"):
    double  m_JP0D4;
    Len     m_JP0D3;
    Len2    m_JP0D2;
    double  m_JP1D4;
    Len     m_JP1D3;
    Len2    m_JP1D2;
    Len3    m_JP1D1;
    Len4    m_JP1D0;
    double  m_KPD3;
    Len     m_KPD2;
    Len2    m_KPD1;

    // Default Ctor is deleted:
    SpherSegm() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr SpherSegm
    (
      bool    a_facing_up,
      Len     a_xb,          // Base center
      Len     a_yb,          //
      Len     a_zb,          //
      double  a_cosA,        // As in "RotationShell"
      double  a_sinA,        //
      double  a_cosP,        //
      double  a_sinP,        //
      Len     a_d,           // Base diameter
      Len     a_h,           // Height
      Density a_rho,         // Propellant Density (may be 0)
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Non-Default Ctors, Simple Cases:                                      //
    //=======================================================================//
    // Rotation axis coinciding with OX, so alpha=0; in this case, we can assume
    // psi=0 as well:
    //
    constexpr SpherSegm
    (
      bool    a_facing_up,
      Len     a_xb,       // Base center X co-ord
      Len     a_d,        // Base diameter
      Len     a_h,        // Height
      Density a_rho,      // 0 may be OK
      Mass    a_empty_mass = ME::UnKnownMass
    );

    // As above, but with d/2 = h, ie a HemiSpehere:
    //
    constexpr SpherSegm
    (
      bool    a_facing_up,
      Len     a_xb,       // Base center X co-ord
      Len     a_d,        // Base diameter
      Density a_rho,      // 0 may be OK
      Mass    a_empty_mass = ME::UnKnownMass
    );

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr std::pair<Len,Vel> PropLevelOfVol(Vol a_v, VolRate a_v_dot) const;

  private:
    //-----------------------------------------------------------------------//
    // For the Low-Facing "SpherSegm":                                       //
    //-----------------------------------------------------------------------//
    constexpr std::pair<Len,Vel> PropLevelOfVolLow
      (Vol a_v, VolRate a_v_dot) const;

    //-----------------------------------------------------------------------//
    // For the Up-Facing "SpherSegm":                                        //
    //-----------------------------------------------------------------------//
    constexpr std::pair<Len,Vel> PropLevelOfVolUp
      (Vol a_v, VolRate a_v_dot) const;

  public:
    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    // Similar method to that of "TrCone":
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

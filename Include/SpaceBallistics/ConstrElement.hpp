// vim:ts=2:et
//===========================================================================//
//                             "ConstrElement.hpp":                          //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include <cmath>
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // "ConstrElement": Base Class for Construction Elements:                  //
  //=========================================================================//
  // General characteristics of standard construction elements:
  // Currently, 2D Shells (Truncated Cones and Spherical Segments), 3D Volumes
  // (of the same shape as Shells) and Point Masses are supported.  For them, 
  // we can compute Centers of Masses (CoM) and Moments of Inertia (MoI) wrt
  // the OX, OY and OZ axes.
  // XXX:
  // Here we do not explicitly specify the co-ord system OXYZ (typically a cer-
  // tain embedded one) in which the CoM  and the MoI (wrt to the OY axis)  are
  // given. It is assumed that the co-ord system OXYZ is uniquely determined by
  // the context:
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

  protected:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    Area  m_surfArea;    // Side Surface Area (w/o Bases)
    Vol   m_vol;         // Notional volume (with imaginary Bases added)
    Len   m_CoM [3];     // (X,Y,Z) co-ords of the Center of Masses
    Mass  m_mass;        // Mass
    MoI   m_MoIs[3];     // Moments of Inertia wrt the OX, OY and OZ axes
    bool  m_massIsFinal; // If not set, "m_mass" and "m_MoIs" are not valid yet

  public:
    //=======================================================================//
    // Default Ctor:                                                         //
    //=======================================================================//
    // NB: All data fields are initialised to 0s, not to NaNs, so the "empty"
    // obj (constructed by the Default Ctor) can be used as an initial value
    // for summation ("+", "+="):
    //
    constexpr ConstrElement()
    : m_surfArea   (0.0),
      m_vol        (0.0),
      m_CoM        {0.0_m,    0.0_m,    0.0_m},
      m_mass       (0.0),
      m_MoIs       {MoI(0.0), MoI(0.0), MoI(0.0)},
      m_massIsFinal(false)
    {}

    //=======================================================================//
    // Copy Ctor, Assignment, Equality:                                      //
    //=======================================================================//
    // They are trivial. Also, the Dtor is trivial, and is auto-generated:
    //
    constexpr ConstrElement            (ConstrElement const& a_right) = default;
    constexpr ConstrElement& operator= (ConstrElement const& a_right) = default;

    constexpr bool           operator==(ConstrElement const& a_right) const
      = default;
    constexpr bool           operator!=(ConstrElement const& a_right) const
      = default;

    //=======================================================================//
    // "SetMass":                                                            //
    //=======================================================================//
    // XXX: A "setter" for the Mass (invoked when we know the final mass, sets
    // it and pro-rates the mass-dependent params). After that, the mass  will
    // be final and cannot be changed anymore:
    //
    constexpr void SetMass(Mass a_final_mass)
    {
      // Cannot adjust the Mass once it has been finalised:
      assert(!m_massIsFinal);

      // Bot the curr and the new mass must be STRICTLY positive (if the former
      // is not, the scaling factor would be infinite; if the latter is not, it
      // just makes no sense):
      assert(IsPos(m_mass) && IsPos(a_final_mass));

      // If all is OK:
      double scaleCoeff = double(a_final_mass / m_mass);
      assert(scaleCoeff > 0.0);

      // Set the Mass and adjust the MoIs:
      m_mass        = a_final_mass;
      m_MoIs[0]    *= scaleCoeff;
      m_MoIs[1]    *= scaleCoeff;
      m_MoIs[2]    *= scaleCoeff;
      m_massIsFinal = true;
    }

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    // Geometrical params are always available, even if the Mass is not final
    // yet:
    constexpr Area GetSurfArea() const { return m_surfArea; }
    constexpr Vol  GetVol()      const { return m_vol;      }

    using     Point = decltype(m_CoM);
    using     MoIs  = decltype(m_MoIs);

    constexpr Point const& GetEmptyCoM() const { return m_CoM; }

    // But the Mass and MoIs  may or may not be available.
    // If not, assert failure will be signaled (NB: in C++ >= 14, "assert" is
    // allowed in "constexpr" functions, whereas throwing exceptions is not):
    //
    constexpr Mass GetEmptyMass() const
    {
      assert(m_massIsFinal);
      return m_mass;
    }

    constexpr MoIs const& GetEmptyMoIs() const
    {
      assert(m_massIsFinal);
      return m_MoIs;
    }

    //=======================================================================//
    // Addition:                                                             //
    //=======================================================================//
    // The summands must both have final or non-final masses. Only one summand
    // is allowed to have zero mass, otherwise we cannot compute the CoM:
    //
    constexpr ConstrElement& operator+= (ConstrElement const& a_right)
    {
      assert(m_massIsFinal == a_right.m_massIsFinal    &&
             !(IsZero(m_mass) && IsZero(a_right.m_mass)));

      // Masses, SurfAreas, Vols and MoIs are directly-additive:
      Mass    m0  = m_mass;
      m_surfArea += a_right.m_surfArea;
      m_vol      += a_right.m_vol;
      m_mass     += a_right.m_mass;
      m_MoIs[0]  += a_right.m_MoIs[0];
      m_MoIs[1]  += a_right.m_MoIs[1];
      m_MoIs[2]  += a_right.m_MoIs[2];

      // For the CoM, do the weighted avg (but the total mass must be non-0):
      assert(IsPos(m_mass));
      double mu0 = double(m0             / m_mass);
      double mu1 = double(a_right.m_mass / m_mass);
      m_CoM[0]   = mu0 * m_CoM[0]  + mu1 * a_right.m_CoM[0];
      m_CoM[1]   = mu0 * m_CoM[1]  + mu1 * a_right.m_CoM[1];
      m_CoM[2]   = mu0 * m_CoM[2]  + mu1 * a_right.m_CoM[2];
      return *this;
    }

    constexpr ConstrElement operator+ (ConstrElement const& a_right) const
    {
      ConstrElement res = *this;
      res += a_right;
      return res;
    }
  };

  //=========================================================================//
  // "PointMass":                                                            //
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
    {
      assert(IsPos(a_mass));

      // Initialising the Base Class's Flds:
      m_surfArea = Area(0.0);
      m_vol      = Vol (0.0);
      m_CoM[0]   = a_x0;
      m_CoM[1]   = a_y0;
      m_CoM[2]   = a_z0;

      m_mass     = a_mass;
      Len2 x2    = Sqr(a_x0);
      Len2 y2    = Sqr(a_y0);
      Len2 z2    = Sqr(a_z0);
      m_MoI[0]   = a_mass * (y2 + z2); // Distance^2 to the OX axis
      m_MoI[1]   = a_mass * (x2 + z2); // Distance^2 to the OY axis
      m_MoI[2]   = a_mass * (x2 + y2); // Distance^2 to the OZ axis

      // The point-mass is considered to be final:
      m_massIsFinal = true;
    }
  };

  //=========================================================================//
  // "TrCone": Truncated Cone Shell (w/o and with Propellant):               //
  //=========================================================================//
  // Object  : Truncated (in general) conical shell (side surface only), with
  //           SurfDensity assumed to be 1 (this may later be changed by cal-
  //           ling "ProRateMass" on the object constructed).
  // Geometry: Base diameters "d0" and "d1", height "h";  may have  d0<=>d1 ;
  //           either "d0" or "d1" may be be 0 (ie full cone), but  not both.
  // Location: XXX: Either "z0" or "y0" must be 0. The cone axis is lying in
  //           the XY or XZ plane, resp., at the angle  "alpha" to the posit-
  //           ive direction of the OX axis; we assume that |alpha| < Pi/2;
  //           the "d0" base diameter corresponds to the base center (x0,y0)
  //           (if z0=0, ie XY plane) or (x0,z0) (if y0=0, ie XZ plane); the
  //           "d1" diameter corresponds to  the other end of the cone axis
  //           segment (with X = x1 > x0);
  //           when alpha=0, the axis of the Spherical Segment coincides with
  //           OX, which is the most common case.
  // Mass:     If the "mass" param is set, its value is used  for the element
  //           mass which is in this case final. If this param is 0, the mass
  //           is set under the assumtion of SutfDensity=1, and it is NOT fi-
  //           nal yet (needs to be set later via "SetMass").
  // Return value: The resulting "ConstrElement" object.
  //
  class TrCone final: public ConstrElement
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds: Truncated Cone's Geometry:                                 //
    //-----------------------------------------------------------------------//
    Len       m_r;        // Left  (upper) base radius
    Len       m_R;        // Right (lower) base radius
    Len       m_h;        // Height
    double    m_cosA;     // cos(alpha)
    double    m_cosA2;    // cos(alpha)^2
    double    m_sinA;     // sin(alpha)
    double    m_sinA2;    // sin(alpha)^2
    Len       m_base1[3]; // Right (lower) base center
    bool      m_inXY;     // If false, then inXZ holds (but both may be true)
    bool      m_inXZ;     //

    // Propellant-related flds:
    Density   m_rho;      // Propellant Density
    Mass      m_propCap;  // Propellant Mass Capacity

    // We also memoise pre-computed coeffs of the Cardano formula and MoIs
    // computation in the 3D case (with Propellant):
    Len       m_deltaR;
    Len3      m_clVol;
    Len2      m_Rh;
    Len6      m_Rh3;
    // The following coeffs has different dimensions, so cannot make an array
    // of them -- they need to be stored individually:
    double    m_JP05;   // coeff(JP0, l^5)
    Len       m_JP04;   // coeff(JP0, l^4)
    Len2      m_JP03;   // coeff(JP0, l^3)
    double    m_JP15;   // coeff(JP1, l^5)
    Len       m_JP14;   // coeff(JP1, l^4)
    Len2      m_JP13;   // coeff(JP1, l^3)
    Len3      m_JP12;   // coeff(JP1, l^2)
    Len4      m_JP11;   // coeff(JP1, l)
    double    m_KP4;    // coeff(KP,  l^4)
    Len       m_KP3;    // coeff(KP,  l^3)
    Len2      m_KP2;    // coeff(KP,  l^2)

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    constexpr TrCone
    (
      Len     a_x0,       // Left  (upper) base center
      Len     a_y0,       //
      Len     a_z0,       //
      double  a_alpha,    // Dimension-less, in (-Pi/2 .. Pi/2)
      Len     a_d0,       // Left  (upper) base diameter
      Len     a_d1,       // Right (lower) base diameter
      Len     a_h,        // Height
      Density a_rho,      // Propellant Density (may be 0)
      Mass    a_empty_mass = UnKnownMass
    )
    {
      //---------------------------------------------------------------------//
      // This Class's Flds Initialisation:                                   //
      //---------------------------------------------------------------------//
      assert( IsPos(a_h)   && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
            !(IsZero(a_d0) && IsZero(a_d1)) && Abs(a_alpha) < Pi_2<double> &&
            ! IsNeg(a_rho) && !IsNeg(a_rho));

      // Is the rotation axis lying in the XY plane or XZ plane (or both, in
      // which case "alpha" must be 0):
      m_inXY     = IsZero(a_z0);
      m_inXZ     = IsZero(a_y0);
      assert(  m_inXY || m_inXZ );
      assert(!(m_inXY && m_inXZ) || a_alpha == 0.0);

      // The over-all sizes:
      m_r        = a_d0 / 2.0; // Left  (upper, smaller "x") base radius
      m_R        = a_d1 / 2.0; // Right (lower, larger  "x") base radius
      m_h        = a_h;

      // Rotation axis orientation angle:
      m_cosA     = Cos(a_alpha);
      m_cosA2    = Sqr(m_cosA);
      m_sinA     = Sin(a_alpha);
      m_sinA2    = Sqr(m_sinA);

      // Co-ords of the right (lower) base center:
      m_base1[0] = a_x0            + cosA * m_h;
      m_base1[1] = (m_inXY ? (a_y0 + sinA * m_h) : 0.0_m);
      m_base1[2] = (m_inXZ ? (a_z0 + sinA * m_h) : 0.0_m);

      //---------------------------------------------------------------------//
      // Parent Class's Flds Initialisation:                                 //
      //---------------------------------------------------------------------//
      m_deltaR   = m_R - m_r;
      Len2 h2    = Sqr(m_h);
      Len  s     = SqRt(Sqr(m_deltaR) + h2);
      double a   = double(m_deltaR /  m_h);
      double a2  = Sqr(a);
      double a3  = a2 * a;
      double a4  = Sqr(a2);
      Len2   R2  = Sqr(m_R);
      Len3   R3  = R2 * m_R;
      Len4   R4  = Sqr(R2);

      // Side Surface Area and Notional Volume:
      m_surfArea = Pi<double>     * s   * (m_R + m_r);
      m_vol      = Pi<double>/3.0 * m_h * (R2  + m_R * m_r + r2);

      // The mass may or may not be given. If it is given,  calculate the Surf-
      // ace Density; otherwise, assume the SurfDens to be 1.0.   Then set the
      // Mass and MoIs:
      m_massIsFinal = IsPos(a_empty_mass);
      SurfDens surfDens =
        m_massIsFinal   ? (a_empty_mass / m_surfArea) : SurfDens(1.0);
      m_mass            =
        m_massIsFinal   ?  a_empty_mass               : (m_surfArea * surfDens);

      // Components of the "empty" MoI per SurfDensity in the "embedded" axes:
      Len4 J0    =  Pi<double> * h2 * s * (m_r / 2.0  + m_R / 6.0);
      Len4 J1    =  Pi<double>/4.0  * s * (m_R + m_r) * (R2 + r2);
      Len3 K     = -Pi<double>/3.0  * s *  m_h * (2.0 * m_r + m_R);

      Len4 Jx    =  m_sinA2 * J0    + (1.0 + m_cosA2) * J1  +
                    m_yz1 * (m_yz1  * m_surfArea + 2.0 * m_sinA * K);
      assert(IsPos(Jx));

      Len4 Jin   =  m_cosA2 * J0  + (1.0 + m_sinA2) * J1  +
                    m_x1  * (m_x1 * m_surfArea  + 2.0 * m_cosA * K);
      assert(IsPos(Jin));

      Len4 Jort  =  J0  + J1    + (Sqr(m_x1)    + Sqr(m_yz1)) * m_surfArea +
                    2.0 * (cosA * m_x1  +  sinA * m_yz1) * K;
      assert(IsPos(Jort));

      m_MoIs[0]  =  m_surfDens * Jx;
      m_MoIs[1]  =  m_surfDens * (m_inXY ? Jin  : Jort);
      m_MoIs[2]  =  m_surfDens * (m_inXY ? Jort : Jin);

      // The "empty" CoM:
      // "xiC" is the "local" co-ord of the CoM (xi=0 at x=x1), so it is < 0:
      Len  xiC   =  K / m_surfArea;
      assert(IsNeg(xiC));
      m_CoM[0]   =  m_base1[0]           + cosA * xiC;
      m_CoM[1]   =  m_inXY ? (m_base1[1] + sinA * xiC) : 0.0_m;
      m_CoM[2]   =  m_inXZ ? (m_base1[2] + sinA * xiC) : 0.0_m;

      // Propellant:
      m_rho      =  a_rho;
      m_propCap  =  m_vol * m_rho;

      // Memoised coeffs of the Cardano formula (used in "GetPropParams"):
      m_clVol    =  (3.0 / Pi<double>) * h2 * m_deltaR;
      m_Rh       =  m_R * m_h;
      m_Rh3      =  Cube (m_Rh);

      // Memoised coeffs for MoIs of the 3D Propellant (see "GetPropParams"):
      m_JP05     =  Pi<double> / 5.0 * a2;
      m_JP04     = -Pi<double> / 2.0 * a  * m_R;
      m_JP03     =  Pi<double> / 3.0 * R2;

      m_JP15     =  Pi<double> /20.0 * a4;
      m_JP14     = -Pi<double> / 4.0 * a3 * m_R;
      m_JP13     =  Pi<double> / 2.0 * a2 * R2;
      m_JP12     = -Pi<double> / 2.0 * a  * R3;
      m_JP11     =  Pi<double> / 4.0 * R4;

      m_KP4      = -Pi<double> / 4.0 * a2;
      m_KP3      =  Pi<double> * 2.0 / 3.0 * a * m_R;
      m_KP2      = -Pi<double> / 2.0 * R2;
    }

    //-----------------------------------------------------------------------//
    // Non-Default Ctor, Simple Common Case:                                 //
    //-----------------------------------------------------------------------//
    // Assume that the rotation axis coinsides with Ox, then y0=z0=alpha=0:
    //
    constexpr TrCone
    (
      Len     a_x0,       // Left  (upper) base center
      Len     a_d0,       // Left  (upper) base diameter
      Len     a_d1,       // Right (lower) base diameter
      Len     a_h,        // Height
      Mass    a_empty_mass = UnKnownMass
    )
    : TrCone(a_x0, 0.0_m, 0.0_m, 0.0, a_d0, a_d1, a_h, a_empty_mass)
    {}

    //-----------------------------------------------------------------------//
    // "GetPropParams":                                                      //
    //-----------------------------------------------------------------------//
    // Fills in the CoM and MoI of the Propellant filling this TrCone, with the
    // current volume "m_prop_vol". If the "InclEmpty" flag is set,  the result
    // also incorporates the "empty" CoM and MoI:
    //
    constexpr void GetPropParams
    (
      Mass  a_prop_mass;  // Curr mass of propellant in this TrCone
      bool  a_with_empty, // Incorporate "empty" CoM and MoI in the results?
      Len   a_com [3],    // Output
      MoI   a_mois[3]     // Output
    )
    const
    {
      // Checks:
      assert(a_com != nullptr && a_mois != nullptr);

      // Check the limits of the Propellant Mass:
      assert(!IsNeg(a_prop_mass) && a_prop_mass <= m_propCap);

      // Also, there is little point in calling this method if the "empty" mass
      // and MoIs atr not final yet:
      assert(m_massIsFinal);

      // The Propellant Volume: Make sure it is within the limits to avoid
      // rounding errors:
      Vol propVol = std::min(a_prop_mass / m_rho, m_vol);

      // The Propellant Level. XXX: We assume that the propellant surface is
      // always orthogonal to the TrCone rotation axis:
      Len l =
        IsZero(m_deltaR):
        ? // R==r: The simplest and the most common case: A Cylinder:
          double(propVol / m_vol) * m_h

        : // General case: Solving a cubic equation by the Cardano formula;
          // since Vol'(l) > 0, there is only 1 root:
          (CbRt(m_clVol * propVol - m_Rh3) + m_Rh) / m_deltaR;

      // Mathematically, we must have 0 <= l <= h; but enforce that interval
      // explicitly to prevent rounding errors:
      l = std::max(std::min(l, m_h), 0.0_m);

      // Radius of the propellant surface (at the level "l"). Again, enforce
      // its range:
      double x  = double(l / m_h);
      assert(0 <= x && x <= 1.0);     // 0 --> R;  1 --> r:
      Len    rl = (1.0-x) * m_R + x * m_r;

      // MoIs of the Priopellant: similar to the "empty" case in the Ctor, but
      // in 3D:
      // NB: The formulas for Jx, Jin, Jort and MoIs are similar to those  the
      // 2D case (with SurfArea -> PropVol), but the "low-level" J0, J1, K are
      // somewhat different. BEWARE: use "propVol", not the total "m_vol"!
      //
      Len2 l2    = Sqr (l);
      Len3 l3    = l2 * l;
      Len5 JP0   = (   (m_JP05 * l + m_JP04) * l + m_JP03) * l3;
      Len5 JP1   = (((((m_JP15 * l + m_JP14) * l + m_JP13) * l) + m_JP12)
                               * l + m_JP11) * l;
      Len4 KP    = (  ((m_KP4  * l + m_KP3 ) * l + m_KP2 ) * l2;
      assert(!IsPos(KP));

      Len5 JPx   =  m_sinA2 * JP0  + (1.0    + m_cosA2)     * JP1   +
                    m_yz1 * (m_yz1 * propVol + 2.0 * m_sinA * K);
      assert(IsPos(JPx));

      Len5 JPin  =  m_cosA2 * JP0  + (1.0    + m_sinA2)     * JP1   +
                    m_x1  * (m_x1  * propVol + 2.0 * m_cosA * KP);
      assert(IsPos(JPin));

      Len5 JPort =  JP0 + JP1   + (Sqr(m_x1)  + Sqr(m_yz1)) * propVol +
                    2.0 * (cosA * m_x1 + sinA * m_yz1) * KP;
      assert(IsPos(JPort));

      a_mois[0]  =  m_rho * JPx;
      a_mois[1]  =  m_rho * (m_inXY ? JPin  : JPort);
      a_mois[2]  =  m_rho * (m_inXY ? JPort : JPin);

      // CoM of the Propellant:
      // "xiC" is the "local" co-ord of the CoM (xi=0 at x=x1), so it is < 0.
      // NB:   if propVol=0, then l=0 and KP=0 as well, so the generic formula
      // for "xiC" does not work anymore; but the limit is obviously xiC=0:
      Len xiC    = (LIKELY(!IsZero(propVol))) ? KP : 0.0_m;
      assert(!IsPos(xiC));
      a_com[0]   =  m_base1[0]           + cosA * xiC;
      a_com[1]   =  m_inXY ? (m_base1[1] + sinA * xiC) : 0.0_m;
      a_com[2]   =  m_inXZ ? (m_base1[2] + sinA * xiC) : 0.0_m;

      // Finally, if "empty" CoM and MoIs are to be taken into account, do so
      // (we made sure they are final yet):
      //
      if (a_with_empty)
      {
        // MoIs are simpli additive:
        a_mois[0] += m_MoIs[0];
        a_mois[1] += m_MoIs[1];
        a_mois[2] += m_MoIs[2];

        // CoM co-ords are weighted averages:
        Mass   totalMass = m_mass + a_prop_mass;
        assert(IsPos(totalMass));
        double fP        = double(a_prop_mass / totalMass); // Prop  Mass Frac
        double fE        = 1.0 - fP;                        // Empty Mass Frac
        a_com[0]  = fP * a_com[0] + fE * m_CoM[0];
        a_com[1]  = fP * a_com[1] + fE * m_CoM[1];
        a_com[2]  = fP * a_com[2] + fP * m_CoM[2];
      }
      // All Done!
    }
  };

  //=========================================================================//
  // "SpherSegm":                                                            //
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
  // (x0, y0, z0) are the co-ords of the base center (NOT of the pole!), and
  // either "z0" or "y0" must be 0;  in the former case, the Segment's rotat-
  // ion axis is assumed to lie in the XY plane, in the latter -- in XZ.
  // The spherical segment may be facing towards the positive direction   of
  // the OX axis (ie the X-coord of the pole is > x0), or otherise, as given
  // by the "is_pos_facing" flag.
  // When alpha=0, the axis of the Spherical Segment coincides with OX, which
  // is the most common case:
  //
/*
  class SpherSegm
  {

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    SpherSegm
    (
      bool    a_facing_up,
      Len     a_x0,           // X-coord of the base center (NOT of the pole!)
      Len     a_y0,
      Len     a_z0,
      double  a_alpha,        // Dimension-less, in (-Pi/2 .. Pi/2)
      Len     a_d,            // Base diameter
      Len     a_h,            // Height: 0 < a_h <= a_d/2.0
      Mass    a_empty_mass    // 0 if the mass is not known yet
    )
    {
      // NB: We must always have 0 < h <= R:
      assert(IsPos(a_d) && IsPos(a_h) && a_h <= a_d/2.0 &&
             Abs(a_alpha) < M_PI_2);

      // The mass may or may not be given:
      bool hasMass = IsPos(a_empty_mass);

      // At least one of "z0", "y0" must be 0:
      bool inXY    = IsZero(a_z0);
      bool inXZ    = IsZero(a_y0);
      assert(inXY || inXZ);

      double cosA  = Cos(a_alpha);
      double sinA  = Sin(a_alpha);
      double tanA  = sinA / cosA;

      // X- and Z-coords of the Pole; the Y-coord is not used:
      double sgn   = a_is_pos_facing ? 1.0 : -1.0;
      auto   xP    = a_x0 +  sgn  * cosA * a_h;
      auto   zP    = inXZ ? (a_z0 + sgn  * sinA * a_h) : 0.0_m;
      auto   r     = a_d / 2.0;                    // Base   radius
      auto   R     = (Sqr(r) / a_h + a_h) / 2.0;   // Sphere radius
      // Therefore, R >= r >= h, but we cannot formally assert this due to pos-
      // sible rounding effects...

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
      auto surfDens = hasMass ? (a_empty_mass / S) : SurfDens(1.0);
      auto M        = hasMass ?  a_empty_mass      : (S * surfDens);
      auto MoIY     = L4 * surfDens;
      auto MoIX     = MoIY;
      auto MoIZ     = MoIY;

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

      return ConstrElement(M, S, V, xC, yC, zC, MoIX, MoIY, MoIZ, hasMass);
    }

    //-----------------------------------------------------------------------//
    // Shell -- Spherical Segment:                                           //
    //-----------------------------------------------------------------------//
    // "d" is the base diameter (NOT the sphere diameter!), "h" is the height;
    // h <= d/2:
    constexpr static Area SideSurfArea_SpherSegm(Len a_d, Len a_h)
    {
      assert(IsPos(a_h) && a_h <= a_d/2.0);
      return M_PI * (Sqr(a_d/2.0) + Sqr(a_h));
    }

    constexpr static Vol  Volume_SpherSegm      (Len a_d, Len a_h)
    {
      assert(IsPos(a_h) && a_h <= a_d/2.0);
      return M_PI * a_h * (Sqr(a_d)/8.0 + Sqr(a_h)/6.0);
    }
*/
}
// End namespace SpaceBallistics

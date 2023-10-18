// vim:ts=2:et
//===========================================================================//
//                              "ConstrElement.h":                           //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include <cmath>
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // Construction Element:                                                   //
  //=========================================================================//
  // General characteristics of standard construction elements:
  // Currently, 2D Shells (Truncated Cones and Spherical Segments), 3D Volumes
  // (of the same shape as Shells) and Point Masses are supported.  For them, 
  // we can compute Centers of Masses (CoM) and Moments of Inertia (MoI) wrt
  // the OY axis (ALWAYS!).
  // XXX:
  // Here we do not explicitly specify the co-ord system OXYZ (typically a cer-
  // tain embedded one) in which the CoM  and the MoI (wrt to the OY axis)  are
  // given. It is assumed that the co-ord system OXYZ is uniquely determined by
  // the context:
  //
  class ConstrElement
  {
  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    Mass  m_mass;         // Mass
    Area  m_surfArea;     // Surface Area
    Vol   m_vol;          // Notional volume (with imaginary bases added)
    Len   m_CoM[3];       // (X,Y,Z) co-ords of the Center of Masses
    MoI   m_MoIY;         // Moment of Inertia wrt the OY axis
    bool  m_massIsFinal;  // If not set, "m_mass" and "m_MoY" are not valid

  public:
    //=======================================================================//
    // Ctors:                                                                //
    //=======================================================================//
    // Default Ctor:
    ConstrElement();

    // Copy Ctor:
    ConstrElement(ConstrElement const& a_right);

  private:
    // Non-Default Ctor: For internal use only. NB: The "a_mass" param may be
    // final or provisional, according to the  "a_mass_is_final" param:
    ConstrElement
      (Mass a_mass,  Area a_surf_area, Vol a_vol,
       Len  a_xc,    Len  a_yc,        Len a_zc,
       MoI  a_moi_y, bool a_mass_is_final);

  public:
    //-----------------------------------------------------------------------//
    // "SetMass":                                                            //
    //-----------------------------------------------------------------------//
    // XXX: A "setter" for the Mass (invoked when we know the final mass, sets
    // it and pro-rates the mass-dependent params). After that, the mass  will
    // be final and cannot be changed anymore:
    //
    void SetMass(Mass a_mass);

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    // Geometrical params are always available, even if the Mass is not final
    // yet:
    Area         GetSurfArea() const { return m_surfArea; }
    Vol          GetVol()      const { return m_vol;      }
    typedef Len  Point[3];
    Point const& GetCoM()      const { return m_CoM;      }

    // But the Mass and MoI may or may not be available:
    Mass GetMass() const
    {
      if (UNLIKELY(!m_massIsFinal))
        throw std::invalid_argument("ConstrElement::Mass: Mass not set yet");
      return m_mass;
    }

    MoI GetMoIY() const
    {
      if (UNLIKELY(!m_massIsFinal))
        throw std::invalid_argument("ConstrElement::MoIY: Mass not set yet");
      return m_MoIY;
    }

    //=======================================================================//
    // Generators:                                                           //
    //=======================================================================//
    // NB: They are "static" factory functions used instead of various non-def-
    // ault ctors, for the sake of clarity:
    //-----------------------------------------------------------------------//
    // "ShellTrCone":                                                        //
    //-----------------------------------------------------------------------//
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
    static ConstrElement ShellTrCone
    (
      Len     a_x0,
      Len     a_y0,
      Len     a_z0,
      double  a_alpha,        // Dimension-less, in (-Pi/2 .. Pi/2)
      Len     a_d0,           // Base diameter at x0
      Len     a_d1,           // Base diameter at the other section (x1 > x0)
      Len     a_h,            // Must be > 0
      Mass    a_mass          // 0 if the mass is not known yet
    );

    //-----------------------------------------------------------------------//
    // "ShellSpherSegm":                                                     //
    //-----------------------------------------------------------------------//
    // Similar to "ShellTrCone", but for a Spherical Segment of the base diam-
    // eter "d" and the height (from the base plane to the pole) "h", where we
    // assume h <= d/2 (the equality corresponds to the case of a HemiSphere).
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
    static ConstrElement ShellSpherSegm
    (
      bool    a_is_pos_facing,
      Len     a_x0,           // X-coord of the base center (NOT of the pole!)
      Len     a_y0,
      Len     a_z0,
      double  a_alpha,        // Dimension-less, in (-Pi/2 .. Pi/2)
      Len     a_d,            // Base diameter
      Len     a_h,            // Height: 0 < a_h <= a_d/2.0
      Mass    a_mass          // 0 if the mass is not known yet
    );

    //-----------------------------------------------------------------------//
    // "PointMass":                                                          //
    //-----------------------------------------------------------------------//
    // A psoitive mass concentrated in the (x0, y0, z0) point:
    //
    static ConstrElement PointMass(Len a_x0, Len a_y0, Len a_z0, Mass a_mass);

    //=======================================================================//
    // Addition:                                                             //
    //=======================================================================//
    // NB: The operands must both have "m_massIsFinal" set or unset, otherwise
    // an exception is thrown:
    //
    ConstrElement& operator+=(ConstrElement const& a_right);

    ConstrElement  operator+ (ConstrElement const& a_right) const
    {
      ConstrElement res = *this;
      res += a_right;
      return res;
    }

    //=======================================================================//
    // UTILS: Side Surface Areas and Volumes of Construction Elements:       //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Shell -- Truncated Cone:                                              //
    //-----------------------------------------------------------------------//
    // "d0" and "d1" are the diameters of the bases, "h" is the height:
    //
    constexpr static Area SideSurfArea_TrCone(Len a_d0, Len a_d1, Len a_h)
    {
      assert(!(IsNeg(a_d0) || IsNeg(a_d1) || (IsZero(a_d0) && IsZero(a_d1))) &&
             IsPos(a_h));
      return M_PI_2 * a_h * (a_d0 + a_d1);
    }

    constexpr static Vol  Volume_TrCone      (Len a_d0, Len a_d1, Len a_h)
    {
      assert(!(IsNeg(a_d0) || IsNeg(a_d1) || (IsZero(a_d0) && IsZero(a_d1))) &&
             IsPos(a_h));
      return (M_PI/12.0) * a_h * (Sqr(a_d0) + a_d0 * a_d1 + Sqr(a_d1));
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
  };
}
// End namespace SpaceBallistics

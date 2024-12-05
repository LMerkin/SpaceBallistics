// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/PhysForces/EarthRotationModel.h":            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include "SpaceBallistics/Utils.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "EarthRotationModel" Class:                                             //
  //=========================================================================//
  // An object of this class  describes  the orientation  of the Earth rotation
  // axis and the RA origin (True Dynamical Equinox), for a given Time (Epoch),
  // relative to the GCRS. NB: For the Nutations part, DE440T numerical data are
  // used.
  // The orientatation is represented as a certain rotation matrix  (expressing
  // GCRS co-ords via the CRS-of-Epoch co-ords) and changes very slowly with the
  // Epoch, due to Precession and Nutation. If the Epoch is J2000.0, the matrix
  // is equal to the Unit Matrix I within the accuracy implemented here.
  // The "EarthRotationModel" object constructed then acts as a function trans-
  // forming GTRS co-ords (for a given TT) to GCRS.
  //
  class EarthRotationModel
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // "M" is the (CRS_Epoch -> GCRS) "direct" conversion matrix, ie the
    // Precession-Nutation matrx (P*N) such that
    // r_GCRS = M * r_[CRS_Epoch] ;
    // XXX: BiasCorrection and PolarMotion effects are currently deemed to be
    // too small, and not taken into account:
    double  m_M  [3][3];

    // The Inverse of "M", for (GCRS -> CRS_Epoch) "inverse" conversions, s.t.
    // r_[CRS_Epoch] = invM * r_GCRS, invM = invN * invP:
    double m_invM[3][3];

    // The Equation of the Origins: Ee = ERA - GAST:
    Angle  m_Ee;

    // DeltaT = TT - UT1:
    Time   m_DeltaT;

    // Default Ctor is deleted: without the Epoch, this model makes no sense:
    EarthRotationModel() = delete;

  public:
    //-----------------------------------------------------------------------//
    // Obliquity of the Ecliptic @ J2000.0:                                  //
    //-----------------------------------------------------------------------//
    constexpr static Angle_arcSec Eps0 = 84381.406_arcSec;

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing the model for a given Year (may be fractional); the exact
    // TimeScale is presumably TT, but that does not matter here  because the
    // result depends on the Epoch VERY WEAKLY:
    // NB: It is not a "constexpr" because requires DE440T invocations. Then
    // there is no point in making any methods of this class "constexpr", ei-
    // ther:
    EarthRotationModel(Time_jyr a_epoch);

    //-----------------------------------------------------------------------//
    // ITRS -> GCRS Conversion (of any "Vector3D"):                          //
    //-----------------------------------------------------------------------//
    template<typename  DQ, Body B = Body::UNDEFINED>
    Vector3D<DQ, GCRS, B>  ToGCRS
    (
      TT                           a_tt,
      Vector3D<DQ, ITRS, B> const& a_terr
    )
    const
    {
      // NB: We must use GAST, not just ERA, because the former is consistent
      // with the origin of RA (Dynamic Equinox of "a_epoch"):
      // res = m_M * R3(-GAST) * a_terr:
      Angle    gast   =  GAST(a_tt);
      double   T[3][3];
      MkMtsR3(-gast, T, nullptr);

      DQ       tmp[3];
      MVMult3(T,   a_terr.GetArr(), tmp);

      Vector3D<DQ, GCRS, B> res;
      MVMult3(m_M, tmp,  res.GetArr());
      return  res;
    }

    //-----------------------------------------------------------------------//
    // GCRS <-> ITRS Conversion (of any "Vector3D"):                         //
    //-----------------------------------------------------------------------//
    template<typename  DQ, Body B = Body::UNDEFINED>
    Vector3D<DQ, ITRS, B>  ToITRS
    (
      TT                           a_tt,
      Vector3D<DQ, GCRS, B> const& a_geo
    )
    const
    {
      // NB: Again, using GAST:
      // res = R3(GAST) * m_invM * a_geo:
      Angle   gast   =  GAST(a_tt);
      double  T[3][3];
      MkMtsR3(gast, T, nullptr);

      DQ      tmp [3];
      MVMult3(m_invM,  a_geo.GetArr(), tmp);

      Vector3D<DQ, ITRS, B> res;
      MVMult3(T, tmp,  res.GetArr());
      return  res;
    }

    //-----------------------------------------------------------------------//
    // Auxiliary Function: GAST (Greenwich Apparent Siderial Time) from TT:  //
    //-----------------------------------------------------------------------//
    // NB: GAST is an a Angle, not Time:
    //
    Angle GAST(TT a_tt) const;

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    Angle GetEe()     const { return m_Ee;     }
    Time  GetDeltaT() const { return m_DeltaT; }
  };
}
// End namespace SpaceBallistics

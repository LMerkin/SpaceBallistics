// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/PhysForces/EarthRotationModel.h":            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include <utility>

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
  public:
    //-----------------------------------------------------------------------//
    // Obliquity of the Ecliptic @ J2000.0:                                  //
    //-----------------------------------------------------------------------//
    constexpr static Angle_arcSec Eps0 = 84381.406_arcSec;

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    Time_jyr  m_ermEpoch;     // Epoch of Precession-Nutation

    // "PN" is the (CRS_Epoch -> GCRS) "direct" conversion matrix, ie the
    // Precession-Nutation matrx (P*N) such that
    // r_GCRS = PN * r_[CRS_Epoch] ;
    // XXX: BiasCorrection and PolarMotion effects are currently deemed to be
    // too small, and not taken into account:
    Mtx33     m_PN;

    // The Inverse of "PN", for (GCRS -> CRS_Epoch) "inverse" conversions, s.t.
    // invPN = N^(-1) * P^(-1),
    // r_[CRS_Epoch] = invPN * r_GCRS:
    Mtx33     m_invPN;

    // The Equation of the Origins: Ee = ERA - GAST:
    Angle     m_Ee;

    // DeltaT = TT - UT1:
    Time      m_DeltaT;

    // Default Ctor is deleted: without the Epoch, this model makes no sense:
    EarthRotationModel() = delete;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // Constructing the model for a given Year (may be fractional); the exact
    // TimeScale is presumably TT, but that does not matter here  because the
    // result depends on the Epoch VERY WEAKLY:
    // NB: They are not "constexpr"s because require DE440T invocations:
    //
    EarthRotationModel(Time_jyr  a_erm_epoch);
    EarthRotationModel(TT        a_erm_epoch);

    // The Ctor with an integral YearNumber arg uses the Analytical Nutations
    // Model (IAU2000B), and is therefore a "constexpr".  Intended for use in
    // "GeoCDynEqFixCOS";  then "EarthRotationModel.hpp"  must  be incuded to
    // provide the body of this Ctor:
    //
    constexpr EarthRotationModel(int a_erm_epoch_year);

    //-----------------------------------------------------------------------//
    // ITRS -> GCRS Conversion (of any "Vector3D"):                          //
    //-----------------------------------------------------------------------//
    template<typename  DQ, Body B = Body::UNDEFINED>
    constexpr Vector3D<DQ, GCRS, B>  ToGCRS
    (
      TT                           a_tt,
      Vector3D<DQ, ITRS, B> const& a_terr
    )
    const
    {
      // NB: We must use GAST, not just ERA, because the former is consistent
      // with the origin of RA (Dynamic Equinox of "m_ermEpoch):
      // res  = m_PN * R3(-GAST) * a_terr:
      Mtx33 T = Mtx33::MkR3(-GAST(a_tt));

      DQ tmp[3];
      T.MVMult(a_terr.GetArr(), tmp);

      Vector3D<DQ, GCRS, B> res;
      m_PN.MVMult(tmp,   res.GetArr());
      return res;
    }

    //-----------------------------------------------------------------------//
    // GCRS -> ITRS Conversion (of any "Vector3D"):                          //
    //-----------------------------------------------------------------------//
    template<typename  DQ, Body B = Body::UNDEFINED>
    constexpr Vector3D<DQ, ITRS, B>  ToITRS
    (
      TT                           a_tt,
      Vector3D<DQ, GCRS, B> const& a_geo
    )
    const
    {
      // NB: Again, using GAST:
      // res  = R3(GAST) * m_invPN * a_geo:
      Mtx33 T = Mtx33::MkR3(GAST(a_tt));

      DQ tmp[3];
      m_invPN.MVMult(a_geo.GetArr(), tmp);

      Vector3D<DQ,  ITRS, B> res;
      T.MVMult(tmp, res.GetArr());
      return res;
    }

    //-----------------------------------------------------------------------//
    // Auxiliary Function: GAST (Greenwich Apparent Siderial Time) from TT:  //
    //-----------------------------------------------------------------------//
    // NB: GAST is an a Angle, not Time:
    //
    constexpr Angle GAST(TT a_tt) const
    {
      // UT1_Since_Epoch = TT_Since_Epoch  - DeltaT:
      Time ut1  = a_tt.GetTimeSinceEpoch() - m_DeltaT;

      // EarthRotationAngle (ERA) is a linear function of UT1_Since_Epoch expr-
      // essed in JDs:
      double t  = double(To_Time_day(ut1) / 1.0_day);
      Angle ERA =
        TwoPi<double> *
        Angle(0.779'057'273'264 + 1.002'737'811'911'354'480 * t);

      // Finally, GAST = ERA + Ee:
      return ERA + m_Ee;
    }

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    constexpr Time_jyr     GetERMEpoch() const { return m_ermEpoch; }
    constexpr Angle        GetEe      () const { return m_Ee;       }
    constexpr Time         GetDeltaT  () const { return m_DeltaT;   }
    constexpr Mtx33 const& GetPN      () const { return m_PN;       }
    constexpr Mtx33 const& GetInvPN   () const { return m_invPN;    }

    // Furthermore, the columns of PN are actually the unit vectors  of the
    // "GeoCDynEqFixCOS" (with the Dynamic Equator and Mean Ecliptic of the
    // "ERMEpoch") in the GCRS system, so we can return them one-by-one as
    // GCRS vectors (assuming they are pos vectors of unit length):
    //
    // Col0: The X axis of "GeoCDynEqFixCOS": Points to the Equinox of the
    // "ERMEpoch":
    constexpr PosKV_GCRS<> GetGeoCDynEqFixX() const
    {
      return  PosKV_GCRS<>
              {m_PN(0,0) * 1.0_km, m_PN(1,0) * 1.0_km, m_PN(2,0) * 1.0_km};
    }

    // Col1: The Y axis of "GeoCDynEqFixCOS": Just complements the X and Z
    //       axes to the right-oriented XYZ frame:
    constexpr PosKV_GCRS<> GetGeoCDynEqFixY() const
    {
      return  PosKV_GCRS<>
              {m_PN(0,1) * 1.0_km, m_PN(1,1) * 1.0_km, m_PN(2,1) * 1.0_km};
    }

    // Col2: The Z axis of "GeoCDynEqFixCOS": Points to the (Equatorial) North
    //       Pole  of the  "ERMEpoch":
    //
    constexpr PosKV_GCRS<> GetGeoCDynEqFixZ() const
    {
      return  PosKV_GCRS<>
              {m_PN(0,2) * 1.0_km, m_PN(1,2) * 1.0_km, m_PN(2,2) * 1.0_km};
    }

  private:
    //-----------------------------------------------------------------------//
    // Utils (for internal use):                                             //
    //-----------------------------------------------------------------------//
    // Time from the Epoch (J2000.0) to "a_t" RELATIVE to Julian Centuries, as
    // a mere "double": For use in various polynomial expansions:
    //
    constexpr static double GetJCYsSinceEpoch(Time_jyr a_t);

    // Earth Precession Mtx:
    constexpr static Mtx33 MkPrecMtx(double a_T);

    // Earth Nutations Angles: (d(Psi), d(Eps), cos(eps)) using the Analytical
    // Model (IAU2000B), hence "constexpr":
    constexpr static std::pair<Angle,Angle> GetNutAnglesAnalyt(double a_T);

    // Earth Nutation Matrix, via Nutation Angles. Also returns Cos(eps):
    constexpr static std::pair<Mtx33,double> MkNutMtx
      (double a_T, Angle a_dpsi, Angle a_deps);

    // DeltaT = TT - UT1:
    constexpr static Time GetDeltaT(Time_jyr a_t);

  public:
    //-----------------------------------------------------------------------//
    // For External Testability:                                             //
    //-----------------------------------------------------------------------//
    // Earth Nutations Angles: (d(Psi), d(Eps), cos(eps)) using DE440T;  thus,
    // NON-"constexpr":
    //
    static           std::pair<Angle,Angle> GetNutAnglesDE440T(Time_jyr a_t);

    // An overload of "GetNutAnglesAnalyt" for external testability;  requires
    // including the  "EarthRotationModel.hpp" for use:
    //
    constexpr static std::pair<Angle,Angle> GetNutAnglesAnalyt(Time_jyr a_t);
  };
}
// End namespace SpaceBallistics

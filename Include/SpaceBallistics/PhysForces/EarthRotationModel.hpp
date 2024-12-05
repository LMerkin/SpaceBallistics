// vim:ts=2:et
//===========================================================================//
//            "SpaceBallistics/PhysForces/EarthRotationModel.hpp":           //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysForces/EarthRotationModel.h"
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "EarthRotationModel" Non-Default Ctor:                                  //
  //=========================================================================//
  EarthRotationModel::EarthRotationModel(Time_jyr a_epoch)
  {
    // Time from the Epoch (J2000.0) in Julian Years:
    double T  = (a_epoch - 2000.0_jyr).Magnitude() / 100.0;
    double T2 = T * T;

    //-----------------------------------------------------------------------//
    // Precession:                                                           //
    //-----------------------------------------------------------------------//
    Angle omegaA = To_Angle(Eps0 - 0.027754_arcSec * T + 0.0512623_arcSec * T2);
    Angle psiA   = To_Angle(    5038.481507_arcSec * T - 1.0790069_arcSec * T2);
    Angle chiA   = To_Angle(      10.556403_arcSec * T - 2.3814292_arcSec * T2);

    double R1pEps0[3][3], R3mPsiA[3][3], R1mOmgA[3][3], R3pChiA[3][3],
           R1mEps0[3][3], R3pPsiA[3][3], R1pOmgA[3][3], R3mChiA[3][3],
           Tmp0   [3][3], Tmp1   [3][3], P      [3][3], invP   [3][3];

    MkMtsR1(To_Angle(Eps0), R1pEps0, R1mEps0);
    MkMtsR3(-psiA,          R3mPsiA, R3pPsiA);
    MkMtsR1(-omegaA,        R1mOmgA, R1pOmgA);
    MkMtsR3( chiA,          R3pChiA, R3mChiA);

    // P    = R3(chiA)  * R1(-omegaA) * R3(-psiA)  * R1(eps0):
    MtxMult3(R3pChiA,       R1mOmgA, Tmp0);
    MtxMult3(Tmp0,          R3mPsiA, Tmp1);
    MtxMult3(Tmp1,          R1pEps0, P);

    // invP = R1(-eps0) * R3(psiA)    * R1(omegaA) * R3(-chiA);
    MtxMult3(R1mEps0,       R3pPsiA, Tmp0);
    MtxMult3(Tmp0,          R1pOmgA, Tmp1);
    MtxMult3(Tmp1,          R3mChiA, invP);

    //-----------------------------------------------------------------------//
    // Nutation:                                                             //
    //-----------------------------------------------------------------------//
    // Use DE440T data if available, otherwise assume N=Id.
    //
    Angle  epsA  =  To_Angle(Eps0 - 46.836769_arcSec * T);  // T2 term is tiny
    Angle  dPsi,    dEps;         // Initially 0; properly inited below
    double N[3][3], invN[3][3];   // Initialised below

    if (DE440T::Bits::FromY <= a_epoch && a_epoch < DE440T::Bits::ToY)
    {
      // Strictly speaking, we must interpret "a_epoch" as TT and convert it to
      // TDB:
      TDB     tdb{TT{a_epoch}};
      Angle   nutAngles[2];       // [dPsi, dEps]
      DE440T::GetEarthNutations(tdb, nutAngles);

      dPsi      = nutAngles[0];
      dEps      = nutAngles[1];
      Angle eps = epsA + dEps;

      // N    = R1(-eps) * R3(-dPsi) * R1(epsA);
      double R1mEps[3][3], R3mDPsi[3][3], R1pEpsA[3][3],
             R1pEps[3][3], R3pDPsi[3][3], R1mEpsA[3][3];

      MkMtsR1(-eps,     R1mEps,  R1pEps);
      MkMtsR3(-dPsi,    R3mDPsi, R3pDPsi);
      MkMtsR1( epsA,    R1pEpsA, R1mEpsA);

      MtxMult3(R1mEps,  R3mDPsi, Tmp0);
      MtxMult3(Tmp0,    R1pEpsA, N);

      // invN = R1(-epsA) * R3(dPsi) * R1(eps):
      MtxMult3(R1mEpsA, R3pDPsi, Tmp0);
      MtxMult3(Tmp0,    R1pEps,  invN);
    }
    else
      // N = invN = Id:
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        N[i][j] = invN[i][j] = (i==j) ? 1.0 : 0.0;

    //-----------------------------------------------------------------------//
    // Over-All:                                                             //
    //-----------------------------------------------------------------------//
    // XXX: Bias Correction and Polar Motion are considered to be too small,
    // and not used. Thus:
    MtxMult3(P,    N,    m_M);
    MtxMult3(invN, invP, m_invM);

    // Check:
    DEBUG_ONLY
    (
      MtxMult3(m_M, m_invM, Tmp0);
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        assert(DimTypes::Bits::CEMaths::ApproxEqual
                (Tmp0[i][j], (i==j) ? 1.0 : 0.0));
    )
    //-----------------------------------------------------------------------//
    // An Approximate Equation of the Origins Ee = ERA - GAST:               //
    //-----------------------------------------------------------------------//
    m_Ee = dPsi * Cos(double(epsA)) / 15.0;

    //-----------------------------------------------------------------------//
    // DeltaT = TT - UT1:                                                    //
    //-----------------------------------------------------------------------//
    // We will use the value from the Maple-generated Cubic Spline approximation
    // using the following table. XXX: "DeltaT" is linearly exptrapolated beyond
    // the range 1620..2024:
    // |---------------|
    // | Year | ΔT (s) |
    // |------|--------|
    // | 1620 | 124.0  |
    // | 1700 |   9.2  |
    // | 1750 |  13.4  |
    // | 1800 |  14.2  |
    // | 1850 |   7.3  |
    // | 1900 |  -2.8  |
    // | 1950 |  29.1  |
    // | 2000 |  63.83 |
    // | 2010 |  66.07 |
    // | 2015 |  68.10 |
    // | 2020 |  70.10 |
    // | 2024 |  72.58 |
    // |---------------|
    // XXX: For compactness, we use the dimension-less "y" below. This is OK as
    // it is not exported outside:
    //
    double y = a_epoch.Magnitude();
    m_DeltaT =
      Time
      (
        (y < 1700.0)
        ? 3253.47435750767  - 1.93177429475782   * y
          + 0.0000776209835559110 * Cube(y - 1620.0)
        :
        (y < 1750.0)
        ? 759.667397823397  - 0.441451410484351  * y
          + 0.0186290360534185    * Sqr (y - 1700.0)
          - 0.000162400156874629  * Cube(y - 1700.0)
        :
        (y < 1800.0)
        ? -342.639282021111 + 0.203451018297778  * y
          - 0.00573098747777592   * Sqr (y - 1750.0)
          + 0.0000396393422364073 * Cube(y - 1750.0)
        :
        (y < 1850.0)
        ? 144.434792872167  - 0.0723526627067592 * y
          + 0.000214913857685184  * Sqr (y - 1800.0)
          - 0.0000305572120710000 * Cube(y - 1800.0)
        :
        (y < 1900.0)
        ? 525.374679820870  - 0.280040367470741  * y
          - 0.00436866795296482   * Sqr (y - 1850.0)
          + 0.000118589506047593  * Cube(y - 1850.0)
        :
        (y < 1950.0)
        ? -330.576851920473 + 0.172514132589723  * y
          + 0.0134197579541741    * Sqr (y - 1900.0)
          - 0.0000822008121193708 * Cube(y - 1900.0)
        :
        (y < 2000.0)
        ? -1721.96848236811 + 0.897983837111850  * y
          + 0.00108963613626847   * Sqr (y - 1950.0)
          - 0.000103146257570110  * Cube(y - 1950.0)
        :
        (y < 2010.0)
        ? -402.871037925750 + 0.233350518962875  * y
          - 0.0143823024992480    * Sqr (y - 2000.0)
          + 0.00134472506029605   * Cube(y - 2000.0)
        :
        (y < 2015.0)
        ? -635.665194004126 + 0.349121987066729  * y
          + 0.0259594493096334    * Sqr (y - 2010.0)
          - 0.00291676934459586   * Cube(y - 2010.0)
        :
        (y < 2020.0)
        ? -717.666940326524 + 0.389958779318374  * y
          - 0.0177920908593045    * Sqr (y - 2015.0)
          + 0.00396006699912594   * Cube(y - 2015.0)
        :
          -958.166649232744 + 0.509042895659774  * y
          + 0.0416089141275846    * Sqr (y - 2020.0)
          - 0.00346740951063205   * Cube(y - 2020.0)
      );
  }

  //=========================================================================//
  // "GAST":                                                                 //
  //=========================================================================//
  Angle EarthRotationModel::GAST(TT a_tt) const
  {
    // UT1_Since_Epoch = TT_Since_Epoch  - DeltaT:
    Time ut1  = a_tt.GetTimeSinceEpoch() - m_DeltaT;

    // EarthRotationAngle (ERA) is a linear function of UT1_Since_Epoch
    // expressed in JDs:
    double t  = To_Time_day(ut1).Magnitude();
    Angle ERA =
      TwoPi<double> * Angle(0.779'057'273'264 + 1.002'737'811'911'354'480 * t);

    // Finally, GAST = ERA + Ee:
    return ERA + m_Ee;
  }
}
// End namespace SpaceBallistics

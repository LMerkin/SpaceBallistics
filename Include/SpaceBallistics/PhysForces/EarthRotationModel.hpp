// vim:ts=2:et
//===========================================================================//
//            "SpaceBallistics/PhysForces/EarthRotationModel.hpp":           //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysForces/EarthRotationModel.h"
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/Utils.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // Time from the Epoch (J2000.0) RELATIVE to Julian Centuries:             //
  //=========================================================================//
  constexpr double EarthRotationModel::GetJCYsSinceEpoch(Time_jyr a_t)
  {
    // XXX: "a_t" is actually the Year Number; it is OK to use the following
    // "low-accuracy" formula (instead of exact Epoch = JD_TT 2451545.0):
    return double((a_t - 2000.0_jyr) / Time_jyr(100.0));
  }
    
  //=========================================================================//
  // Precession Mtx and its Inverse:                                         //
  //=========================================================================//
  constexpr std::pair<Mtx33,Mtx33> EarthRotationModel::MkPrecMts(double a_T)
  {
    double T  = a_T;
    double T2 = T * T;
    double T3 = T * T2;

    //-----------------------------------------------------------------------//
    // Precession Angles (Capitain, Wallace & Chapront 2003):                //
    //-----------------------------------------------------------------------//
    Angle psi1 =
      To_Angle(       5038.481507_arcSec * T - 1.0790069_arcSec * T2
                    - 0.001214045_arcSec * T3);
    Angle eps1 =
      To_Angle(Eps0 - 0.027754_arcSec    * T + 0.0512623_arcSec * T2
                    - 0.00772503_arcSec  * T3);
    Angle chi  =      // Planetary Precession:
      To_Angle(       10.556403_arcSec   * T - 2.3814292_arcSec * T2
                    - 0.00121197_arcSec  * T3);

    Mtx33 R1pEps0, R3pPsi1, R1pEps1,  R3pChi,
          R1mEps0, R3mPsi1, R1mEps1,  R3mChi,
          Tmp0,    Tmp1,    P,        invP;

    // Elementary Rotation Matrices:
    MkMtsR1(To_Angle(Eps0), &R1pEps0, &R1mEps0);
    MkMtsR3(psi1,           &R3pPsi1, &R3mPsi1);
    MkMtsR1(eps1,           &R1pEps1, &R1mEps1);
    MkMtsR3(chi,            &R3pChi,  &R3mChi);

    //-----------------------------------------------------------------------//
    // P    = R1(-eps0) * R3(psi1) * R1(eps1) * R3(-chi):                    //
    //-----------------------------------------------------------------------//
    MtxMult33(R1mEps0,      R3pPsi1,  &Tmp0);
    MtxMult33(Tmp0,         R1pEps1,  &Tmp1);
    MtxMult33(Tmp1,         R3mChi,   &P);

    //-----------------------------------------------------------------------//
    // invP = R3(chi) * R1(-eps1) * R3(-psi1) * R1(eps0):                    //
    //-----------------------------------------------------------------------//
    MtxMult33(R3pChi,       R1mEps1,  &Tmp0);
    MtxMult33(Tmp0,         R3mPsi1,  &Tmp1);
    MtxMult33(Tmp1,         R1pEps0,  &invP);

    return std::make_pair(P, invP);
  }

  //=========================================================================//
  // Nutation Matrices as Functions of d(Psi) and d(Eps):                    //
  //=========================================================================//
  constexpr std::tuple<Mtx33,Mtx33,double> EarthRotationModel::MkNutMts
    (double a_T, Angle a_dpsi, Angle a_deps)
  {
    Angle eps  =
      To_Angle(((0.00200340_arcSec * a_T - 0.0001831_arcSec) * a_T
               - 46.836769_arcSec) * a_T + Eps0);
    Angle epsD = eps + a_deps;

    Mtx33 R1pEpsD,  R3pDPsi,  R1pEps,
          R1mEpsD,  R3mDPsi,  R1mEps,
          Tmp,      N,        invN;

    MkMtsR1(epsD,   &R1pEpsD, &R1mEpsD);
    MkMtsR3(a_dpsi, &R3pDPsi, &R3mDPsi);
    MkMtsR1(eps,    &R1pEps,  &R1mEps );

    // N    = R1(-eps) * R3(dPsi) * R1(epsD):
    MtxMult33(R1mEps,  R3pDPsi, &Tmp);
    MtxMult33(Tmp,     R1pEpsD, &N);

    // invN = R1(-epsD) * R3(-dPsi) * R1(eps);
    MtxMult33(R1mEpsD, R3mDPsi, &Tmp);
    MtxMult33(Tmp,     R1pEps,  &invN);

    // NB: Cos(eps) == R1pEps(1,1), but it's better compute it anew:
    return std::make_tuple(N, invN, Cos(double(eps)));
  }

  //=========================================================================//
  // Analytical Nutation Angles (IAU2000B Model): Returns (d(psi),d(eps)):   //
  //=========================================================================//
  constexpr std::pair<Angle,Angle> EarthRotationModel::GetNutAnglesAnalyt
    (double a_T)
  {
    //-----------------------------------------------------------------------//
    // Compute the Delaunay Angles:                                          //
    //-----------------------------------------------------------------------//
    double T  = a_T;
    double T2 = T * T;
    double T3 = T * T2;

    struct DelaunayCoeffs  // Up to T^3
    {
      Angle_deg     m_c0;
      Angle_arcSec  m_c1;
      Angle_arcSec  m_c2;
      Angle_arcSec  m_c3;
    };
    constexpr DelaunayCoeffs DCs[5]
    {
      // [0]: Mean Anomaly of the Moon (l):
      { 134.96340251_deg,
        1717915923.2178_arcSec, 31.8792_arcSec,   0.051635_arcSec },

      // [1]: Mean Anomaly of the Sun  (l'):
      { 357.52910918_deg,
        129596581.0481_arcSec,  -0.5532_arcSec,  -0.000136_arcSec },

      // [2]: Mean longitude of the Moon from the Ascending Node (F):
      { 93.27209062_deg,
        1739527262.8478_arcSec, -12.7512_arcSec, -0.001037_arcSec },

      // [3]: Mean Elongation of the Moon from the Sun (D):
      { 297.85019547_deg,
        1602961601.2090_arcSec, -6.3706_arcSec,   0.006593_arcSec },

      // [4] Mean Longitude of the Moon Ascending Node (Omega):
      { 125.04455501_deg,
        - 6962890.5431_arcSec,   7.4722_arcSec,   0.007702_arcSec }
    };

    // Delanauy Angles:
    Angle DAs[5];
    for (int j = 0; j < 5; ++j)
    {
      DelaunayCoeffs const& dcj   = DCs[j];
      DAs[j] = To_Angle(dcj.m_c0) +
               To_Angle(dcj.m_c1  * T + dcj.m_c2 * T2 + dcj.m_c3 * T3);
    }

    //-----------------------------------------------------------------------//
    // Coeffs for d(phi) and d(eps):                                         //
    //-----------------------------------------------------------------------//
    // Luni-Solar nutation coefficients, unit 1e-7 arcsec:
    // 5 modes, then longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin):
    //
    constexpr int    NMs     = 77;
    constexpr double NutCoeffs[NMs][11]
    {
      { 0,  0,  0,  0,  1, -172064161, -174666,  33386, 92052331,  9086, 15377},
      { 0,  0,  2, -2,  2,  -13170906,   -1675, -13696,  5730336, -3015, -4587},
      { 0,  0,  2,  0,  2,   -2276413,    -234,   2796,   978459,  -485,  1374},
      { 0,  0,  0,  0,  2,    2074554,     207,   -698,  -897492,   470,  -291},
      { 0,  1,  0,  0,  0,    1475877,   -3633,  11817,    73871,  -184, -1924},
      { 0,  1,  2, -2,  2,    -516821,    1226,   -524,   224386,  -677,  -174},
      { 1,  0,  0,  0,  0,     711159,      73,   -872,    -6750,     0,   358},
      { 0,  0,  2,  0,  1,    -387298,    -367,    380,   200728,    18,   318},
      { 1,  0,  2,  0,  2,    -301461,     -36,    816,   129025,   -63,   367},
      { 0, -1,  2, -2,  2,     215829,    -494,    111,   -95929,   299,   132},
      { 0,  0,  2, -2,  1,     128227,     137,    181,   -68982,    -9,    39},
      {-1,  0,  2,  0,  2,     123457,      11,     19,   -53311,    32,    -4},
      {-1,  0,  0,  2,  0,     156994,      10,   -168,    -1235,     0,    82},
      { 1,  0,  0,  0,  1,      63110,      63,     27,   -33228,     0,    -9},
      {-1,  0,  0,  0,  1,     -57976,     -63,   -189,    31429,     0,   -75},
      {-1,  0,  2,  2,  2,     -59641,     -11,    149,    25543,   -11,    66},
      { 1,  0,  2,  0,  1,     -51613,     -42,    129,    26366,     0,    78},
      {-2,  0,  2,  0,  1,      45893,      50,     31,   -24236,   -10,    20},
      { 0,  0,  0,  2,  0,      63384,      11,   -150,    -1220,     0,    29},
      { 0,  0,  2,  2,  2,     -38571,      -1,    158,    16452,   -11,    68},
      { 0, -2,  2, -2,  2,      32481,       0,      0,   -13870,     0,     0},
      {-2,  0,  0,  2,  0,     -47722,       0,    -18,      477,     0,   -25},
      { 2,  0,  2,  0,  2,     -31046,      -1,    131,    13238,   -11,    59},
      { 1,  0,  2, -2,  2,      28593,       0,     -1,   -12338,    10,    -3},
      {-1,  0,  2,  0,  1,      20441,      21,     10,   -10758,     0,    -3},
      { 2,  0,  0,  0,  0,      29243,       0,    -74,     -609,     0,    13},
      { 0,  0,  2,  0,  0,      25887,       0,    -66,     -550,     0,    11},
      { 0,  1,  0,  0,  1,     -14053,     -25,     79,     8551,    -2,   -45},
      {-1,  0,  0,  2,  1,      15164,      10,     11,    -8001,     0,    -1},
      { 0,  2,  2, -2,  2,     -15794,      72,    -16,     6850,   -42,    -5},
      { 0,  0, -2,  2,  0,      21783,       0,     13,     -167,     0,    13},
      { 1,  0,  0, -2,  1,     -12873,     -10,    -37,     6953,     0,   -14},
      { 0, -1,  0,  0,  1,     -12654,      11,     63,     6415,     0,    26},
      {-1,  0,  2,  2,  1,     -10204,       0,     25,     5222,     0,    15},
      { 0,  2,  0,  0,  0,      16707,     -85,    -10,      168,    -1,    10},
      { 1,  0,  2,  2,  2,      -7691,       0,     44,     3268,     0,    19},
      {-2,  0,  2,  0,  0,     -11024,       0,    -14,      104,     0,     2},
      { 0,  1,  2,  0,  2,       7566,     -21,    -11,    -3250,     0,    -5},
      { 0,  0,  2,  2,  1,      -6637,     -11,     25,     3353,     0,    14},
      { 0, -1,  2,  0,  2,      -7141,      21,      8,     3070,     0,     4},
      { 0,  0,  0,  2,  1,      -6302,     -11,      2,     3272,     0,     4},
      { 1,  0,  2, -2,  1,       5800,      10,      2,    -3045,     0,    -1},
      { 2,  0,  2, -2,  2,       6443,       0,     -7,    -2768,     0,    -4},
      {-2,  0,  0,  2,  1,      -5774,     -11,    -15,     3041,     0,    -5},
      { 2,  0,  2,  0,  1,      -5350,       0,     21,     2695,     0,    12},
      { 0, -1,  2, -2,  1,      -4752,     -11,     -3,     2719,     0,    -3},
      { 0,  0,  0, -2,  1,      -4940,     -11,    -21,     2720,     0,    -9},
      {-1, -1,  0,  2,  0,       7350,       0,     -8,      -51,     0,     4},
      { 2,  0,  0, -2,  1,       4065,       0,      6,    -2206,     0,     1},
      { 1,  0,  0,  2,  0,       6579,       0,    -24,     -199,     0,     2},
      { 0,  1,  2, -2,  1,       3579,       0,      5,    -1900,     0,     1},
      { 1, -1,  0,  0,  0,       4725,       0,     -6,      -41,     0,     3},
      {-2,  0,  2,  0,  2,      -3075,       0,     -2,     1313,     0,    -1},
      { 3,  0,  2,  0,  2,      -2904,       0,     15,     1233,     0,     7},
      { 0, -1,  0,  2,  0,       4348,       0,    -10,      -81,     0,     2},
      { 1, -1,  2,  0,  2,      -2878,       0,      8,     1232,     0,     4},
      { 0,  0,  0,  1,  0,      -4230,       0,      5,      -20,     0,    -2},
      {-1, -1,  2,  2,  2,      -2819,       0,      7,     1207,     0,     3},
      {-1,  0,  2,  0,  0,      -4056,       0,      5,       40,     0,    -2},
      { 0, -1,  2,  2,  2,      -2647,       0,     11,     1129,     0,     5},
      {-2,  0,  0,  0,  1,      -2294,       0,    -10,     1266,     0,    -4},
      { 1,  1,  2,  0,  2,       2481,       0,     -7,    -1062,     0,    -3},
      { 2,  0,  0,  0,  1,       2179,       0,     -2,    -1129,     0,    -2},
      {-1,  1,  0,  1,  0,       3276,       0,      1,       -9,     0,     0},
      { 1,  1,  0,  0,  0,      -3389,       0,      5,       35,     0,    -2},
      { 1,  0,  2,  0,  0,       3339,       0,    -13,     -107,     0,     1},
      {-1,  0,  2, -2,  1,      -1987,       0,     -6,     1073,     0,    -2},
      { 1,  0,  0,  0,  2,      -1981,       0,      0,      854,     0,     0},
      {-1,  0,  0,  1,  0,       4026,       0,   -353,     -553,     0,  -139},
      { 0,  0,  2,  1,  2,       1660,       0,     -5,     -710,     0,    -2},
      {-1,  0,  2,  4,  2,      -1521,       0,      9,      647,     0,     4},
      {-1,  1,  0,  1,  1,       1314,       0,      0,     -700,     0,     0},
      { 0, -2,  2, -2,  1,      -1283,       0,      0,      672,     0,     0},
      { 1,  0,  2,  2,  1,      -1331,       0,      8,      663,     0,     4},
      {-2,  0,  2,  2,  2,       1383,       0,     -2,     -594,     0,    -2},
      {-1,  0,  0,  0,  2,       1405,       0,      4,     -610,     0,     2},
      { 1,  1,  2, -2,  2,       1290,       0,      0,     -556,     0,     0}
    };
    //-----------------------------------------------------------------------//
    // Compute d(phi), d(eps):                                               //
    //-----------------------------------------------------------------------//
    Angle_arcSec dPhi, dEps;   // Initially 0s

    for (int i = 0; i < NMs; ++i)
    {
      Angle theta;             // Initially 0
      for (int j = 0; j < 5; ++j)
        theta += NutCoeffs[i][j] * DAs[j];

      double sinTh = Sin(double(theta));
      double cosTh = Cos(double(theta));

      // Longitude: (sin, t*sin, cos):
      dPhi += Angle_arcSec
              ((NutCoeffs[i][5] + NutCoeffs[i][6] * T) * sinTh +
                NutCoeffs[i][7]                        * cosTh);

      // Obliquity: (cos, t*cos, sin):
      dEps += Angle_arcSec
              ((NutCoeffs[i][8] + NutCoeffs[i][9] * T) * cosTh +
                NutCoeffs[i][10]                       * sinTh);
    }
    // Apply the Scaling Factors:
    dPhi *= 1e-7;
    dEps *= 1e-7;

    // Convert the results to radians:
    return std::make_pair(To_Angle(dPhi), To_Angle(dEps));
  }

  //-------------------------------------------------------------------------//
  // With the "Time_jyr" arg:                                                //
  //-------------------------------------------------------------------------//
  constexpr std::pair<Angle,Angle> EarthRotationModel::GetNutAnglesAnalyt
    (Time_jyr a_t)
  { return GetNutAnglesAnalyt(GetJCYsSinceEpoch(a_t)); }

  //=========================================================================//
  // DE440T Nutation Angles: Returns (d(psi),d(eps)):                        //
  //=========================================================================//
  // NON-"constexpr"!
  //
  std::pair<Angle,Angle> EarthRotationModel::GetNutAnglesDE440T(Time_jyr a_t)
  {
    // Strictly speaking, we must interpret "a_t" as TT and convert it to TDB:
    TDB     tdb{TT(a_t)};
    Angle   nutAngles[2];       // [dPsi, dEps]
    DE440T::GetEarthNutations(tdb, nutAngles);
    return std::make_pair(nutAngles[0], nutAngles[1]);
  }

  //=========================================================================//
  // DeltaT = TT - UT1:                                                      //
  //=========================================================================//
  constexpr Time EarthRotationModel::GetDeltaT(Time_jyr a_t)
  {
    // We will use the value from the Maple-generated Cubic Spline approximation
    // using the following table. XXX: "DeltaT" is linearly extrapolated beyond
    // the range 1620..2024:
    struct    DeltaTCoeffs
    {
      Time_jyr  m_year;
      Time      m_c0;
      Time      m_c1;
      Time      m_c2;
      Time      m_c3;
    };
    constexpr int NDTCs  = 11;

    constexpr DeltaTCoeffs DTCs[NDTCs]
    {
      { 1620.0_jyr,    124.0_sec,   -1.93177429475782_sec,
        0.0_sec, 0.0000776209835559110_sec
      },
      { 1700.0_jyr,      9.2_sec,   -0.441451410484351_sec,
        0.0186290360534185_sec,     -0.000162400156874629_sec
      },
      { 1750.0_jyr,     13.4_sec,    0.203451018297778_sec,
        -0.00573098747777592_sec,    0.0000396393422364073_sec
      },
      { 1800.0_jyr,  14.2_sec,      -0.0723526627067592_sec,
        0.000214913857685184_sec,   -0.0000305572120710000_sec
      },
      { 1850.0_jyr,   7.3_sec,      -0.280040367470741_sec,
        -0.00436866795296482_sec,    0.000118589506047593_sec
      },
      { 1900.0_jyr,  -2.8_sec,       0.172514132589723_sec,
        0.0134197579541741_sec,     -0.0000822008121193708_sec
      },
      { 1950.0_jyr,  29.1_sec,       0.897983837111850_sec,
        0.00108963613626847_sec,    -0.000103146257570110_sec
      },
      { 2000.0_jyr,  63.83_sec,      0.233350518962875_sec,
        -0.0143823024992480_sec,     0.00134472506029605_sec
      },
      { 2010.0_jyr,  66.07_sec,      0.349121987066729_sec,
        0.0259594493096334_sec,     -0.00291676934459586_sec
      },
      { 2015.0_jyr,  68.10_sec,      0.389958779318374_sec,
        -0.0177920908593045_sec,     0.00396006699912594_sec
      },
      { 2020.0_jyr,  70.10_sec,      0.509042895659774_sec,
        0.0416089141275846_sec,     -0.00346740951063205_sec
      }
    };
    for (int i = 0; i < NDTCs; ++i)
      if (i == NDTCs-1 || a_t < DTCs[i+1].m_year)
      {
        // Use this entry:
        DeltaTCoeffs const& dtci = DTCs[i];
        double y = double((a_t - dtci.m_year) / 1.0_jyr);
        return ((dtci.m_c3 * y + dtci.m_c2) * y + dtci.m_c1) * y + dtci.m_c0;
      }
    // This point should not be reached:
    assert(false);
    return 0.0_sec;
  }

  //=========================================================================//
  // "EarthRotationModel" Non-Default Ctor, with DE440T Nutations Model:     //
  //=========================================================================//
  // So this Ctor is NON-"constexpr"!
  //
  EarthRotationModel::EarthRotationModel(Time_jyr a_erm_epoch /* Assumed TT */)
  {
    m_ermEpoch = a_erm_epoch;

    //-----------------------------------------------------------------------//
    // Precession Matrices:                                                  //
    //-----------------------------------------------------------------------//
    double      T  = GetJCYsSinceEpoch(a_erm_epoch);
    auto [P, invP] = MkPrecMts(T);

    //-----------------------------------------------------------------------//
    // Nutation Angles and Matrices:                                         //
    //-----------------------------------------------------------------------//
    // Use DE440T data if available, otherwise use the analytical formulas:
    //
    auto [dPsi, dEps] =
      (DE440T::Bits::FromY <= a_erm_epoch && a_erm_epoch < DE440T::Bits::ToY)
      ? GetNutAnglesDE440T(a_erm_epoch)
      : GetNutAnglesAnalyt(T);

    //-----------------------------------------------------------------------//
    // Nutation Matrices:                                                    //
    //-----------------------------------------------------------------------//
    auto [N, invN, cosEps] = MkNutMts(T, dPsi, dEps);

    //-----------------------------------------------------------------------//
    // Over-All "PN" and "invPN" Matrices:                                   //
    //-----------------------------------------------------------------------//
    // XXX: Bias Correction and Polar Motion Matrices are considered to be too
    // small, and not used. Thus:
    MtxMult33(P,    N,    &m_PN);
    MtxMult33(invN, invP, &m_invPN);

    // Check:
    DEBUG_ONLY
    (
      Mtx33 Tmp;
      MtxMult33(m_PN,  m_invPN, &Tmp);
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        assert(DimTypes::Bits::CEMaths::ApproxEqual
                (Tmp(i,j), (i==j) ? 1.0 : 0.0));
    )
    //-----------------------------------------------------------------------//
    // An Approximate Equation of the Origins Ee = ERA - GAST:               //
    //-----------------------------------------------------------------------//
    m_Ee     = dPsi  *  cosEps / 15.0;
    m_DeltaT = GetDeltaT(a_erm_epoch);

    // All Done!
  }

  //=========================================================================//
  // As above, but the ERM Epoch is given by TT:                             //
  //=========================================================================//
  EarthRotationModel::EarthRotationModel(TT a_erm_epoch)
  : EarthRotationModel
    (TT::EpochJYr + To_Time_jyr(a_erm_epoch.GetTimeSinceEpoch()))
  {}

  //=========================================================================//
  // "EarthRotationModel" Non-Default "constexpr" Ctor:                      //
  //=========================================================================//
  // Using the Analytical Nutations Model:
  //
  constexpr EarthRotationModel::EarthRotationModel(int a_erm_epoch_year)
  {
    m_ermEpoch = Time_jyr(double(a_erm_epoch_year));

    //-----------------------------------------------------------------------//
    // Precession Matrices:                                                  //
    //-----------------------------------------------------------------------//
    double      T  = GetJCYsSinceEpoch(m_ermEpoch);
    auto [P, invP] = MkPrecMts(T);

    //-----------------------------------------------------------------------//
    // Nutation Angles and Matrices:                                         //
    //-----------------------------------------------------------------------//
    // Use analytic exprs for dPsi, dEps:
    auto         nutAngles = GetNutAnglesAnalyt(T);
    Angle dPsi = nutAngles.first;
    Angle dEps = nutAngles.second;

    auto [N, invN, cosEps] = MkNutMts(T, dPsi, dEps);

    //-----------------------------------------------------------------------//
    // Over-All "PN" and "invPN" Matrices:                                   //
    //-----------------------------------------------------------------------//
    // XXX: Bias Correction and Polar Motion Matrices are considered to be too
    // small, and not used. Thus:
    MtxMult33(P,    N,    &m_PN);
    MtxMult33(invN, invP, &m_invPN);

    // Check:
    DEBUG_ONLY
    (
      Mtx33 Tmp;
      MtxMult33(m_PN,  m_invPN, &Tmp);
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        assert(DimTypes::Bits::CEMaths::ApproxEqual
                (Tmp(i,j), (i==j) ? 1.0 : 0.0));
    )
    //-----------------------------------------------------------------------//
    // An Approximate Equation of the Origins Ee = ERA - GAST:               //
    //-----------------------------------------------------------------------//
    m_Ee     = dPsi  * cosEps / 15.0;
    m_DeltaT = GetDeltaT(m_ermEpoch);

    // All Done!
  }
}
// End namespace SpaceBallistics

// vim:ts=2:et
//===========================================================================//
//                       "Src/Missions/RTLS1-OptGen.cpp":                    //
//           Return-to-(Launch/Landing)-Site: Parametric Optimisation        //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <vector>

namespace SpaceBallistics
{
//===========================================================================//
// "SetCtlParams":                                                           //
//===========================================================================//
// Physical Params:
// [PropMassS,       CoastDur,  BBBurnDur,      EntryBurnQ,    EntryBurnDur,
//  EntryBurnThrtAL, LandBurnH, LandBurnThrt1L, BBBurnSinTheta[1..4]]
//
// However, "a_opt_params_n" are all normalised to 0..1:
//---------------------------------------------------------------------------//
// Individual Args Form:                                                     //
//---------------------------------------------------------------------------//
// NB: "PropMassS" cannot be set here; it can  only be set when the "RTLS1" obj
// is constructed!
//
void RTLS1::SetCtlParams
(
  double a_coastDurN,     double a_bbBurnDurN,
  double a_entryBurnQN,   double a_entryBurnDurN,
  double a_entryBurnThrtAN,
  double a_landBurnHN,    double a_landBurnThrt1N,
  double a_sinTheta0,     double a_sinTheta1,
  double a_sinTheta2,     double a_sinTheta3
)
{
  // NB: "m_coastDur" is set relative to the "flight time to the ballistic top",
  // assuming a const gravity acceleration (as @ "m_hS"),  since it would prob-
  // ably be too late and sub-optimal to start the BBBurn during the descending
  // branch of the ballistic trajectory:
  AccK gS    = K / Sqr(R + m_hS);
  VelK VrS   = m_VS * Cos(m_psiS);
  assert(IsPos(VrS));
  Time toTop = VrS / gS;

  assert(0.0       <= a_coastDurN       && a_coastDurN       <= 1.0);
  m_coastDur        = a_coastDurN       *  toTop;

  assert(0.0       <= a_bbBurnDurN      && a_bbBurnDurN      <= 1.0);
  m_bbBurnDur       = a_bbBurnDurN      *  MaxBBBurnDur;

  assert(0.0       <= a_entryBurnQN     && a_entryBurnQN     <= 1.0);
  m_entryBurnQ      = a_entryBurnQN     *  MaxEntryBurnQ;

  assert(0.0       <= a_entryBurnDurN   && a_entryBurnDurN   <= 1.0);
  m_entryBurnDur    = a_entryBurnDurN   *  MaxEntryBurnDur;

  // NB: Take the actual Stage1 throttling range and the number of
  // engines into account:
  assert(0.0       <= a_entryBurnThrtAN && a_entryBurnThrtAN <= 1.0);
  m_entryBurnThrtAL =
    Base::m_minThrtL1 / double(NE) * (1.0 - a_entryBurnThrtAN) +
    a_entryBurnThrtAN;

  assert(0.0       <= a_landBurnHN      && a_landBurnHN      <= 1.0);
  m_landBurnH       = a_landBurnHN      *  MaxLandBurnH;

  assert(0.0       <= a_landBurnThrt1N && a_landBurnThrt1N   <= 1.0);
  // NB: Take the actual Stage1 throttling range into account:
  m_landBurnThrt1L  =
    Base::m_minThrtL1 * (1.0 - a_landBurnThrt1N) + a_landBurnThrt1N;

  // sin(theta) coeffs:
  assert(0.0       <= a_sinTheta0       && a_sinTheta0       <= 1.0);
  m_bbBurnSinTheta[0] = 2.0 * a_sinTheta0 - 1.0;

  assert(0.0       <= a_sinTheta1       && a_sinTheta1       <= 1.0);
  m_bbBurnSinTheta[1] = 2.0 * a_sinTheta1 - 1.0;

  assert(0.0       <= a_sinTheta2       && a_sinTheta2       <= 1.0);
  m_bbBurnSinTheta[2] = 2.0 * a_sinTheta2 - 1.0;

  assert(0.0       <= a_sinTheta3       && a_sinTheta3       <= 1.0);
  m_bbBurnSinTheta[3] = 2.0 * a_sinTheta3 - 1.0;
  // All Done!
}

//---------------------------------------------------------------------------//
// Vector Form:                                                              //
//---------------------------------------------------------------------------//
void RTLS1::SetCtlParams(std::vector<double> const& a_opt_params_n)
{
  int    np = int(a_opt_params_n.size());
  assert(NP - 3 <= np && np <= NP);

  // NB: As above: Skipping  a_opt_params[0]  here, since "PropMassS" is not a
  // Ctl Param -- it can only be set when the "RTLS1" obj is constructed:
  SetCtlParams
  (
    a_opt_params_n[1],                      // coastDurN
    a_opt_params_n[2],                      // bbBurnDurN
    a_opt_params_n[3],                      // entryBurnQN
    a_opt_params_n[4],                      // entryBurnDurN
    a_opt_params_n[5],                      // entryBurnThrtAL
    a_opt_params_n[6],                      // landBurnHN
    a_opt_params_n[7],                      // landBurnThrtN

    // sin(theta) coeffs: NB: the default normalised coeffs are equal to 0.5,
    // which will translate into 0 (see the above function):
    a_opt_params_n[8],                      // sinTheta0
    (np >= 10) ? a_opt_params_n[ 9] : 0.5,  // sinTheta1
    (np >= 11) ? a_opt_params_n[10] : 0.5,  // sinTheta2
    (np == 12) ? a_opt_params_n[11] : 0.5   // sinTheta3
  );
}

//===========================================================================//
// "OptRes" Non-Default Ctor:                                                //
//===========================================================================//
RTLS1::OptRes::OptRes(RTLS1 const& a_rtls)
: m_propMassS     (a_rtls.Base::m_propMass1),
  m_coastDur      (a_rtls.m_coastDur),
  m_bbBurnDur     (a_rtls.m_bbBurnDur),
  m_bbBurnSinTheta{a_rtls.m_bbBurnSinTheta[0],
                   a_rtls.m_bbBurnSinTheta[1],
                   a_rtls.m_bbBurnSinTheta[2],
                   a_rtls.m_bbBurnSinTheta[3]},
  m_entryBurnQ    (a_rtls.m_entryBurnQ),
  m_entryBurnDur  (a_rtls.m_entryBurnDur),
  m_landBurnH     (a_rtls.m_landBurnH),
  m_landBurnThrt1L(a_rtls.m_landBurnThrt1L)
{
  static_assert(RTLS1::NS == 4);
}

//===========================================================================//
// "OptRes" Output:                                                          //
//===========================================================================//
std::ostream& operator<< (std::ostream& a_os, RTLS1::OptRes const& a_res)
{
  a_os
    << "\tpropMassS      = "   << a_res.m_propMassS
    << "\n\tcoastDur       = " << a_res.m_coastDur
    << "\n\tbbBurnDur      = " << a_res.m_bbBurnDur;

  for (int i = 0; i < RTLS1::NS; ++i)
    a_os << "\n\tsinTheta["   << i << "]   = " << a_res.m_bbBurnSinTheta[i];

  a_os
    << "\n\tentryBurnQ     = " << a_res.m_entryBurnQ
    << "\n\tentryBurnDur   = " << a_res.m_entryBurnDur
    << "\n\tlandBurnH      = " << a_res.m_landBurnH
    << "\n\tlandBurnThrt1L = " << a_res.m_landBurnThrt1L
    << std::endl;

  return a_os;
}

//===========================================================================//
// "FindOptimalAscentCtls"                                                   //
//===========================================================================//
std::pair<std::optional<RTLS1::OptRes>,
          std::optional<RTLS1::Base::RunRes>>  // If the FinalRun is performed
RTLS1::FindOptimalReturnCtls
(
  std::string const& a_config_ini,
  std::ostream*      a_os       // May be NULL
)
{
  //-------------------------------------------------------------------------//
  // Open and Parse the Config.ini File:                                     //
  //-------------------------------------------------------------------------//
  boost::property_tree::ptree  pt;
  boost::property_tree::ini_parser::read_ini(a_config_ini, pt);

  // "Fixed" params for the "prototype" "RTLS1" obj:
  Mass   maxFullMass1           (pt.get<double>("LV.MaxStartMass1"));
  double fullK1           =      pt.get<double>("LV.K1");
  double fullPropRem1     =      pt.get<double>("LV.PropRem1");
  Time   IspSL1                 (pt.get<double>("LV.IspSL1" ));
  Time   IspVac1                (pt.get<double>("LV.IspVac1"));
  ForceK thrustVacI1      = Mass(pt.get<double>("LV.ThrustVacI1")) * g0K;
  double minThrtL1        =      pt.get<double>("LV.MinThrtL1");
  Len    diam                   (pt.get<double>("LV.Diameter"));

  // Return-to-Launch/Landing-Site Params: Considered to be "fixed" here, but
  // may eventually be made variable in the "outer loop".   They describe the
  // conds @ Stage1 Separation.
  // NB: (HS, LS) are the Altitude and Down-Range Distance (of Stage1 @ Separ-
  // ation Time) from the Landing Site, NOT from the Launch Site (though they
  // are assumed to be close from to each other):
  LenK   hS                     (pt.get<double>("RTLS.hS"));
  LenK   lS                     (pt.get<double>("RTLS.lS"));
  VelK   VS                     (pt.get<double>("RTLS.VS"));
  Angle  psiS                   (pt.get<double>("RTLS.psiS"));

  // "Technical" Params:
  Time   odeIntegrStep          (pt.get<double>("Technical.ODEIntegrStep"));
  bool   withFinalRun     =      pt.get<bool>  ("Technical.WithFinalRun");
  int    finalRunLogLevel =      pt.get<int>   ("Technical.FinalRunLogLevel");
  int    optMaxEvals      =      pt.get<int>   ("Technical.OptMaxEvals");
  int    optLogLevel      =      pt.get<int>   ("Technical.OptLogLevel");

  // NOMAD-Specific Params:
  int    optSeed          =      pt.get<int>   ("Technical.NOMADSeed");
  bool   stopIfFeasible   =      pt.get<bool>  ("Technical.NOMADStopIfFeasible");
  double useVNS           =      pt.get<double>("Technical.NOMADUseVNS");
  if (useVNS < 0.0 || useVNS >= 1.0)      // 0: VNS not used
    throw std::invalid_argument("NOMADUseVNS: Must be in [0..1)");
  bool   useMT            =      pt.get<bool>  ("Technical.NOMADUseMT");

  // XXX: The optimisation params are currently NOT configurable, with the
  // exception of the degree of "sin(theta)" expansion (0..3):
  int    sinThetaDeg      =      pt.get<int>   ("Opt.sinThetaDeg");
  if (sinThetaDeg < 0     || sinThetaDeg > 3)
    throw std::invalid_argument("sinThetaDeg: Must be in [0..3]");

  double minRelPropMassS  =      pt.get<double>("Opt.minRelPropMassS");
  if (minRelPropMassS  <= 0 || minRelPropMassS  > 1)
    throw std::invalid_argument("minRelPropMassS: Must be in (0..1]");

  double minRelBBBurnDur  =      pt.get<double>("Opt.minRelBBBurnDur");
  if (minRelBBBurnDur <= 0 || minRelBBBurnDur > 1)
    throw std::invalid_argument("minRelBBBurnDur: Must be in (0..1]");

  // The initial vals of all (NORMALISED) params:  0.5:
  std::vector<double> initParamsN(size_t(NP - NS + 1 + sinThetaDeg), 0.5);

  // The Limits:
  // TODO: Possibly add the Acceleration Limit at Landing as well:
  LenK     landDLLimit(pt.get<double>("Opt.LandDLLimit"));
  VelK     landVLimit (pt.get<double>("Opt.LandVLimit")) ;
  Pressure QLimit     (pt.get<double>("Opt.QLimit"  ));

  //-------------------------------------------------------------------------//
  // Create the "prototype" "RTLS1" obj:                                     //
  //-------------------------------------------------------------------------//
  RTLS1 proto
  (
    maxFullMass1,        fullK1,  fullPropRem1, IspSL1, IspVac1, thrustVacI1,
    minThrtL1,           diam,
    RTLS1::MaxPropMassS, hS,      lS,       VS, psiS,
    odeIntegrStep,       a_os,    optLogLevel
  );

  proto.SetCtlParams(initParamsN);

  //-------------------------------------------------------------------------//
  // Actually run the Optimiser:                                             //
  //-------------------------------------------------------------------------//
  // XXX: For the moment, only the NOMAD Optimiser is available.
  // The results are returned via the "initParamsN":
  bool ok =
    RunNOMAD
    (
      &proto,          maxFullMass1,    fullK1,         fullPropRem1, diam,
      &initParamsN,
      minRelPropMassS, minRelBBBurnDur, landDLLimit,    landVLimit,   QLimit,
      optMaxEvals,     optSeed,         stopIfFeasible, useVNS,       useMT
    );
  if (!ok)
      return std::make_pair(std::nullopt, std::nullopt);

  //-------------------------------------------------------------------------//
  // Post-Processing:                                                        //
  //-------------------------------------------------------------------------//
  // In order to get the final "OptRes",  put  the "initVals" (now containing
  // the opt args found) into the "proto" obj, and extract the "optRes"  from
  // the "proto":
  //
  proto.SetCtlParams(initParamsN);

  // The Optimisation Result:
  std::optional<OptRes> optRes = OptRes(proto);

  if (withFinalRun)
  {
    //-----------------------------------------------------------------------//
    // Perform the Final Run on the optimised params:                        //
    //-----------------------------------------------------------------------//
    // "proto" is already set up with those params, except for the LogLevel:
    //
    proto.m_logLevel    = finalRunLogLevel;
    Base::RunRes runRes = proto.Run();
    return std::make_pair(optRes, runRes);
  }
  // Otherwise, there is no "runRes":
  return std::make_pair(optRes, std::nullopt);
}
}
// End namespace SpaceBallistics

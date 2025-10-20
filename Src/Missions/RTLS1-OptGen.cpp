// vim:ts=2:et
//===========================================================================//
//                       "Src/Missions/RTLS1-Optim.cpp":                     //
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
// Dimensioned (Physical) Params:
// [PropMassS, CoastDur,      BBBurnDur, EntryBurnQ, EntryBurnDur,
//  LandBurnH, LandBurnThrtL, BBBurnSinTheta[1..4]]
//
// However, "a_opt_params_n" are all normalised to 0..1:
//
void RTLS1::SetCtlParams(std::vector<double> const& a_opt_params_n)
{
  int    np = int(a_opt_params_n.size());
  assert(NP - 3 <= np && np <= NP);

  // XXX: Scaling vals are currently consts, and are set here.
  constexpr Mass     MaxPropMassS     = 30000.0_kg;
  constexpr Time     MaxCoastDur      =   300.0_sec;
  constexpr Time     MaxBBBurnDur     =    20.0_sec;  // Actually too much!
  constexpr Pressure MaxEntryBurnQ     (30000.0);
  constexpr Time     MaxEntryBurnDur  =    10.0_sec;  // Actually too much!
  constexpr LenK     MaxLandBurnH     =     5.0_km;   // OK or too low?

  double propMassSN     = a_opt_params_n[0];
  assert(0.0     <= propMassSN    && propMassSN    <= 1.0);
  m_propMassS     = propMassSN    *  MaxPropMassS;

  double coastDurN      = a_opt_params_n[1];
  assert(0.0     <= coastDurN     && coastDurN     <= 1.0);
  m_coastDur      = coastDurN     *  MaxCoastDur;

  double bbBurnDurN     = a_opt_params_n[2];
  assert(0.0     <= bbBurnDurN    && bbBurnDurN    <= 1.0);
  m_bbBurnDur     = bbBurnDurN    *  MaxBBBurnDur;

  double entryBurnQN    = a_opt_params_n[3];
  assert(0.0     <= entryBurnQN   && entryBurnQN   <= 1.0);
  m_entryBurnQ    = entryBurnQN   *  MaxEntryBurnQ;

  double entryBurnDurN  = a_opt_params_n[4];
  assert(0.0     <= entryBurnDurN && entryBurnDurN <= 1.0);
  m_entryBurnDur  = entryBurnDurN *  MaxEntryBurnDur;

  double landBurnHN     = a_opt_params_n[5];
  assert(0.0     <= landBurnHN    && landBurnHN    <= 1.0);
  m_landBurnH     = landBurnHN    *  MaxLandBurnH;

  double landBurnThrtN = a_opt_params_n[6];
  assert(0.0     <= landBurnThrtN && landBurnThrtN <= 1.0);
  // NB: Take the actual Stage1 throttling range into account:
  m_landBurnThrtL = Base::m_minThrtL1 * (1.0 - landBurnThrtN) + landBurnThrtN;

  for (int i = 0; i < NS; ++i)
  {
    double si   = (7 + i < np) ? a_opt_params_n[size_t(7 + i)] : 0.5;
    assert(0.0 <= si && si <= 1.0);
    // XXX: However, we assume that the actual coeffs in the "sin(theta)"
    // expansion are in the range  [-1..1]:
    m_bbBurnSinTheta[i] = 2.0 * si - 1.0;
  }
  // All Done!
}

//===========================================================================//
// "FindOptimalAscentCtls"                                                   //
//===========================================================================//
std::pair<std::optional<RTLS1::OptRes>,
          std::optional<RTLS1::RunRes>> // If the FinalRun is performed
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
  Mass   fullMass1              (pt.get<double>("LV.StartMass1"));
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
  Mass   propMassS              (pt.get<double>("RTLS.PropMassS"));
  LenK   hS                     (pt.get<double>("RTLS.HS"));
  LenK   lS                     (pt.get<double>("RTLS.LS"));
  VelK   VrS                    (pt.get<double>("RTLS.VrS"));
  VelK   VhorS                  (pt.get<double>("RTLS.VhorS"));

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

  // XXX: The optimisation params are currently NOT configurable, with the
  // exception of the degree of "sin(theta)" expansion (0..3):
  int    sinThetaDeg      =      pt.get<int>   ("Opt.sinThetaDeg");
  if (sinThetaDeg < 0 || sinThetaDeg > 3)
    throw std::invalid_argument("sinThetaDeg: Must be in [0..3]");

  // The initial vals of all (NORMALISED) params:  0.5:
  std::vector<double> initParamsN(size_t(8 + sinThetaDeg), 0.5);

  // The Limits:
  // TODO: Possibly add the Acceleration Limit at Landing as well:
  LenK     maxLandMissL(pt.get<double>("Opt.MaxLandMissL"));
  VelK     maxLandV    (pt.get<double>("Opt.MaxLandV")) ;
  Pressure QLimit      (pt.get<double>("Opt.QLimit"  ));

  //-------------------------------------------------------------------------//
  // Create the "prototype" "RTLS1" obj:                                     //
  //-------------------------------------------------------------------------//
  RTLS1 proto
  (
    fullMass1, fullK1, fullPropRem1, IspSL1, IspVac1, thrustVacI1,
    minThrtL1, diam,
    propMassS, hS,     lS,      VrS, VhorS,
    odeIntegrStep,     a_os,    optLogLevel
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
      &proto,      &initParamsN, maxLandMissL,   maxLandV,  QLimit,
      optMaxEvals, optSeed,      stopIfFeasible, useVNS
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
    proto.m_logLevel = finalRunLogLevel;
    RunRes runRes    = proto.Run();
    return std::make_pair(optRes, runRes);
  }
  // Otherwise, there is no "runRes":
  return std::make_pair(optRes, std::nullopt);
}
}
// End namespace SpaceBallistics

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
// Dimensioned (Physical) Params:
// [PropMassS, CoastDur,      BBBurnDur, EntryBurnQ, EntryBurnDur,
//  LandBurnH, LandBurnThrtL, BBBurnSinTheta[1..4]]
//
// However, "a_opt_params_n" are all normalised to 0..1:
//---------------------------------------------------------------------------//
// Individual Args Form:                                                     //
//---------------------------------------------------------------------------//
void RTLS1::SetCtlParams
(
  double a_propMassSN,    double a_coastDurN,
  double a_bbBurnDurN,    double a_entryBurnQN,
  double a_entryBurnDurN, double a_landBurnHN,
  double a_landBurnThrtN,
  double a_sinTheta0,     double a_sinTheta1,
  double a_sinTheta2,     double a_sinTheta3
)
{
  assert(0.0     <= a_propMassSN    && a_propMassSN    <= 1.0);
  m_propMassS     = a_propMassSN    *  MaxPropMassS;

  assert(0.0     <= a_coastDurN     && a_coastDurN     <= 1.0);
  m_coastDur      = a_coastDurN     *  MaxCoastDur;

  assert(0.0     <= a_bbBurnDurN    && a_bbBurnDurN    <= 1.0);
  m_bbBurnDur     = a_bbBurnDurN    *  MaxBBBurnDur;

  assert(0.0     <= a_entryBurnQN   && a_entryBurnQN   <= 1.0);
  m_entryBurnQ    = a_entryBurnQN   *  MaxEntryBurnQ;

  assert(0.0     <= a_entryBurnDurN && a_entryBurnDurN <= 1.0);
  m_entryBurnDur  = entryBurnDurN *  MaxEntryBurnDur;

  assert(0.0     <= a_landBurnHN    && a_landBurnHN    <= 1.0);
  m_landBurnH     = a_landBurnHN    *  MaxLandBurnH;

  assert(0.0     <= a_landBurnThrtN && a_landBurnThrtN <= 1.0);
  // NB: Take the actual Stage1 throttling range into account:
  m_landBurnThrtL =
    Base::m_minThrtL1 * (1.0 - a_landBurnThrtN) + a_landBurnThrtN;

  // sin(theta) coeffs:
  assert(0.0     <= a_sinTheta0     && a_sinTheta0    <= 1.0);
  m_bbBurnSinTheta[0] = 2.0 * a_sinTheta0 - 1.0;

  assert(0.0     <= a_sinTheta1     && a_sinTheta1    <= 1.0);
  m_bbBurnSinTheta[1] = 2.0 * a_sinTheta1 - 1.0;

  assert(0.0     <= a_sinTheta2     && a_sinTheta2    <= 1.0);
  m_bbBurnSinTheta[2] = 2.0 * a_sinTheta2 - 1.0;

  assert(0.0     <= a_sinTheta3     && a_sinTheta3    <= 1.0);
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

  SetCtlParams
  (
    a_opt_params_n[0],                      // propMassSN
    a_opt_params_n[1],                      // coastDurN
    a_opt_params_n[2],                      // bbBurnDurN
    a_opt_params_n[3],                      // entryBurnQN
    a_opt_params_n[4],                      // entryBurnDurN
    a_opt_params_n[5],                      // landBurnHN
    a_opt_params_n[6],                      // landBurnThrtN

    // sin(theta) coeffs: NB: the default normalised coeffs are equal to 0.5,
    // which will translate into 0 (see the avove function):
    a_opt_params_n[7],                      // sinTheta0
    (np >=  9) ? a_opt_params_n[ 8] : 0.5,  // sinTheta1
    (np >= 10) ? a_opt_params_n[ 9] : 0.5,  // sinTheta2
    (np == 11) ? a_opt_params_n[10] : 0.5   // sinTheta3
  );
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
    proto.m_logLevel    = finalRunLogLevel;
    Base::RunRes runRes = proto.Run();
    return std::make_pair(optRes, runRes);
  }
  // Otherwise, there is no "runRes":
  return std::make_pair(optRes, std::nullopt);
}
}
// End namespace SpaceBallistics

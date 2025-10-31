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
// "RTLS1" Ctor from a Proto and Normalised Ctl Params:                      //
//===========================================================================//
// Params:
// [
//  PropMassS,  CoastDur,
//  BBBurnDur,  BBBurnThrtL0,   BBBurnThrtL1,    BBBurnTheta0, BBBurnTheta1,
//  EntryBurnQ, EntryBurnDur,   EntryBurnThrtL0, EntryBurnThrtL1,
//  LandBurnH,  LandBurnThrtL0, LandBurnThrtL1
// ]:
RTLS1::RTLS1
(
  RTLS1  const&               a_proto,
  Pressure                    a_Q_limit,
  double a_propMassSN,        double a_coastDurN,
  double a_bbBurnDurN,        double a_bbBurnThrtL0N,
  double a_bbBurnThrtL1N,     double a_bbBurnTheta0N,
  double a_bbBurnTheta1N,
  double a_entryBurnQN,       double a_entryBurnDurN,
  double a_entryBurnThrtL0N,  double a_entryBurnThrtL1N,
  double a_landBurnHN,        double a_landBurnThrtL0N,
  double a_landBurnThrtL1N
)
: //-------------------------------------------------------------------------//
  // Construct the "RTLS1" Obj from "Proto":                                 //
  //-------------------------------------------------------------------------//
  RTLS1
	(
    // Stage Params: NB: Here we use the stored "m_fpl*" params:
    a_proto.m_fplMass1,
    a_proto.m_fplK1,
    a_proto.m_fplPropRem1,
    a_proto.Base::m_IspSL1,
    a_proto.Base::m_IspVac1,
    a_proto.Base::m_thrustVacI1,
    a_proto.Base::m_minThrtL1,
    a_proto.m_diam,
		// IMPORTANT: Specifying the "PropMassS"!
    a_propMassSN * MaxPropMassS,
    // Mission Params:
    a_proto.m_hS,
    a_proto.m_lS,
    a_proto.m_VS,
    a_proto.m_psiS,
    // Integration and Output Params:
    a_proto.Base::m_odeIntegrStep,
    a_proto.Base::m_os,
    a_proto.Base::m_logLevel
	)
{
  assert(!IsNeg(a_Q_limit) && 0.0 < a_propMassSN && a_propMassSN <= 1.0);
  //-------------------------------------------------------------------------//
	// Now set the Ctl Params:                                                 //
	//-------------------------------------------------------------------------//
  // NB: "m_coastDur" is set relative to the "flight time to the ballistic top",
  // assuming a const gravity acceleration (as @ "m_hS"),  since it would prob-
  // ably be too late and sub-optimal to start the BBBurn during the descending
  // branch of the ballistic trajectory:
  AccK gS    = K / Sqr(R + m_hS);
  VelK VrS   = m_VS * Cos(m_psiS);
  assert(IsPos(VrS));
  Time toTop = VrS / gS;

  //------------//
  // Coast:     //
  //------------//
  assert(0.0       <= a_coastDurN        && a_coastDurN      <= 1.0);
  m_coastDur        = a_coastDurN        *  toTop;

  //------------//
  // BBBurn:    //
  //------------//
  assert(0.0       <= a_bbBurnDurN       && a_bbBurnDurN     <= 1.0);
  m_bbBurnDur       = a_bbBurnDurN       *  MaxBBBurnDur;

  assert(0.0       <= a_bbBurnThrtL0N    && a_bbBurnThrtL0N  <= 1.0);
  m_bbBurnThrtL0    =
    Base::m_minThrtL1 * (1.0 - a_bbBurnThrtL0N) +   a_bbBurnThrtL0N;
  assert(Base::m_minThrtL1  <= m_bbBurnThrtL0   &&  m_bbBurnThrtL0 <= 1.0);

  // NB: "m_bbBurnThrtL1" must be in [Base::m_minThrtL1 .. m_bbBurnThrtL0]:
  assert(0.0       <= a_bbBurnThrtL1N    && a_bbBurnThrtL1N  <= 1.0);
  m_bbBurnThrtL1    =
    Base::m_minThrtL1 * (1.0 - a_bbBurnThrtL1N) +
    m_bbBurnThrtL0    * a_bbBurnThrtL1N;
  assert(Base::m_minThrtL1  <= m_bbBurnThrtL1   &&
         m_bbBurnThrtL1     <= m_bbBurnThrtL0);

  // NB: "theta"s are centered on Pi:
  assert(0.0       <= a_bbBurnTheta0N    && a_bbBurnTheta0N  <= 1.0);
  m_bbBurnTheta0    = PI * (a_bbBurnTheta0N + 0.5);

  assert(0.0       <= a_bbBurnTheta1N    && a_bbBurnTheta1N  <= 1.0);
  m_bbBurnTheta1    = PI * (a_bbBurnTheta1N + 0.5);

  //------------//
  // EntryBurn: //
  //------------//
  assert(0.0       <= a_entryBurnQN      && a_entryBurnQN      <= 1.0);
  m_entryBurnQ      = a_entryBurnQN      *  a_Q_limit;

  assert(0.0       <= a_entryBurnDurN    && a_entryBurnDurN    <= 1.0);
  m_entryBurnDur    = a_entryBurnDurN    *  MaxEntryBurnDur;

  assert(0.0       <= a_entryBurnThrtL0N && a_entryBurnThrtL0N <= 1.0);
  m_entryBurnThrtL0 =
    Base::m_minThrtL1 * (1.0 - a_entryBurnThrtL0N) +  a_entryBurnThrtL0N;
  assert(Base::m_minThrtL1  <= m_entryBurnThrtL0   && m_entryBurnThrtL0 <= 1.0);

  // NB: "m_entryBurnThrtL1" must be in
  // [Base::m_minThrtL1 .. m_entryBurnThrtL0]:
  assert(0.0       <= a_entryBurnThrtL1N && a_entryBurnThrtL1N <= 1.0);
  m_entryBurnThrtL1 =
    Base::m_minThrtL1 * (1.0 - a_entryBurnThrtL1N) +
    m_entryBurnThrtL0 * a_entryBurnThrtL1N;
  assert(Base::m_minThrtL1  <= m_entryBurnThrtL1   &&
         m_entryBurnThrtL1  <= m_entryBurnThrtL0);

  //------------//
  // LandBurn:  //
  //------------//
  assert(0.0       <= a_landBurnHN      && a_landBurnHN      <= 1.0);
  m_landBurnH       = a_landBurnHN      *  MaxLandBurnH;

  assert(0.0       <= a_landBurnThrtL0N && a_landBurnThrtL0N <= 1.0);
  m_landBurnThrtL0  =
    Base::m_minThrtL1 * (1.0 - a_landBurnThrtL0N) + a_landBurnThrtL0N;
  assert(Base::m_minThrtL1  <= m_landBurnThrtL0  && m_landBurnThrtL0 <= 1.0);

  // NB: "m_landBurnThrtL1" must be in
  // [Base::m_minThrtL1 .. m_landBurnThrtL0]:
  assert(0.0       <= a_landBurnThrtL1N && a_landBurnThrtL1N <= 1.0);
  m_landBurnThrtL1 =
    Base::m_minThrtL1 * (1.0 - a_landBurnThrtL1N) +
    m_landBurnThrtL0  * a_landBurnThrtL1N;
  assert(Base::m_minThrtL1  <= m_landBurnThrtL1   &&
         m_landBurnThrtL1   <= m_landBurnThrtL0);

  // All Done!
}

//===========================================================================//
// "OptRes" Non-Default Ctor:                                                //
//===========================================================================//
RTLS1::OptRes::OptRes(RTLS1 const& a_rtls)
: m_propMassS      (a_rtls.Base::m_propMass1),
  m_coastDur       (a_rtls.m_coastDur),
  m_bbBurnDur      (a_rtls.m_bbBurnDur),
  m_bbBurnThrtL0   (a_rtls.m_bbBurnThrtL0),
  m_bbBurnThrtL1   (a_rtls.m_bbBurnThrtL1),
  m_bbBurnTheta0   (a_rtls.m_bbBurnTheta0),
  m_bbBurnTheta1   (a_rtls.m_bbBurnTheta1),
  m_entryBurnQ     (a_rtls.m_entryBurnQ),
  m_entryBurnDur   (a_rtls.m_entryBurnDur),
  m_entryBurnThrtL0(a_rtls.m_entryBurnThrtL0),
  m_entryBurnThrtL1(a_rtls.m_entryBurnThrtL1),
  m_landBurnH      (a_rtls.m_landBurnH),
  m_landBurnThrtL0 (a_rtls.m_landBurnThrtL0),
  m_landBurnThrtL1 (a_rtls.m_landBurnThrtL1)
{}

//===========================================================================//
// "OptRes" Output:                                                          //
//===========================================================================//
std::ostream& operator<< (std::ostream& a_os, RTLS1::OptRes const& a_res)
{
  return
    a_os
      <<   "\tpropMassS       = " << a_res.m_propMassS
      << "\n\tcoastDur        = " << a_res.m_coastDur
      << "\n\tbbBurnDur       = " << a_res.m_bbBurnDur
      << "\n\tbbBurnThrtL0    = " << a_res.m_bbBurnThrtL0
      << "\n\tbbBurnThrtL1    = " << a_res.m_bbBurnThrtL1
      << "\n\tbbBurnTheta0    = " << To_Angle_deg(a_res.m_bbBurnTheta0)
      << "\n\tbbBurnTheta1    = " << To_Angle_deg(a_res.m_bbBurnTheta1)
      << "\n\tentryBurnQ      = " << a_res.m_entryBurnQ
      << "\n\tentryBurnDur    = " << a_res.m_entryBurnDur
      << "\n\tentryBurnThrtL0 = " << a_res.m_entryBurnThrtL0
      << "\n\tentryBurnThrtL1 = " << a_res.m_entryBurnThrtL1
      << "\n\tlandBurnH       = " << a_res.m_landBurnH
      << "\n\tlandBurnThrtL0  = " << a_res.m_landBurnThrtL0
      << "\n\tlandBurnThrtL1  = " << a_res.m_landBurnThrtL1
      << std::endl;
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
  Mass   maxFPLMass1            (pt.get<double>("LV.MaxStartMass1"));
  double fplK1            =      pt.get<double>("LV.K1");
  double fplPropRem1      =      pt.get<double>("LV.PropRem1");
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
  std::vector<double> initParamsN(NP, 0.5);

  // The Limits:
  // TODO: Possibly add the Acceleration Limit at Landing as well:
  LenK     landDLLimit    (pt.get<double>("Opt.LandDLLimit"));
  VelK     landVLimit     (pt.get<double>("Opt.LandVLimit")) ;
  double   landAccGLimit = pt.get<double>("opt.LandAccGLimit");
  Pressure QLimit         (pt.get<double>("Opt.QLimit"  ));
  double   longGLimit    = pt.get<double>("Opt.LongGLimit");

  //-------------------------------------------------------------------------//
  // Create the "prototype" "RTLS1" obj:                                     //
  //-------------------------------------------------------------------------//
  RTLS1 proto
  (
    maxFPLMass1,         fplK1, fplPropRem1, IspSL1, IspVac1, thrustVacI1,
    minThrtL1,           diam,
    RTLS1::MaxPropMassS, hS,    lS,      VS, psiS,
    odeIntegrStep,       a_os,  optLogLevel
  );

  //-------------------------------------------------------------------------//
  // Actually run the Optimiser:                                             //
  //-------------------------------------------------------------------------//
  // XXX: For the moment, only the NOMAD Optimiser is available.
  // The results are returned via the "initParamsN":
  bool ok =
    RunNOMAD
    (
      &proto,          &initParamsN,
      minRelPropMassS, minRelBBBurnDur, landDLLimit,    landVLimit,
      landAccGLimit,   QLimit,          longGLimit,
      optMaxEvals,     optSeed,         stopIfFeasible, useVNS,  useMT
    );
  if (!ok)
      return std::make_pair(std::nullopt, std::nullopt);

  //-------------------------------------------------------------------------//
  // Post-Processing:                                                        //
  //-------------------------------------------------------------------------//
  // In order to get the final "OptRes",  put  the "initVals"  (now containing
  // the opt args found) into a new "RTLS1" obj, and extract the "optRes" from
  // it:
  RTLS1 rtls
  (
    proto,
    QLimit,
    initParamsN[ 0],  // propMassSN
    initParamsN[ 1],  // coastDurN
    initParamsN[ 2],  // bbBurnDurN
    initParamsN[ 3],  // bbBurnThrtL0N
    initParamsN[ 4],  // bbBurnThrtL1N
    initParamsN[ 5],  // bbBurnTheta0N
    initParamsN[ 6],  // bbBurnTheta1N
    initParamsN[ 7],  // entryBurnQN
    initParamsN[ 8],  // entryBurnDurN
    initParamsN[ 9],  // entryBurnThrtL0N,
    initParamsN[10],  // entryBurnThrtL1N,
    initParamsN[11],  // landBurnHN
    initParamsN[12],  // landBurnThrtL0N
    initParamsN[13]   // landBurnThrtL1N
  );
  static_assert(13 == NP-1);

  // The Optimisation Result:
  std::optional<OptRes> optRes = OptRes(rtls);

  if (withFinalRun)
  {
    //-----------------------------------------------------------------------//
    // Perform the Final Run on the optimised params (in "rtls"):            //
    //-----------------------------------------------------------------------//
    rtls.m_logLevel     = finalRunLogLevel;
    Base::RunRes runRes = rtls.Run();
    return std::make_pair(optRes, runRes);
  }
  // Otherwise, there is no "runRes":
  return std::make_pair(optRes, std::nullopt);
}
}
// End namespace SpaceBallistics

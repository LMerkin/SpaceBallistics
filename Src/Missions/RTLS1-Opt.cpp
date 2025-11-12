// vim:ts=2:et
//===========================================================================//
//                       "Src/Missions/RTLS1-OptGen.cpp":                    //
//           Return-to-(Launch/Landing)-Site: Parametric Optimisation        //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include "SpaceBallistics/Missions/MkNOMADParams.hpp"
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
//  EntryBurnQ, EntryBurnDur,   EntryBurnThrtL0, EntryBurnThrtL1
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
  double a_entryBurnThrtL0N,  double a_entryBurnThrtL1N
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
    // IMPORTANT: Specifying the "PropMassS" here, relative to the (previously
    // estimated) value stored in "a_proto", in the (1 +- PropMassSRelVar) ra-
    // nge, AND making sure it is not below the UnSpendable Remnant:
    std::max
    (
      a_proto.Base::m_propMass1 *
        (1.0 - PropMassSRelVar + 2.0 * a_propMassSN * PropMassSRelVar),
      a_proto.m_unSpendable1
    ),
    // Mission Params:
    a_proto.m_hS,
    a_proto.m_lS,
    a_proto.m_VS,
    a_proto.m_psiS,
    a_proto.m_dVhor,
    // Integration and Output Params:
    a_proto.Base::m_odeIntegrStep,
    a_proto.Base::m_os,
    a_proto.Base::m_logLevel
  )
{
  assert(!IsNeg(a_Q_limit) && 0.0 <= a_propMassSN && a_propMassSN <= 1.0);
  Mass propMassS = Base::m_propMass1;

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
  // Throttling Levels:
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

  // "Theta" Levels:
  // IMPORTANT: "theta"s are assumed to be in [min*Pi, Pi] for the BNBurn.
  // (*) there is no point in using theta > Pi: in that case, we would give the
  //     stage a down-ward impulse which is clearly counter-productive;
  // (*) theta < Pi is sometimes useful (a small upward impulse MAY be helpful),
  //     but not too much: if theta < 3/4*Pi, then the vertical impulse is high-
  //     er than the horizontal one (which provides the actual Boost-Back);
  //     that's why we limit "theta"s from below by min*Pi, where min > 3/4:
  // So: theta{0,1} in [MinBBBurnThetaPiFrac .. 1] * Pi:
  //
  assert(0.0       <= a_bbBurnTheta0N    && a_bbBurnTheta0N  <= 1.0);
  m_bbBurnTheta0    =
    PI * (MinBBBurnThetaPiFrac * (1.0 - a_bbBurnTheta0N) + a_bbBurnTheta0N);

  assert(0.0       <= a_bbBurnTheta1N    && a_bbBurnTheta1N  <= 1.0);
  m_bbBurnTheta1    =
    PI * (MinBBBurnThetaPiFrac * (1.0 - a_bbBurnTheta1N) + a_bbBurnTheta1N);

  double absCosTheta0   = Abs(Cos(m_bbBurnTheta0));
  double absCosTheta1   = Abs(Cos(m_bbBurnTheta1));
  double minAbsCosTheta = std::min(absCosTheta0, absCosTheta1);
  double maxAbsCosTheta = std::max(absCosTheta0, absCosTheta1);

  // Duration:
  // Estimate the "mid-level" "BBurnDur", compute the range of vals and the act-
  // ual val in that range.
  // The range of possible "BBurnDur"s depends on the "dVhor" velocity component
  // via the integral of motion in the horizontal direction, with variable mass.
  // It also depends on the assumed values of the Horizontal Thrust:
  // (*) Only "BBBurnEngPart" of the engines are burning;
  // (*) they may be further  throttled;
  // (*) the thrust vector may not be fully horizontal:
  //
  VelK   W       = m_IspVac1 * g0K;
  VelK   WhorMax = W * maxAbsCosTheta;
  VelK   WhorMin = W * minAbsCosTheta;
  assert(IsPos(WhorMin) && WhorMin <= WhorMax);

  MassRate muMax = m_thrustVacI1 * BBBurnEngPart * m_bbBurnThrtL0 / W;
  MassRate muMin = m_thrustVacI1 * BBBurnEngPart * m_bbBurnThrtL1 / W;
  assert(IsPos(muMin) &&  muMin <= muMax);
  MassRate muAvg = 0.5 * (muMin +  muMax);

  Time   tBBBMin = (propMassS / muMax) * (1.0 - Exp(- m_dVhor / WhorMax));
  Time   tBBBMax = (propMassS / muMin) * (1.0 - Exp(- m_dVhor / WhorMin));
  assert(tBBBMin <= tBBBMax);

  // Make the range wider:
  tBBBMin /= (1.0 + BBBurnRelDurVar);
  tBBBMax *= (1.0 + BBBurnRelDurVar);
  assert(tBBBMin <  tBBBMax);

  // Constrain the "tBBBMax" (and even the "tBBBMin") by the spendable
  // propellant (@ the "muAvg" rate). XXX: We do not make any reserves here for
  // other burns (EntryBurn and LandBurn):
  Time   tConstr =  (propMassS - Base::m_unSpendable1) / muAvg;
  assert(IsPos(tConstr));

  tBBBMin = std::min(tBBBMin, tConstr);
  tBBBMax = std::min(tBBBMax, tConstr);
  assert(tBBBMin <= tBBBMax);

  // And now, the actual "m_bbBurnDur" will be selected in that range (XXX which
  // is not static -- it depends on other params!):
  assert      (0.0 <= a_bbBurnDurN && a_bbBurnDurN <= 1.0);
  m_bbBurnDur =
    tBBBMin * (1.0  - a_bbBurnDurN) + tBBBMax * a_bbBurnDurN;

  if (a_proto.m_os != nullptr && a_proto.m_logLevel >= 3)
    *a_proto.m_os
      << "# PropMassS="    << propMassS  .Magnitude()
      << " kg, bbBurnDur=" << m_bbBurnDur.Magnitude()
      << " in ["           << tBBBMin    .Magnitude()
      << " .. "            << tBBBMax    .Magnitude()
      << "] sec"           << std::endl;

  //------------//
  // EntryBurn: //
  //------------//
  assert(0.0       <= a_entryBurnQN      && a_entryBurnQN      <= 1.0);
  m_entryBurnQ      = a_entryBurnQN      * a_Q_limit;

  assert(0.0       <= a_entryBurnDurN    && a_entryBurnDurN    <= 1.0);
  m_entryBurnDur    = a_entryBurnDurN    * MaxEntryBurnDur;

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

  //-----------//
  // LandBurn: //
  //-----------//
  // Don't set the LandBurn params yet -- they are set dynamically or externally
  // at a later stage; for now, they are NANs:
  m_landBurnH     = LenK(NAN);
  m_landBurnGamma = NAN;

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
  m_landBurnGamma  (a_rtls.m_landBurnGamma)
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
      << "\n\tlandBurnGamma   = " << a_res.m_landBurnGamma
      << std::endl;
}

//===========================================================================//
// "RTLS1::NOMADMainEvaluator": Helper Class used in NOMAD Optimisation:     //
//===========================================================================//
class RTLS1::NOMADMainEvaluator final: public NOMAD::Evaluator
{
private:
  //=========================================================================//
  // Data Flds:                                                              //
  //=========================================================================//
  // LV Params:
  RTLS1 const*  const m_proto;

  // Optimisation Limits (Constraints):
  LenK          const m_landDLLimit;
  Pressure      const m_QLimit;
  double        const m_longGLimit;

public:
  //=========================================================================//
  // Non-Default Ctor, Dtor:                                                 //
  //=========================================================================//
  NOMADMainEvaluator
  (
    std::shared_ptr<NOMAD::EvalParameters> const& a_params,
    // LV Proto:
    RTLS1 const*                                  a_proto,
    // Optimisation Limits (Constraints):
    LenK                                          a_land_dL_limit,
    Pressure                                      a_Q_limit,
    double                                        a_longG_limit
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_landDLLimit   (a_land_dL_limit),
    m_QLimit        (a_Q_limit),
    m_longGLimit    (a_longG_limit)
  {
    assert(m_proto != nullptr);

    if (!(IsPos(m_landDLLimit) && IsPos(m_QLimit) && IsPos(m_longGLimit)))
      throw std::invalid_argument
            ("RTLS1::NOMADMainEvaluator::Ctor: Invalid Limit(s)");
  }

  ~NOMADMainEvaluator() override {}

  //=========================================================================//
  // "eval_x": The Actual Evaluation Method for NOMAD:                       //
  //=========================================================================//
  bool eval_x
  (
    NOMAD::EvalPoint&    a_x,  // Not "const" -- the result is also set here
    NOMAD::Double const& UNUSED_PARAM(a_hMax),
    bool&                a_countEval
  )
  const override
  {
    // XXX: For the moment, all params are optimised over:
    assert(a_x.size() == NP);

    //-----------------------------------------------------------------------//
    // For Thread-Safety, construct a new "RTLS1" obj from "Proto":          //
    //-----------------------------------------------------------------------//
    RTLS1 rtls
    (
      *m_proto,
      m_QLimit,
      a_x[ 0].todouble(),    // propMassSN
      a_x[ 1].todouble(),    // coastDurN
      a_x[ 2].todouble(),    // bbBurnDurN
      a_x[ 3].todouble(),    // bbBurnThrtLN0
      a_x[ 4].todouble(),    // bbBurnThrtLN1
      a_x[ 5].todouble(),    // bbBurnTheta0
      a_x[ 6].todouble(),    // bbBurnTheta1
      a_x[ 7].todouble(),    // entryBurnQN
      a_x[ 8].todouble(),    // entryBurnDurN
      a_x[ 9].todouble(),    // entryBurnThrtLN0
      a_x[10].todouble()     // entryBurnThrtLN1
    );
    static_assert(10 == NP-1);

    //-----------------------------------------------------------------------//
    // Run the Integrator!                                                   //
    //-----------------------------------------------------------------------//
    // NB: Exceptions are handled inside "Run":
    RunRes res = rtls.Run(RTLS1::FlightMode::Coast);

    // IMPORTANT: NOMAD allows us to indicate that evaluation has failed:
    if (res.m_rc == RunRC::Error)
      return false;

    //-----------------------------------------------------------------------//
    // Push the results back to NOMAD in string form:                        //
    //-----------------------------------------------------------------------//
    char  buff[512];
    char* curr = buff;

    // The Objective Function Value (to me minimised): It is "propMassS" itself
    // (provided that all constraints are satisfied):
    curr += sprintf(curr, "%.16e",  rtls.m_propMass1.Magnitude());

    // Constraints:
    // Down-Range Miss (the aiming point is 0):
    LenK dL = Abs(res.m_LT);
    curr += sprintf(curr, " %.16e", (dL         - m_landDLLimit).Magnitude());

    // The MaxQ encountered:
    assert(!IsNeg(res.m_maxQ));
    curr += sprintf(curr, " %.16e", (res.m_maxQ - m_QLimit).Magnitude());

    // The Max LongG encountered:
    assert(!IsNeg(res.m_maxLongG));
    curr += sprintf(curr, " %.16e",  res.m_maxLongG - m_longGLimit);

    // Done!
    assert(size_t(curr - buff) <= sizeof(buff));

    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 2)
#   pragma omp critical(Output)
    {
      (*m_proto->m_os)   << "# PARAMS:";
      for (int i = 0; i < int(a_x.size()); ++i)
        (*m_proto->m_os) << "  " << a_x[size_t(i)].todouble();

      *(m_proto->m_os)   << "\n# RES   :  " << buff << std::endl << std::endl;
    }

    // Set the results back in "a_x":
    a_x.setBBO(buff);
    a_countEval = true;
    return        true;
  }
};

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

  // The initial vals of all (NORMALISED) params: 0.5:
  std::vector<double> initParamsN(NP, 0.5);

  // The Limits: XXX: Landing Velocity and Acceleration Limits are treated
  // SEPARETELY:
  LenK     landDLLimit    (pt.get<double>("Opt.LandDLLimit"));
  Pressure QLimit         (pt.get<double>("Opt.QLimit"  ));
  double   longGLimit    = pt.get<double>("Opt.LongGLimit");

  //-------------------------------------------------------------------------//
  // Create the "prototype" "RTLS1" obj:                                     //
  //-------------------------------------------------------------------------//
  // NB: "propMassS" and "dVhor" are estimates only (to facilitate the subsequ-
  // ent optimisation process), which are pre-computed and memoised in the
  // "proto":
  Mass me1 = maxFPLMass1 *  (1.0  - fplK1);         // Empty Mass
  Mass mr1 = maxFPLMass1 *  fplK1 * fplPropRem1;    // UnSpendable Remnant

  auto [dVhor, propMassS] = MkEstimates(hS, lS, VS, psiS, me1, mr1, IspVac1);

  if (a_os != nullptr && optLogLevel >= 2)
    *a_os << "# Estimated propMassS = " << propMassS.Magnitude()
          << " kg, dVHor=" << dVhor.Magnitude() << " km/sec" << std::endl;

  RTLS1 proto
  (
    maxFPLMass1,   fplK1, fplPropRem1, IspSL1, IspVac1, thrustVacI1,
    minThrtL1,     diam,
    propMassS,     hS,    lS,     VS,  psiS,   dVhor,
    odeIntegrStep, a_os,  optLogLevel
  );

  //-------------------------------------------------------------------------//
  // Create the NOMAD Optimiser:                                             //
  //-------------------------------------------------------------------------//
  // Create the Main NOMAD obj:
  NOMAD::MainStep opt;

  // Create and set the NOMAD params:
  // There are 3 Constraints: (LandMissL, MaxQ, MaxLongG),
  // but their actual vals are not required at this point;
  // XXX: the Landing Velocity is NOT a formal constraint; we make it near-0 by
  // calibrating the LandBurnH and LandBurnGamma, but do not include  it in the
  // list of constraints to avoid the instabilities in NOMAD:
  //
  constexpr int NC = 3;

  // LoBounds and UpBounds are All-0s and All-1s, resp.:
  std::vector<double> loBounds(initParamsN.size(), 0.0);
  std::vector<double> upBounds(initParamsN.size(), 1.0);

  std::shared_ptr<NOMAD::AllParameters> params =
    MkNOMADParams(initParamsN, loBounds, upBounds,       NC,
                  optMaxEvals, optSeed,  stopIfFeasible, useVNS, useMT);

  opt.setAllParameters(params);

  // Create the "NOMADMainEvaluator":
  std::unique_ptr<NOMADMainEvaluator> ev
    (new NOMADMainEvaluator
    (
      params->getEvalParams(),
      &proto,
      landDLLimit, QLimit, longGLimit
    ));
  opt.setEvaluator(std::move(ev));

  //-------------------------------------------------------------------------//
  // RUN the Optimiser:                                                      //
  //-------------------------------------------------------------------------//
  opt.start();
  opt.run();
  opt.end();

  //-------------------------------------------------------------------------//
  // Extract the results:                                                    //
  //-------------------------------------------------------------------------//
  std::vector<NOMAD::EvalPoint> feasPts;
  (void) NOMAD::CacheBase::getInstance()->findBestFeas(feasPts);

  if (feasPts.empty())
    // No feasible solution has been found, "intParamsN" is invalid:
    return std::make_pair(std::nullopt, std::nullopt);

  // If OK: Put the best solution back to "initParamsN":
  NOMAD::EvalPoint bestF = feasPts[0];

  for (int unsigned j = 0; j < initParamsN.size(); ++j)
    initParamsN[j] = bestF[j].todouble();

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
    initParamsN[10]   // entryBurnThrtL1N,
  );
  static_assert(10 == NP-1);

  if (withFinalRun)
  {
    //-----------------------------------------------------------------------//
    // Perform the Final Run on the optimised params (in "rtls"):            //
    //-----------------------------------------------------------------------//
    rtls.m_logLevel     = finalRunLogLevel;
    Base::RunRes runRes = rtls.Run(RTLS1::FlightMode::Coast);

    // NB: "optRes" needs to reflect the dynamic effects of the Funal Run:
    std::optional<OptRes> optRes = OptRes(rtls);
    return std::make_pair(optRes,  runRes);
  }

  // Otherwise, there is no "runRes", and "optRes" is taken w/o the dynamic
  // effects:
  std::optional<OptRes>   optRes = OptRes(rtls);
  return   std::make_pair(optRes, std::nullopt);
}

//===========================================================================//
// "MkEstimates" (to narrow-down the optimisation domain):                   //
//===========================================================================//
// Returns (dVhor, propMassS) estimates.
// XXX: They are computed W/O taking any atmospheric drag effects into account!
//
std::pair<VelK, Mass> RTLS1::MkEstimates
(
  LenK  a_hS,
  LenK  a_lS,
  VelK  a_VS,
  Angle a_psiS,
  Mass  a_empty_mass,
  Mass  a_unspendable_mass,
  Time  a_IspVac1
)
{
  // The Horizontal Delta-V required to fly from the Separation point to (0,0):
  // First, the Separation Conds:
  VelK  VrS   = a_VS * Sin(a_psiS);
  VelK  VhorS = a_VS * Cos(a_psiS);
  if (!(IsPos(a_hS) && IsPos(a_lS) && IsPos(VrS)  && IsPos(VhorS) &&
        IsPos(a_empty_mass)        && IsPos(a_unspendable_mass)))
    throw std::invalid_argument("RTLS1::MkEstimates: Invalid Param(s)");

  // XXX: Using some approximate average "gFall":
  AccK  gFall = 0.5 * (g0K  + K  / Sqr(R + a_hS));

  VelK  dVhor, dVr;  // Not known yet
  Time  tFall;       // ditto

  // Proceed with a small number of iterations:
  for (int i = 0; i < 3; ++i)
  {
    VelK Vr  = VrS + dVr;

    // The time of the ballistic flight to the surface from the Separation
    // Point;
    tFall    = Vr / gFall + SqRt(Sqr(Vr / gFall) + 2.0 * a_hS / gFall);
    assert(IsPos(tFall));

    // Then, assuming the constant horozontal velocity, the following one will
    // be required to get to l=0 at the time of h=0:
    VelK  Vback = a_lS / tFall;
    assert(IsPos(Vback));

    // And including the existing "VhorS" we need to compensate, the following
    // is the total Horizontal DeltaV:
    dVhor = VhorS + Vback;

    // XXX: The actuall Boost-Back Burn may also give us a small vertical "dVr",
    // which we don't know yet,  but assume it to be no more than ~10% of the
    // "dvHor" -- this is an estimate only:
    dVr   = 0.1 * dVhor;
  }
  // So: Got "dVhor", "dVr" estimates. Then the total DeltaV required for fly-
  // back and soft landing must also include the "velocity equivalent" of the
  // altitude "hS":
  //
  VelK dV = SqRt(Sqr(dVr) + Sqr(dVhor) + K * (2.0 / R - 2.0 / (R + a_hS)));

  // Then the corresp estimated "propMassS" (INCLUDING the unspendable remnant),
  // by the Tsiolkovsky formula:
  Mass propMassS =
    (a_empty_mass +  a_unspendable_mass) *
    Exp(double(dV / (g0K * a_IspVac1))) - a_empty_mass;

  assert (propMassS > a_unspendable_mass);
  return std::make_pair(dVhor, propMassS);
}
}
// End namespace SpaceBallistics

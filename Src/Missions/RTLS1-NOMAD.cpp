// vim:ts=2:et
//===========================================================================//
//                       "Src/Missions/RTLS1-NOMAD.cpp":                     //
//              Return-to-(Launch/Landing)-Site: NOMAD Interface             //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"

// XXX: NOMAD headers produce tons of warnings, suppress them all:
#ifdef   __clang__
#pragma  clang diagnostic push
#pragma  clang diagnostic ignored "-Wunused-parameter"
#pragma  clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma  clang diagnostic ignored "-Wcast-qual"
#pragma  clang diagnostic ignored "-Wcast-align"
#pragma  clang diagnostic ignored "-Wold-style-cast"
#pragma  clang diagnostic ignored "-Wsuggest-override"
#pragma  clang diagnostic ignored "-Wsuggest-destructor-override"
#pragma  clang diagnostic ignored "-Woverloaded-virtual"
#pragma  clang diagnostic ignored "-Wshadow"
#pragma  clang diagnostic ignored "-Wextra-semi"
#pragma  clang diagnostic ignored "-Wunused-member-function"
#pragma  clang diagnostic ignored "-Winconsistent-missing-destructor-override"
#pragma  clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
#pragma  clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
#pragma  clang diagnostic ignored "-Wdeprecated-dynamic-exception-spec"
#pragma  clang diagnostic ignored "-Wreserved-identifier"
#pragma  clang diagnostic ignored "-Wheader-hygiene"
#pragma  clang diagnostic ignored "-Wsign-conversion"
#pragma  clang diagnostic ignored "-Warray-bounds"
#pragma  clang diagnostic ignored "-Wmissing-noreturn"
#pragma  clang diagnostic ignored "-Wexit-time-destructors"
#pragma  clang diagnostic ignored "-Wglobal-constructors"
#pragma  clang diagnostic ignored "-Wundefined-func-template"
#else
#pragma  GCC   diagnostic push
#pragma  GCC   diagnostic ignored "-Wunused-parameter"
#pragma  GCC   diagnostic ignored "-Woverloaded-virtual="
#pragma  GCC   diagnostic ignored "-Wnon-virtual-dtor"
#pragma  GCC   diagnostic ignored "-Wcast-qual"
#pragma  GCC   diagnostic ignored "-Warray-bounds"
#pragma  GCC   diagnostic ignored "-Wold-style-cast"
#pragma  GCC   diagnostic ignored "-Wcast-align"
#endif
#include <Nomad/nomad.hpp>
#include <Cache/CacheBase.hpp>
#ifdef   __clang__
#pragma  clang diagnostic pop
#else
#pragma  GCC   diagnostic pop
#endif

#include "SpaceBallistics/Missions/MkNOMADParams.hpp"

namespace SpaceBallistics
{
//===========================================================================//
// "RTLS1::NOMADEvaluator": Helper Class used in NOMAD Optimisation:         //
//===========================================================================//
class RTLS1::NOMADEvaluator final: public NOMAD::Evaluator
{
private:
	//=========================================================================//
	// Data Flds:																															 //
	//=========================================================================//
	RTLS1 const*      const m_proto;
	// Optimisation Limits (Constraints):
	LenK              const m_landDLLimit;
  VelK              const m_landVLimit;
  Pressure          const m_QLimit;

public:
  //=========================================================================//
  // Non-Default Ctor, Dtor:                                                 //
  //=========================================================================//
  NOMADEvaluator
  (
    std::shared_ptr<NOMAD::EvalParameters> const& a_params,
    RTLS1 const*                                  a_proto,
    // Optimisation Limits (Constraints):
    LenK                                          a_land_dL_limit,
    VelK                                          a_land_V_limit,
    Pressure                                      a_Q_limit
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_landDLLimit   (a_land_dL_limit),
    m_landVLimit    (a_land_V_limit),
    m_QLimit        (a_Q_limit)
  {
    assert(m_proto != nullptr);
    if (!(IsPos(m_landDLLimit) && IsPos(m_landVLimit) && IsPos(m_QLimit)))
      throw std::invalid_argument
            ("RTLS1::NOMADEvaluator::Ctor: Invalid Limit(s)");
  }

  ~NOMADEvaluator() override {}

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
    assert(8 <= a_x.size() && a_x.size() <= NP);

    //-----------------------------------------------------------------------//
    // For Thread-Safety, construct a new "RTLS1" obj from "m_proto":        //
    //-----------------------------------------------------------------------//
    Ascent2 rtls(*m_proto);

    // Install the curr Ctl params in "rtls":
    rtls.SetCtlParams
    (
      a_x[0].todouble(),                              // propMassSN
      a_x[1].todouble(),                              // coastDurN
      a_x[2].todouble(),                              // bbBurnDurN
      a_x[3].todouble(),                              // entryBurnQN
      a_x[4].todouble(),                              // entryBurnDurN
      a_x[5].todouble(),                              // landBurnHN
      a_x[6].todouble(),                              // landBurnThrtN
      // sin(theta) coeffs:
      a_x[7].todouble(),                              // sinTheta0
      (a_x.size() >=  9) ? a_x[ 8].todouble() : 0.5,  // sinTheta1
      (a_x.size() >= 10) ? a_x[ 9].todouble() : 0.5,  // sinTheta2
      (a_x.size() == 11) ? a_x[10].todouble() : 0.5   // sinTheta3
    );

    //-----------------------------------------------------------------------//
    // Run the Integrator!                                                   //
    //-----------------------------------------------------------------------//
    RunRes res = rtls.Run();  // NB: Exceptions are handled inside "Run"

    // IMPORTANT: NOMAD allows us to indicate that evaluation has failed:
    if (res.m_rc == RunRC::Error)
      return false;

    //-----------------------------------------------------------------------//
    // Push the results back to NOMAD in string form:                        //
    //-----------------------------------------------------------------------//
    char  buff[512];
    char* curr = buff;

    // The Objective Function Value (to me minimised): PropMassS which was pre-
    // viously put into "rtls":
    curr += sprintf(curr, "%.16e", res.m_propMassS);

    // Constraints:

    return true;
  }

public:
};

//===========================================================================//
// "RunNOMAD":                                                               //
//===========================================================================//
bool RTLS1::RunNOMAD
(
  // Main Optimisation Problem Setup:
  RTLS1 const*            a_proto,
  std::vector<double>*    a_init_vals,
  // Optimisation Constraints (Limits):
  LenK                    a_land_dL_limit,
  VelK                    a_land_V_limit,
  Pressure                a_Q_limit,
  // NOMAD Params:
  int                     a_max_evals,
  int                     a_opt_seed,
  bool                    a_stop_if_feasible,
  double                  a_use_vns
)
{
  assert(a_proto != nullptr && a_init_vals != nullptr);

  //-------------------------------------------------------------------------//
  // Generic Case: Perform NOMAD Optimisation:                               //
  //-------------------------------------------------------------------------//
  // Create the Main NOMAD obj:
  NOMAD::MainStep opt;

  // Create and set the NOMAD params:
  // There are 3 Constraints: (LandMissL, LandV, MaxQ), but their actual vals
  // are not required yet;
  // LoBounds and UpBounds are all-0s and all-1s,  resp:
  std::vector<double> loBounds(a_init_vals->size(), 0.0);
  std::vector<double> upBounds(a_init_vals->size(), 1.0);

  std::shared_ptr<NOMAD::AllParameters> params =
    MkNOMADParams(*a_init_vals, loBounds,   upBounds,        3,
                  a_max_evals,  a_opt_seed, a_stop_if_feasible, a_use_vns);

  opt.setAllParameters(params);

  //-----------------------------------------------------------------------//
  // Create and Run the "NOMADEvaluator":                                  //
  //-----------------------------------------------------------------------//
  std::unique_ptr<NOMADEvaluator> ev
    (new NOMADEvaluator
    (
      params->getEvalParams(), a_proto,
      a_land_dL_limit,         a_land_V_limit,  a_Q_limit
    ));
  opt.setEvaluator(std::move(ev));

  // RUN!
  opt.start();
  opt.run();
  opt.end();

  // Extract the results:
  std::vector<NOMAD::EvalPoint> feasPts;
  (void) NOMAD::CacheBase::getInstance()->findBestFeas(feasPts);
  if (feasPts.empty())
    return false;

  // If OK: Put the best solution back to "initVals":
  NOMAD::EvalPoint bestF = feasPts[0];

  for (int unsigned j = 0; j < a_init_vals->size(); ++j)
    (*a_init_vals)[j] = bestF[j].todouble();

  // All Done:
  return true;
}

}
// End namespace SpaceBallistics

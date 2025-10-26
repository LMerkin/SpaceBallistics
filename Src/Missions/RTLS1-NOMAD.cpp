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
#include "MkNOMADParams.hpp"

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
  // LV Params:
	RTLS1 const*      const m_proto;
  Mass              const m_fullMass1;
  double            const m_fullK1;
  double            const m_fullPropRem1;
  Len               const m_diam;
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
    // LV Params:
    RTLS1 const*                                  a_proto,
    Mass                                          a_full_mass1,
    double                                        a_full_k1,
    double                                        a_full_prop_rem1,
    Len                                           a_diam,
    // Optimisation Limits (Constraints):
    LenK                                          a_land_dL_limit,
    VelK                                          a_land_V_limit,
    Pressure                                      a_Q_limit
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_fullMass1     (a_full_mass1),
    m_fullK1        (a_full_k1),
    m_fullPropRem1  (a_full_prop_rem1),
    m_diam          (a_diam),
    m_landDLLimit   (a_land_dL_limit),
    m_landVLimit    (a_land_V_limit),
    m_QLimit        (a_Q_limit)
  {
    assert(m_proto != nullptr   && IsPos(m_fullMass1)   &&
           0.0 < m_fullK1       && m_fullK1       < 1.0 &&
           0.0 < m_fullPropRem1 && m_fullPropRem1 < 1.0 && IsPos(m_diam));

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

    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 2)
    {
      (*m_proto->m_os)   << '#';
      for (int i = 0; i < int(a_x.size()); ++i)
        (*m_proto->m_os) << "  " << a_x[i].todouble();
      (*m_proto->m_os)   << std::endl;
    }

    //-----------------------------------------------------------------------//
    // For Thread-Safety, construct a new "RTLS1" obj:                       //
    //-----------------------------------------------------------------------//
    // (In contrast to "Ascent2", here we do not clone "m_proto" and modify the
    // clone's params; rather, we construct a new obj using the params memoised
    // in this class and in "m_proto");
    //
    // IMPORTANT: HERE "PropMassS" is set!
    Mass propMassS = a_x[0].todouble() * MaxPropMassS;

    RTLS1 rtls
    (
      // Stage Params:
      m_fullMass1,
      m_fullK1,
      m_fullPropRem1,
      m_proto->Base::m_IspSL1,
      m_proto->Base::m_IspVac1,
      m_proto->Base::m_thrustVacI1,
      m_proto->Base::m_minThrtL1,
      m_diam,
      propMassS,          // NB: Variable!
      // Mission Params:
      m_proto->m_hS,
      m_proto->m_lS,
      m_proto->m_VS,
      m_proto->m_phiS,
      // Integration and Output Params:
      m_proto->Base::m_odeIntegrStep,
      m_proto->Base::m_os,
      m_proto->Base::m_logLevel
    );

    // Install the curr Ctl params in "rtls":
    rtls.SetCtlParams
    (
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

    // The Objective Function Value (to me minimised): It is "propMassS" itself
    // (provided that all constraints are satisfied):
    curr += sprintf(curr, "%.16e",  propMassS.Magnitude());

    // Constraints:
    // Down-Range Miss (the aiming point is 0):
    LenK dL = Abs(res.m_LT);
    curr += sprintf(curr, " %.16e", (dL         - m_landDLLimit).Magnitude());

    // The Landing Velocity (the target is 0):
    assert(!IsNeg(res.m_VT));
    curr += sprintf(curr, " %.16e", (res.m_VT   - m_landVLimit).Magnitude());

    // XXX: Possibly TODO: Landing Acceleration as well

    // The MaxQ encountered:
    assert(!IsNeg(res.m_maxQ));
    curr += sprintf(curr, " %.16e", (res.m_maxQ - m_QLimit).Magnitude());

    // Output done!
    assert(size_t(curr - buff) <= sizeof(buff));

    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 2)
    {
#     pragma omp critical(NOMADOutput)
      *(m_proto->m_os) << buff << std::endl;
    }
    // Set the results back in "a_x":
    a_x.setBBO(buff);
    a_countEval = true;
    return        true;
  }

public:
};

//===========================================================================//
// "RunNOMAD":                                                               //
//===========================================================================//
bool RTLS1::RunNOMAD
(
  // LV Params:
  RTLS1 const*            a_proto,
  Mass                    a_full_mass1,
  double                  a_full_k1,
  double                  a_full_prop_rem1,
  Len                     a_diam,
  // Optimisation Params:
  std::vector<double>*    a_init_vals,
  double                  a_min_prop_massSN,
  double                  a_min_bbb_durN,
  LenK                    a_land_dL_limit,
  VelK                    a_land_V_limit,
  Pressure                a_Q_limit,
  // NOMAD Params:
  int                     a_max_evals,
  int                     a_opt_seed,
  bool                    a_stop_if_feasible,
  double                  a_use_vns,
  bool                    a_use_mt
)
{
  assert(a_proto != nullptr      && a_init_vals != nullptr   &&
         0.0 < a_min_prop_massSN && a_min_prop_massSN <= 1.0 &&
         0.0 < a_min_bbb_durN    && a_min_bbb_durN    <= 1.0);

  //-------------------------------------------------------------------------//
  // Generic Case: Perform NOMAD Optimisation:                               //
  //-------------------------------------------------------------------------//
  // Create the Main NOMAD obj:
  NOMAD::MainStep opt;

  // Create and set the NOMAD params:
  // There are 3 Constraints: (LandMissL, LandV, MaxQ), but their actual vals
  // are not required yet;
  // LoBounds and UpBounds are all-0s and all-1s, with some special cases below:
  //
  std::vector<double> loBounds(a_init_vals->size(), 0.0);
  std::vector<double> upBounds(a_init_vals->size(), 1.0);

  // PropMassSN (Idx=0):
  loBounds[0] =
    std::max
    (
      a_min_prop_massSN,
      double(1.01 * a_full_mass1 * a_full_k1 * a_full_prop_rem1 /
             MaxPropMassS)
    );
  (*a_init_vals)[0] = 0.5 * (loBounds[0] + 1.0);

  // BBBurnDurN (Idx=2):
  loBounds[2] = a_min_bbb_durN;
  (*a_init_vals)[2] = 0.5 * (loBounds[2] + 1.0);

  std::shared_ptr<NOMAD::AllParameters> params =
    MkNOMADParams(*a_init_vals, loBounds,   upBounds,        3,
                  a_max_evals,  a_opt_seed, a_stop_if_feasible, a_use_vns,
                  a_use_mt);

  opt.setAllParameters(params);

  //-----------------------------------------------------------------------//
  // Create and Run the "NOMADEvaluator":                                  //
  //-----------------------------------------------------------------------//
  std::unique_ptr<NOMADEvaluator> ev
    (new NOMADEvaluator
    (
      params->getEvalParams(),
      a_proto,         a_full_mass1,   a_full_k1, a_full_prop_rem1, a_diam,
      a_land_dL_limit, a_land_V_limit, a_Q_limit
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

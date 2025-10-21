// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/Ascent2-NOMAD.cpp":                     //
//         Ascent-to-Orbit for a "Model" 2-Stage LV: NOMAD Interface         //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"

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
// "Ascent2::NOMADEvaluator": Helper Class used in NOMAD Optimisation:       //
//===========================================================================//
class Ascent2::NOMADEvaluator final: public NOMAD::Evaluator
{
private:
  //=========================================================================//
  // Data Flds:                                                              //
  //=========================================================================//
  Ascent2 const*      const m_proto;
  std::array<bool,NP> const m_actOpts;     // Flags for Active Opt Params
  // Optimisation Limits (Constraints):
  VelK                const  m_startVLimit;
  Pressure            const  m_QLimit;     // NAN -> disabled
  Pressure            const  m_sepQLimit;  // ditto
  double              const  m_longGLimit; // ditto

public:
  //=========================================================================//
  // Non-Default Ctor, Dtor:                                                 //
  //=========================================================================//
  NOMADEvaluator
  (
    std::shared_ptr<NOMAD::EvalParameters> const& a_params,
    Ascent2             const*                    a_proto,
    // Flags indicating which variables are optimisation args:
    std::array<bool,NP> const&                    a_act_opts,
    // Optimisation Limits (Constraints):
    VelK                                          a_startV_limit,
    Pressure                                      a_Q_limit,
    Pressure                                      a_sepQ_limit,
    double                                        a_longG_limit
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_actOpts       (a_act_opts),
    m_startVLimit   (a_startV_limit),
    m_QLimit        (a_Q_limit),
    m_sepQLimit     (a_sepQ_limit),
    m_longGLimit    (a_longG_limit)
  {
    assert(m_proto != nullptr);

    if (!(IsPos(m_startVLimit) && IsPos(m_QLimit)   &&
          IsPos(m_sepQLimit))  && IsPos(m_longGLimit))
      throw std::invalid_argument
            ("Ascent2::NOMADEvaluator::Ctor: Invalid Limit(s)");
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
    //-----------------------------------------------------------------------//
    // Defaults for all Optimisation Params (unless over-written):           //
    //-----------------------------------------------------------------------//
    // The over-all list of params is:
    // [thrustMult2, bHat2, muHat2, aAoAHat2, bAoAHat2, TGapRel,
    //  thrustMult1, bHat1, muHat1, aAoAHat1, bAoAHat1, alpha1, payLoadMassRel],
    // where all params are Dim-Less:
    //
    double thrustMult2  = m_proto->m_thrustMult2;
    double bHat2        = m_proto->m_bHat2;
    double muHat2       = m_proto->m_muHat2;
    double aAoAHat2     = m_proto->m_aAoAHat2;
    double bAoAHat2     = m_proto->m_bAoAHat2;
    Time   TGap         = m_proto->m_TGap;
    double thrustMult1  = m_proto->m_thrustMult1;
    double bHat1        = m_proto->m_bHat1;
    double muHat1       = m_proto->m_muHat1;
    double aAoAHat1     = m_proto->m_aAoAHat1;
    double bAoAHat1     = m_proto->m_bAoAHat1;
    double alpha1       = m_proto->m_alpha1;
    Mass   payLoadMass  = m_proto->m_payLoadMass;

    //-----------------------------------------------------------------------//
    // Extract the actual params from "a_x":                                 //
    //-----------------------------------------------------------------------//
    int j = 0;
    for (int i = 0; i < NP; ++i)
    {
      if (!m_actOpts[size_t(i)])
        continue;

      // Otherwise: the "i"th over-all idx is mapped to the "j"th idx in "a_x":
      double xj = a_x[size_t(j)].todouble();

      switch (i)
      {
      case  0: thrustMult2 = xj; break;
      case  1: bHat2       = xj; break;
      case  2: muHat2      = xj; break;
      case  3: aAoAHat2    = xj; break;
      case  4: bAoAHat2    = xj; break;
      case  5: TGap        = xj * MaxTGap;                 break;
      case  6: thrustMult1 = xj; break;
      case  7: bHat1       = xj; break;
      case  8: muHat1      = xj; break;
      case  9: aAoAHat1    = xj; break;
      case 10: bAoAHat1    = xj; break;
      case 11: alpha1      = xj; break;
      case 12: payLoadMass = xj * m_proto->m_maxStartMass; break;
      default: assert(false);
      }
      ++j;
    }
    //-----------------------------------------------------------------------//
    // For thread-safety, construct a new "Ascent2" obj from "m_proto":      //
    //-----------------------------------------------------------------------//
    Ascent2 asc(*m_proto);

    // Install the possibly-modified params in "asc":
    asc.ModifyLVParams
    (
      // Thrust and Mass Params:
      thrustMult2, thrustMult1,
      alpha1,      payLoadMass,
      // Ctl Params:
      bHat2,       muHat2,  aAoAHat2,  bAoAHat2, TGap,
      bHat1,       muHat1,  aAoAHat1,  bAoAHat1
    );

    //-----------------------------------------------------------------------//
    // Run the Integrator!                                                   //
    //-----------------------------------------------------------------------//
    RunRes res = asc.Run(); // NB: Exceptions are handled inside "Run"

    // IMPORTANT: NOMAD allows us to indicate that evaluation has failed:
    if (res.m_rc == RunRC::Error)
      return false;

    //-----------------------------------------------------------------------//
    // Push the results back to NOMAD in string form:                        //
    //-----------------------------------------------------------------------//
    char  buff[512];
    char* curr = buff;

    // XXX: The Objective Function Value: StartMass or (-PayLoadMass), both to
    // be minimised. Which one is to be used, depends on whether the PayLoadMass
    // is the list of Optimisation Params:
    bool objIsStartMass = !m_actOpts[NP-1];
    Mass objMass        =
      objIsStartMass    ? res.m_mT : (- asc.m_payLoadMass);
    curr += sprintf(curr, "%.16e", objMass.Magnitude());

    // Constraints:
    // StartV (or a StartV equivalent taking StartH into account):
    VelK VT = std::max(res.m_VT, VelK(0.0));
    assert(!(IsNeg(res.m_maxQ) || IsNeg(res.m_sepQ) || IsNeg(res.m_maxLongG)));

    curr += sprintf(curr, " %.16e", (VT - m_startVLimit).Magnitude());
    curr += sprintf(curr, " %.16e", (res.m_maxQ     - m_QLimit).Magnitude());
    curr += sprintf(curr, " %.16e", (res.m_sepQ     - m_sepQLimit).Magnitude());
    curr += sprintf(curr, " %.16e", (res.m_maxLongG - m_longGLimit));

    // Output done!
    assert(size_t(curr - buff) <= sizeof(buff));

    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 4)
    {
#     pragma omp critical(NOMADOutput)
      *(m_proto->m_os) << buff << std::endl;
    }
    // Set the results back in "a_x":
    a_x.setBBO(buff);
    a_countEval = true;
    return        true;
  }
};

//===========================================================================//
// "RunNOMAD":                                                               //
//===========================================================================//
bool Ascent2::RunNOMAD
(
  // Main Optimisation Problem Setup:
  Ascent2             const*  a_proto,
  std::array<bool,NP> const&  a_act_opts,
  std::vector<double>*        a_init_vals,
  std::vector<double> const&  a_lo_bounds,
  std::vector<double> const&  a_up_bounds,
  // Optimisation Constraints:
  VelK                        a_startV_limit,
  Pressure                    a_Q_limit,
  Pressure                    a_sepQ_limit,
  double                      a_longG_limit,
  // NOMAD Params:
  int                         a_max_evals,
  int                         a_opt_seed,
  bool                        a_stop_if_feasible,
  double                      a_use_vns
)
{
  assert(a_proto != nullptr && a_init_vals != nullptr);
  size_t np       = a_init_vals->size();
  assert(0 < np && np <= NP &&
         a_lo_bounds.size() == np && a_up_bounds.size() == np);

  //-------------------------------------------------------------------------//
  // Generic Case: Perform NOMAD Optimisation:                               //
  //-------------------------------------------------------------------------//
  // Create the Main NOMAD obj:
  NOMAD::MainStep opt;

  // Create and set the NOMAD params:
  // There are 4 Constraints: (StartV, MaxQ, SepQ, MaxLongG), but their actual
  // vals are not required yet:
  std::shared_ptr<NOMAD::AllParameters> params =
    MkNOMADParams(*a_init_vals, a_lo_bounds, a_up_bounds,     4,
                  a_max_evals,  a_opt_seed,  a_stop_if_feasible, a_use_vns);

  opt.setAllParameters(params);

  //-------------------------------------------------------------------------//
  // Create and Run the "NOMADEvaluator":                                    //
  //-------------------------------------------------------------------------//
  std::unique_ptr<NOMADEvaluator> ev
    (new NOMADEvaluator
    (
      params->getEvalParams(),
      a_proto,        a_act_opts,
      a_startV_limit, a_Q_limit, a_sepQ_limit, a_longG_limit
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

  for (int unsigned j = 0; j < unsigned(np); ++j)
    (*a_init_vals)[j] = bestF[j].todouble();

  // The "optMass" found (either the mininimised StartMass, or the maximised
  // PayLoadMass):
  // Mass optMass(bestF.getF(NOMAD::defaultFHComputeType).todouble());

  // All Done:
  return true;
}
}
// End namespace SpaceBallistics

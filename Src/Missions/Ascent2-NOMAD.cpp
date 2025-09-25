// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/Ascent2-Optim.cpp":                     //
//     Ascent-to-Orbit for a "Model" 2-Stage LV: Parametric Optimisation     //
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
#else
#pragma  GCC   diagnostic push
#pragma  GCC   diagnostic ignored "-Wunused-parameter"
#pragma  GCC   diagnostic ignored "-Woverloaded-virtual="
#pragma  GCC   diagnostic ignored "-Wnon-virtual-dtor"
#pragma  GCC   diagnostic ignored "-Wcast-qual"
#pragma  GCC   diagnostic ignored "-Warray-bounds"
#endif
#include <Nomad/nomad.hpp>
#include <Cache/CacheBase.hpp>
#ifdef   __clang__
#pragma  clang diagnostic pop
#else
#pragma  GCC   diagnostic pop
#endif

namespace SpaceBallistics
{
//===========================================================================//
// "Ascent2::NOMADEvaluator": Helper Class used in NOMAD Optimisation:       //
//===========================================================================//
class Ascent2::NOMADEvaluator final: public NOMAD::Evaluator
{
private:
  //-------------------------------------------------------------------------//
  // Data Flds:                                                              //
  //-------------------------------------------------------------------------//
  Ascent2 const*  m_proto;
  bool            m_actOpts[NP];    // Flags for Active Optimisation Params

public:
  //=========================================================================//
  // Non-Default Ctor, Dtor:                                                 //
  //=========================================================================//
  NOMADEvaluator
  (
    Ascent2 const*                                a_proto,
    bool    const                                 a_act_opts[NP],
    std::shared_ptr<NOMAD::EvalParameters> const& a_params
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto)
  {
    assert(m_proto != nullptr && a_act_opts != nullptr);
    for (int i = 0; i < NP; ++i)
      m_actOpts[i]  = a_act_opts[i];
  }

  ~NOMADEvaluator() override {}

  //=========================================================================//
  // "eval_x": The actual evaluation method for NOMAD:                       //
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
      if (!m_actOpts[i])
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
    // For thread-safety, construct a new "Ascent2" obj from "m_proto":
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
  
    // Run the Integrator!
    RunRes res = asc.Run();  // NB: Exceptions are handled inside "Run"
  
    // Analyse the Results:
    if (res.m_rc == RunRC::Error)
      return false;
  
    // StartH:
    LenK hT = std::max(res.m_hT, 0.0_km);
  
    // StartV:
    VelK VT = std::max(res.m_VT, VelK(0.0));
  
    // Push the results back to NOMAD in string form:
    char  buff[512];
    char* curr  = buff;
  
    // XXX: The Objective Function Value: StartMass or (-PayLoadMass), both to
    // be minimised. Which one is to be used, depends on whether the PayLoadMass
    // is the list of Optimisation Params:
    bool objIsStartMass = !m_actOpts[NP-1];
    Mass objMass        =
      objIsStartMass    ? res.m_mT : (- asc.m_payLoadMass);
    curr += sprintf(curr, "%.16e", objMass.Magnitude());
  
    // Constraint0: StartH:
    curr += sprintf(curr, " %.16e ", (hT - Ascent2::MaxStartH).Magnitude());
  
    // Constraint1: StartV:
    curr += sprintf(curr, " %.16e",  (VT - Ascent2::MaxStartV).Magnitude());
  
    // Other Constraints if enabled:
    if (IsFinite(asc.m_QLimit))
      curr +=
        sprintf(curr, " %.16e", (res.m_maxQ - asc.m_QLimit).Magnitude());
  
    if (IsFinite(asc.m_sepQLimit))
      curr +=
        sprintf(curr, " %.16e", (res.m_sepQ - asc.m_sepQLimit).Magnitude());
  
    if (IsFinite(asc.m_longGLimit))
      curr +=
      sprintf(curr, " %.16e",   (res.m_maxLongG - asc.m_longGLimit));
  
    // Output done!
    assert(size_t(curr - buff) <= sizeof(buff));
  
    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 1)
    {
  #   pragma omp critical(NOMADOutput)
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
  Ascent2*                            a_proto,
  bool                const           a_act_opts[NP],
  std::vector<double>*                a_init_vals,
  std::vector<double> const&          a_lo_bounds,
  std::vector<double> const&          a_up_bounds,
  int                                 a_max_evals,
  bool                                a_constr_q,
  bool                                a_constr_sep_q,
  bool                                a_constr_long_g,
  boost::property_tree::ptree const&  a_pt
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

  // Create the Params:
  // Problem Geometry:
  auto params = std::make_shared<NOMAD::AllParameters>();
  params->setAttributeValue("DIMENSION",         np);
  params->setAttributeValue("X0",                NOMAD::Point(*a_init_vals));
  params->setAttributeValue
    ("LOWER_BOUND", NOMAD::ArrayOfDouble(a_lo_bounds));
  params->setAttributeValue
    ("UPPER_BOUND", NOMAD::ArrayOfDouble(a_up_bounds));

  // Stopping Criterion: XXX: Currently, only via the MaxEvals param:
  params->setAttributeValue("MAX_BB_EVAL", a_max_evals);

  // XXX: The following must be compatible with "NOMADEvaluator::eval_x":
  NOMAD::BBOutputTypeList  bbTypes;
  bbTypes.push_back(NOMAD::BBOutputType::OBJ);

  // 2 constraints on the start consitions: h(-T) and V(-T),
  // to be satisfied at the solution point only:
  bbTypes.push_back(NOMAD::BBOutputType::PB);
  bbTypes.push_back(NOMAD::BBOutputType::PB);

  // Possible extra constraints: MaxQ, MaxSepQ, MaxLongG,
  // to be satisfied at the solution point only as well:
  if (a_constr_q)
    bbTypes.push_back(NOMAD::BBOutputType::PB);
  if (a_constr_sep_q)
    bbTypes.push_back(NOMAD::BBOutputType::PB);
  if (a_constr_long_g)
    bbTypes.push_back(NOMAD::BBOutputType::PB);

  params->setAttributeValue("BB_OUTPUT_TYPE", bbTypes );

  // Parallel Evaluation: The number of threads is decided automatically:
  params->setAttributeValue("NB_THREADS_PARALLEL_EVAL", -1);

  // "DIRECTION_TYPE" selects a variant of the optimisation algorithm. Other
  // possible vals include "ORTHO_NP1_NEG", "ORTHO_NP1_QUAD", "NP1_UNI"  and
  // many others:
  params->setAttributeValue("DIRECTION_TYPE",
    NOMAD::DirectionType::ORTHO_2N);

  // Other params:
  params->setAttributeValue("DISPLAY_DEGREE",          2);
  params->setAttributeValue("DISPLAY_ALL_EVAL",        false);
  params->setAttributeValue("DISPLAY_UNSUCCESSFUL",    false);
  params->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES",  false);
  params->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);

  int    optSeed = a_pt.get<int>   ("Technical.NOMADSeed");
  params->setAttributeValue("SEED", optSeed);

  double useVNS  = a_pt.get<double>("Technical.NOMADUseVNS");
  if (useVNS <  0.0 || useVNS >= 1.0)
    throw std::invalid_argument ("NOMADUseVNS: The arg must be in [0..1)");

  if (useVNS != 0.0)
  {
    params->setAttributeValue("VNS_MADS_SEARCH",         true);
    params->setAttributeValue("VNS_MADS_SEARCH_TRIGGER",
                              NOMAD::Double(useVNS));
  }

  // Validate the "params" and install them in the "opt":
  params->checkAndComply();
  opt.setAllParameters(params);

  //-----------------------------------------------------------------------//
  // Create and Run the "NOMADEvaluator":                                  //
  //-----------------------------------------------------------------------//
  std::unique_ptr<NOMADEvaluator> ev
    (new NOMADEvaluator
    (
      a_proto,
      a_act_opts,
      params->getEvalParams()
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

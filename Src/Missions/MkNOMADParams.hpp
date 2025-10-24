// vim:ts=2:et
//===========================================================================//
//                        "Missions/MkNOMADParams.hpp":                      //
//               Common Part of NOMAD Interface Initialisation               //
//===========================================================================//
// XXX: This header can only be included in the context of <Nomad/nomad.hpp>;
// however,  it is not included automatically   because it requires multiple
// pragmas to suppess cimpiler warnings
#pragma  once
#include <vector>

namespace SpaceBallistics
{
//===========================================================================//
// "MkNOMADParams":                                                          //
//===========================================================================//
inline std::shared_ptr<NOMAD::AllParameters> MkNOMADParams
(
  std::vector<double> const& a_init_vals,
  std::vector<double> const& a_lo_bounds,
  std::vector<double> const& a_up_bounds,
  int                        a_n_constrs,
  int                        a_max_evals,
  int                        a_opt_seed,
  bool                       a_stop_if_feasible,
  double                     a_use_vns,
  bool                       a_use_mt
)
{
  assert(a_init_vals.size() == a_lo_bounds.size() &&
         a_init_vals.size() == a_up_bounds.size());

  // Create the Params:
  // Problem Geometry:
  auto params = std::make_shared<NOMAD::AllParameters>();
  params->setAttributeValue("DIMENSION",         a_init_vals.size()),
  params->setAttributeValue("X0",                NOMAD::Point(a_init_vals));
  params->setAttributeValue
    ("LOWER_BOUND", NOMAD::ArrayOfDouble(a_lo_bounds));
  params->setAttributeValue
    ("UPPER_BOUND", NOMAD::ArrayOfDouble(a_up_bounds));

  // Stopping Criterion: XXX: Currently, only via the MaxEvals param:
  params->setAttributeValue("MAX_BB_EVAL", a_max_evals);

  // XXX: The following must be compatible with "NOMADEvaluator::eval_x":
  // The Objective Function to be Minimised:
  NOMAD::BBOutputTypeList  bbTypes;
  bbTypes.push_back(NOMAD::BBOutputType::OBJ);

  // And the Constraints:
  for (int i = 0; i < a_n_constrs; ++i)
    bbTypes.push_back(NOMAD::BBOutputType::PB);

  params->setAttributeValue("BB_OUTPUT_TYPE", bbTypes );

  // Parallel Evaluation: If enabled, the number of threads is decided
  // automatically:
  params->setAttributeValue("NB_THREADS_PARALLEL_EVAL", a_use_mt ? -1 : 1);

  // "DIRECTION_TYPE" selects a variant of the optimisation algorithm. Other
  // possible vals include "ORTHO_NP1_NEG", "ORTHO_NP1_QUAD", "NP1_UNI"  and
  // many others:
  params->setAttributeValue("DIRECTION_TYPE",
    NOMAD::DirectionType::ORTHO_2N);

  // Other params:
  params->setAttributeValue("DISPLAY_DEGREE",          2);    // TODO: Config!
  params->setAttributeValue("DISPLAY_ALL_EVAL",        false);
  params->setAttributeValue("DISPLAY_UNSUCCESSFUL",    false);
  params->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES",  false);
  params->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);

  params->setAttributeValue("SEED",             a_opt_seed);
  params->setAttributeValue("STOP_IF_FEASIBLE", a_stop_if_feasible);

  if (a_use_vns <  0.0 || a_use_vns >= 1.0)
    throw std::invalid_argument ("NOMADUseVNS: The arg must be in [0..1)");

  if (a_use_vns != 0.0)
  {
    assert(0.0 < a_use_vns && a_use_vns < 1.0);
    params->setAttributeValue("VNS_MADS_SEARCH",         true);
    params->setAttributeValue("VNS_MADS_SEARCH_TRIGGER",
                              NOMAD::Double(a_use_vns));
  }
  // Validate the "params" and install them in the "opt":
  params->checkAndComply();
  return params;
}
}
// End namespace SpaceBallistics

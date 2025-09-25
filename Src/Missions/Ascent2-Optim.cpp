// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/Ascent2-Optim.cpp":                     //
//     Ascent-to-Orbit for a "Model" 2-Stage LV: Parametric Optimisation     //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"
#include <boost/property_tree/ini_parser.hpp>
#include <cstdlib>

namespace SpaceBallistics
{
//===========================================================================//
// "ModifyLVParams":                                                         //
//===========================================================================//
// This method allows us to modify consistently the Mass and Thrust params in
// the "Ascent2" object. Intended for use with Optimisation Procedures:
//
void Ascent2::ModifyLVParams
(
  // Mass and Thrust Params:
  double    a_thrustMult2,
  double    a_thrustMult1,
  double    a_alpha1,
  Mass      a_payLoadMass,

  // Ctl Params:
  double    a_bHat2,      double a_muHat2,
  double    a_aAoAHat2,   double a_bAoAHat2,
  Time      a_TGap,
  double    a_bHat1,      double a_muHat1,
  double    a_aAoAHat1,   double a_bAoAHat1
)
{
  // Verify the Mass and Thrust Params (the Ctl Params and Constraints are
  // verified later in "SetCtlParams"):
  if (a_thrustMult2    <= 0.0 || a_thrustMult1    <= 0.0 ||
      a_alpha1         <  1.0 || IsNeg(a_payLoadMass))
    throw std::invalid_argument
          ("Ascent2::ModifyLVParams: Invalid Mass and/or Thrust Params");

  //-------------------------------------------------------------------------//
  // Set the Mass and Thrust Params:                                         //
  //-------------------------------------------------------------------------//
  // Over-All:
  m_alpha1       = a_alpha1;
  m_payLoadMass  = a_payLoadMass;

  // Stage2:
  m_fullMass2    = (m_maxStartMass   - m_fairingMass - m_payLoadMass) /
                    (1.0 + m_alpha1);
  m_emptyMass2   = m_fullMass2   * (1.0 - m_K2);
  m_propMass2    = m_fullMass2   * m_K2;
  m_unSpendable2 = m_propMass2   * m_propRem2;
  m_spendable2   = m_propMass2   - m_unSpendable2;
  m_thrustVacI2 *= a_thrustMult2 / m_thrustMult2;
  m_thrustMult2  = a_thrustMult2;
  m_burnRateI2   = m_thrustVacI2 / (m_IspVac2 * g0K);
  m_T2           = m_spendable2  / m_burnRateI2;

  // Stage1:
  m_fullMass1    = (m_maxStartMass   - m_fairingMass - m_payLoadMass) /
                    (1.0 + m_alpha1) * m_alpha1;
  m_emptyMass1   = m_fullMass1   * (1.0 - m_K1);
  m_propMass1    = m_fullMass1   * m_K1;
  m_unSpendable1 = m_propMass1   * m_propRem1;
  m_spendable1   = m_propMass1   - m_unSpendable1;    
  m_thrustVacI1 *= a_thrustMult1 / m_thrustMult1;
  m_thrustMult1  = a_thrustMult1;
  m_burnRateI1   = m_thrustVacI1 / (m_IspVac1 * g0K);
  m_T1           = m_spendable1  / m_burnRateI1;

  // The above settings also invalidate the curr Ctls, so we need to modify
  // them anyway:
  SetCtlParams(a_bHat2, a_muHat2, a_aAoAHat2, a_bAoAHat2, a_TGap,
               a_bHat1, a_muHat1, a_aAoAHat1, a_bAoAHat1);
}
//===========================================================================//
// "SetCtlParams":                                                           //
//===========================================================================//
void Ascent2::SetCtlParams
(
  double   a_bHat2,     double   a_muHat2,
  double   a_aAoAHat2,  double   a_bAoAHat2,
  Time     a_TGap,
  double   a_bHat1,     double   a_muHat1,
  double   a_aAoAHat1,  double   a_bAoAHat1
)
{
  //-------------------------------------------------------------------------//
  // Check what we got:                                                      //
  //-------------------------------------------------------------------------//
  if (!(0.0 <= a_bHat2       && a_bHat2    <= 1.0  &&
        0.0 <= a_muHat2      && a_muHat2   <= 1.0  &&
        0.0 <= a_aAoAHat2    && a_aAoAHat2 <= 1.0  &&
        0.0 <= a_bAoAHat2    && a_bAoAHat2 <= 1.0  &&
        !IsNeg(a_TGap)                             &&
        0.0 <= a_bHat1       && a_bHat1    <= 1.0  &&
        0.0 <= a_muHat1      && a_muHat1   <= 1.0  &&
        0.0 <= a_aAoAHat1    && a_aAoAHat1 <= 1.0  &&
        0.0 <= a_bAoAHat1    && a_bAoAHat1 <= 1.0))
      throw std::invalid_argument("Ascent2::SetCtlParams: Invalid Param(s)");

  //-------------------------------------------------------------------------//
  // Ballistic Gap:                                                          //
  //-------------------------------------------------------------------------//
  m_TGap    = a_TGap;

  //-------------------------------------------------------------------------//
  // BurnRate Coeffs:                                                        //
  //-------------------------------------------------------------------------//
  // Must have (2*MinThrtL+1)/3 <= muHat <= 1. XXX: We currently only allow
  // "down-throttling", not "up-rating" of the Engines:
  //
  double           muHatLo2 = (2.0 * m_minThrtL2 + 1.0) / 3.0;
  double           muHatLo1 = (2.0 * m_minThrtL1 + 1.0) / 3.0;
  constexpr double muHatUp  =  1.0;

  // Stage2:
  m_muHat2 = a_muHat2;
  double muHat2  = (1.0 - a_muHat2) * muHatLo2 + a_muHat2 * muHatUp;
  assert           (0.0 < muHat2  &&  muHat2 <= 1.0);

  double bHatLo2 = 3.0 * muHat2 * (muHat2 - 1.0);
  assert(bHatLo2 <= 0.0);
  double bHatUp2 = 2.0 * muHat2 * (3.0 * muHat2 - (2.0 + m_minThrtL2));
  bHatUp2  = std::min(bHatUp2, 0.0);

  // NB: We must have bHatLo2 <= bHatUp2, up to rounding errors:
  m_bHat2        = a_bHat2;
  double bHat2   = (1.0 - a_bHat2) * bHatLo2 + a_bHat2 * bHatUp2;
  assert(bHat2  <=  0.0);

  // Then set the actual (dimensioned) coeffs for Stage2 (using the Spendable
  // PropMass2, ie one w/o the Remnants):
  m_T2     = m_spendable2      / (m_burnRateI2 * muHat2);
  m_bMu2   = Sqr(m_burnRateI2) /  m_spendable2 * bHat2;
  m_aMu2   = 3.0 / Cube(m_T2)  *
             (m_spendable2 - m_burnRateI2 * m_T2 - 0.5 * m_bMu2 * Sqr(m_T2));

  // Stage1:
  m_muHat1       =  a_muHat1;
  double muHat1  =  (1.0 - a_muHat1) * muHatLo1 + a_muHat1 * muHatUp;
  assert            (0.0 < muHat1  &&  muHat1 <= 1.0);

  double bHatLo1 =  3.0 * muHat1 * (muHat1 - 1.0);
  assert(bHatLo1 <= 0.0);
  double bHatUp1 =  2.0 * muHat1 * (3.0 * muHat1 - (2.0 + m_minThrtL1));
  bHatUp1  = std::min(bHatUp1, 0.0);

  // NB: We must have bHatLo1 <= bHatUp1, up to rounding errors:
  m_bHat1  = a_bHat1;
  double bHat1   = (1.0 - a_bHat1) * bHatLo1 + a_bHat1 * bHatUp1;
  assert(bHat1  <=  0.0);

  // Then set the actual (dimensioned) coeffs for Stage1 (using the Spendable
  // PropMass1, ie one w/o the Remnants):
  m_T1     = m_spendable1      / (m_burnRateI1 * muHat1);
  m_bMu1   = Sqr(m_burnRateI1) /  m_spendable1 * bHat1;
  m_aMu1   = 3.0 / Cube(m_T1)  *
             (m_spendable1 - m_burnRateI1 * m_T1 - 0.5 * m_bMu1 * Sqr(m_T1));

  //-------------------------------------------------------------------------//
  // AoA Coeffs:                                                             //
  //-------------------------------------------------------------------------//
  // Stage2:
  m_aAoAHat2     = a_aAoAHat2;
  m_bAoAHat2     = a_bAoAHat2;
  m_bAoA2        = 4.0 * m_maxAoA2 / m_T2 * a_bAoAHat2;
  AngAcc aAoALo2 = - m_bAoA2       / m_T2;
  AngAcc aAoAUp2 =
    (a_bAoAHat2  < 0.5)
    ? (m_maxAoA2     / m_T2 - m_bAoA2) / m_T2
    : - Sqr(m_bAoA2) / (4.0 * m_maxAoA2);
  m_aAoA2        =     (1.0 - a_aAoAHat2) * aAoALo2 + a_aAoAHat2 * aAoAUp2;

  // Stage1:
  m_aAoAHat1     = a_aAoAHat1;
  m_bAoAHat1     = a_bAoAHat1;
  m_bAoA1        = 4.0 * m_maxAoA1 / m_T1 * a_bAoAHat1;
  AngAcc aAoALo1 = - m_bAoA1       / m_T1;
  AngAcc aAoAUp1 =
    (a_bAoAHat1  < 0.5)
    ? (m_maxAoA1     / m_T1 - m_bAoA1) / m_T1
    : - Sqr(m_bAoA1) / (4.0 * m_maxAoA1);
  m_aAoA1        =     (1.0 - a_aAoAHat1) * aAoALo1 + a_aAoAHat1 * aAoAUp1;
}

//===========================================================================//
// "NOMADEvaluator" Ctor:                                                    //
//===========================================================================//
Ascent2::NOMADEvaluator::NOMADEvaluator
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

//===========================================================================//
// "NOMADEvaluator::eval_x":                                                 //
//===========================================================================//
bool Ascent2::NOMADEvaluator::eval_x
(
  NOMAD::EvalPoint&    a_x,  // Not "const" -- the result is also set here
  NOMAD::Double const& UNUSED_PARAM(a_hMax),
  bool&                a_countEval
)
const
{
  //-------------------------------------------------------------------------//
  // Defaults for all Optimisation Params (unless over-written):             //
  //-------------------------------------------------------------------------//
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

  //-------------------------------------------------------------------------//
  // Extract the actual params from "a_x":                                   //
  //-------------------------------------------------------------------------//
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

  // XXX: The Objective Function Value: StartMass or (-PayLoadMass), both to be
  // minimised. Which one is to be used, depends on whether the PayLoadMass  is
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
    curr += sprintf(curr, " %.16e", (res.m_maxQ - asc.m_QLimit).Magnitude());

  if (IsFinite(asc.m_sepQLimit))
    curr += sprintf(curr, " %.16e", (res.m_sepQ - asc.m_sepQLimit).Magnitude());

  if (IsFinite(asc.m_longGLimit))
    curr += sprintf(curr, " %.16e", (res.m_maxLongG - asc.m_longGLimit));

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
  return true;
}

//===========================================================================//
// "OptRes" Non-Default Ctor:                                                //
//===========================================================================//
Ascent2::OptRes::OptRes(Ascent2 const& a_asc)
{
  // For Stage2:
  m_thrustMult2 = a_asc.m_thrustMult2;
  m_thrustVacI2 = a_asc.m_thrustVacI2;
  m_burnRateI2  = a_asc.m_burnRateI2;
  m_bHat2       = a_asc.m_bHat2;
  m_muHat2      = a_asc.m_muHat2;
  m_T2          = a_asc.m_T2;
  m_aMu2        = a_asc.m_aMu2;
  m_bMu2        = a_asc.m_bMu2;
  m_aAoAHat2    = a_asc.m_aAoAHat2;
  m_bAoAHat2    = a_asc.m_bAoAHat2;
  m_aAoA2       = a_asc.m_aAoA2;
  m_bAoA2       = a_asc.m_bAoA2;

  // Ballistic Gap:
  m_TGap        = a_asc.m_TGap;

  // For Stage1:
  m_thrustMult1 = a_asc.m_thrustMult1;
  m_thrustVacI1 = a_asc.m_thrustVacI1;
  m_burnRateI1  = a_asc.m_burnRateI1;
  m_bHat1       = a_asc.m_bHat1;
  m_muHat1      = a_asc.m_muHat1;
  m_T1          = a_asc.m_T1;
  m_aMu1        = a_asc.m_aMu1;
  m_bMu1        = a_asc.m_bMu1;
  m_aAoAHat1    = a_asc.m_aAoAHat1;
  m_bAoAHat1    = a_asc.m_bAoAHat1;
  m_aAoA1       = a_asc.m_aAoA1;
  m_bAoA1       = a_asc.m_bAoA1;

  // Over-All:
  m_alpha1      = a_asc.m_alpha1;
  m_payLoadMass = a_asc.m_payLoadMass;
}

//===========================================================================//
// "OptRes" Output:                                                          //
//===========================================================================//
std::ostream& operator<< (std::ostream& a_os, Ascent2::OptRes const& a_res)
{
  return
  a_os
    // Dim-Less Params (NP):
    << "\tthrustMult2    = "     << a_res.m_thrustMult2
    << "\n\tbHat2          = "   << a_res.m_bHat2
    << "\n\tmuHat2         = "   << a_res.m_muHat2
    << "\n\taAoAHat2       = "   << a_res.m_aAoAHat2
    << "\n\tbAoAHat2       = "   << a_res.m_bAoAHat2
    << "\n\tthrustMult1    = "   << a_res.m_thrustMult1
    << "\n\tbHat1          = "   << a_res.m_bHat1
    << "\n\tmuHat1         = "   << a_res.m_muHat1
    << "\n\taAoAHat1       = "   << a_res.m_aAoAHat1
    << "\n\tbAoAHat1       = "   << a_res.m_bAoAHat1
    << "\n\talpha1         = "   << a_res.m_alpha1
    // Dimensioned Params:
    << "\n\tthrustVacI2    = "   << a_res.m_thrustVacI2
    << "\n\tT2             = "   << a_res.m_T2
    << "\n\taMu2           = "   << a_res.m_aMu2
    << "\n\tbMu2           = "   << a_res.m_bMu2
    << "\n\taAoA2          = "   << a_res.m_aAoA2
    << "\n\tbAoA2          = "   << a_res.m_bAoA2
    << "\n\tTGap           = "   << a_res.m_TGap
    << "\n\tthrustVacI1    = "   << a_res.m_thrustVacI1
    << "\n\tT1             = "   << a_res.m_T1
    << "\n\taMu1           = "   << a_res.m_aMu1
    << "\n\tbMu1           = "   << a_res.m_bMu1
    << "\n\taAoA1          = "   << a_res.m_aAoA1
    << "\n\tbAoA1          = "   << a_res.m_bAoA1
    << "\n\tpayLoadMass    = "   << a_res.m_payLoadMass
    // XXX: Ctl Eqs: Similar to "Ascent2::OutputCtls", see there for the
    // comments:
    << "\n\t# tau2 = 0 .. "      << a_res.m_T2.Magnitude()
    << "\n\t# t2  = -"           << a_res.m_T2.Magnitude()
    << " .. 0"
       "\n\tAoA2 := t2 * ("      << a_res.m_aAoA2.Magnitude()
    << " * t2 - ("               << a_res.m_bAoA2.Magnitude()      << "));"
       "\n\tmu2  := "            << a_res.m_burnRateI2.Magnitude() << " + ("
    << a_res.m_bMu2.Magnitude()  << ") * tau2 + ("
    << a_res.m_aMu2.Magnitude()  << ") * tau2^2;"
    << "\n\t# tau1 = 0 .. "      << a_res.m_T1.Magnitude()
    << "\n\tAoA1 := tau1 * ("    << a_res.m_aAoA1.Magnitude()
    << " * tau1 + ("             << a_res.m_bAoA1.Magnitude()      << ")); "
       "\n\tmu1  := "            << a_res.m_burnRateI1.Magnitude() << " + ("
    << a_res.m_bMu1.Magnitude()  << ") * tau1 + ("
    << a_res.m_aMu1.Magnitude()  << ") * tau1^2;"
    << std::endl;
}

//===========================================================================//
// "FindOptimalAscentCtls"                                                   //
//===========================================================================//
std::pair<std::optional<Ascent2::OptRes>,
          std::optional<Ascent2::RunRes>> // If the FinalRun is performed
Ascent2::FindOptimalAscentCtls
(
  std::string const& a_config_ini,
  std::ostream*      a_os       // May be NULL
)
{
  //-------------------------------------------------------------------------//
  // Bounds and Initial Vals for all Params:                                 //
  //-------------------------------------------------------------------------//
  // [thrustMult2, bHat2, muHat2, aAoAHat2, bAoAHat2, TGapRel,
  //  thrustMult1, bHat1, muHat1, aAoAHat1, bAoAHat1, alpha1, payLoadMassRel]:
  constexpr double LoBounds[NP]
    { 0.7,         0.0,   0.0,    0.0,      0.0,      0.0,
      0.7,         0.0,   0.0,    0.0,      0.0,      1.0,    0.02 };

  constexpr double UpBounds[NP]
    { 1.5,         1.0,   1.0,    1.0,      1.0,      1.0,
      1.5,         1.0,   1.0,    1.0,      1.0,      7.0,    0.06 };

  constexpr double InitVals[NP]
    { 1.0,         0.5,   0.8,    0.5,      0.5,      0.1,
      1.0,         0.5,   0.8,    0.5,      0.5,      4.0,    0.04 };

  //-------------------------------------------------------------------------//
  // Open and Parse the Config.ini File:                                     //
  //-------------------------------------------------------------------------//
  boost::property_tree::ptree  pt;
  boost::property_tree::ini_parser::read_ini(a_config_ini, pt);

  //-------------------------------------------------------------------------//
  // Get the LV and Mission params for the "prototype" "Ascent2" obj:        //
  //-------------------------------------------------------------------------//
  // Stage2:
  double    K2          =      pt.get<double>("LV.K2");
  double    PropRem2    =      pt.get<double>("LV.PropRem2");
  Time      IspVac2           (pt.get<double>("LV.IspVac2"));
  ForceK    ThrustVacI2 = Mass(pt.get<double>("LV.ThrustVacI2")) * g0K;
  double    MinThrttL2  =      pt.get<double>("LV.MinThrttL2");
  Angle_deg MaxAoA2           (pt.get<double>("LV.MaxAoA2"));
  // Stage1:
  double    K1          =      pt.get<double>("LV.K1");
  double    PropRem1    =      pt.get<double>("LV.PropRem1");
  Time      IspSL1            (pt.get<double>("LV.IspSL1" ));
  Time      IspVac1           (pt.get<double>("LV.IspVac1"));
  ForceK    ThrustVacI1 = Mass(pt.get<double>("LV.ThrustVacI1")) * g0K;
  double    MinThrttL1  =      pt.get<double>("LV.MinThrttL1");
  Angle_deg MaxAoA1           (pt.get<double>("LV.MaxAoA1"));
  // Over-All:
  Mass   MaxStartMass         (pt.get<double>("LV.MaxStartMass"));
  Mass   FairingMass          (pt.get<double>("LV.FairingMass" ));
  Len    Diam                 (pt.get<double>("LV.Diameter"    ));
  // NB: "alpha1" and "payLoadMass" belong to the "Opt" section (see below).

  // Mission:
  LenK      Perigee           (pt.get<double>("Mission.Perigee"));
  LenK      Apogee            (pt.get<double>("Mission.Apogee" ));
  Angle_deg Incl              (pt.get<double>("Mission.Inclination"));
  Angle_deg LaunchLat         (pt.get<double>("Mission.LaunchLat"));

  //-------------------------------------------------------------------------//
  // Now  get the "Opt" params:                                              //
  //-------------------------------------------------------------------------//
  // (*) If a param is give a numerical value in Config.ini, it is considered
  //     to be "fixed" (ie NOT an optimisation argument), and the corresp bit
  //     in "actOpts" is set to "false";
  // (*) if a param has the special value "OPTIMISE", is will serve as an opt-
  //     imisation argument, so the corresp "actOpts" bit is set to "true";
  // (*) XXX: For the user convenience, in "Config.ini", the values of "TGap"
  //     and "payLoadMass" are given as actual magnitudes (in sec and kg, resp),
  //     rather than by relative vals in [0..1]; whereas NOMAD always operates
  //     with relative vals for the sake of uniformity:
  bool                actOpts[NP];
  std::vector<double> loBounds;
  std::vector<double> upBounds;
  std::vector<double> initVals;

# ifdef  GetOptParam
# undef  GetOptParam
# endif
# define GetOptParam(I, Type, Name, Scale)  \
  std::string Name##Str = pt.get<std::string>("Opt." #Name); \
  Type        Name(NAN);       \
  if (Name##Str == "OPTIMISE") \
  { \
    /* XXX: "InitVals" are Dim-Less, so need to apply the Scale: */ \
    Name       = InitVals[I] * Scale; \
    actOpts[I] = true;         \
    loBounds.push_back(LoBounds[I]); \
    upBounds.push_back(UpBounds[I]); \
    initVals.push_back(InitVals[I]); \
  }    \
  else \
  { \
    /* NB: The value for "Name" already has the right magnitude: */ \
    Name       = Type(std::atof(Name##Str.data())); \
    actOpts[I] = false; \
  } \
  assert(IsFinite(Name));

  GetOptParam( 0, double, thrustMult2, 1.0)
  GetOptParam( 1, double, bHat2,       1.0)
  GetOptParam( 2, double, muHat2,      1.0)
  GetOptParam( 3, double, aAoAHat2,    1.0)
  GetOptParam( 4, double, bAoAHat2,    1.0)
  GetOptParam( 5, Time,   TGap,        MaxTGap)
  GetOptParam( 6, double, thrustMult1, 1.0)
  GetOptParam( 7, double, bHat1,       1.0)
  GetOptParam( 8, double, muHat1,      1.0)
  GetOptParam( 9, double, aAoAHat1,    1.0)
  GetOptParam(10, double, bAoAHat1,    1.0)
  GetOptParam(11, double, alpha1,      1.0)
  GetOptParam(12, Mass,   PayLoadMass, MaxStartMass)

  // So: How many Optimisation Args have we got?
  int np = 0;
  for (int i = 0; i < NP; ++i)
    np += int(actOpts[i]);

  assert(int(loBounds.size()) == np && int(upBounds.size()) == np &&
         int(initVals.size()) == np);

  // Extra Constraints: "inf" and "nan" are also allowed, hence the use of
  // "std::atof":
  Pressure QLimit
           (std::atof(pt.get<std::string>("Opt.QLimit"    ).data()));

  Pressure SepQLimit
           (std::atof(pt.get<std::string>("Opt.SepQLimit" ).data()));

  double   LongGLimit =
            std::atof(pt.get<std::string>("Opt.LongGLimit").data());

  //-------------------------------------------------------------------------//
  // Finally: Technical Params:                                              //
  //-------------------------------------------------------------------------//
  // Generic:
  bool   withFinalRun     = pt.get<bool>  ("Technical.WithFinalRun");
  int    finalRunLogLevel = pt.get<int>   ("Technical.FinalRunLogLevel", 3);
  int    optMaxEvals      = pt.get<int>   ("Technical.OptMaxEvals");

  // NOMAD-specific:
  int    optSeed          = pt.get<int>   ("Technical.OptSeed");
  int    optLogLevel      = pt.get<int>   ("Technical.OptLogLevel");
  double useVNS           = pt.get<double>("Technical.OptUseVNS",      0.0);

  if (useVNS < 0.0 || useVNS >= 1.0)
    throw std::invalid_argument("OptUseVNS: The arg must be in [0..1)");

  //-------------------------------------------------------------------------//
  // Create the "prototype" "Ascent2" obj:                                   //
  //-------------------------------------------------------------------------//
  Ascent2 proto
  (
    K2, PropRem2,          IspVac2, ThrustVacI2, MinThrttL2, MaxAoA2,
    K1, PropRem1,  IspSL1, IspVac1, ThrustVacI1, MinThrttL1, MaxAoA1,
    alpha1,  MaxStartMass, FairingMass,    Diam, PayLoadMass,
    QLimit,  SepQLimit,    LongGLimit,
    Perigee, Apogee, Incl, LaunchLat,      a_os, optLogLevel
  );

  proto.SetCtlParams
  (
    bHat2, muHat2, aAoAHat2, bAoAHat2, TGap,
    bHat1, muHat1, aAoAHat1, bAoAHat1
  );

  // NB: Depending on the "actOpts", "optVal" is the minimised StartMass, or
  // the maximised PayLoadMass. For now, it is undefined:
  Mass optMass(NAN);

  if (np != 0)
  {
    //-----------------------------------------------------------------------//
    // Generic Case: Perform NOMAD Optimisation:                             //
    //-----------------------------------------------------------------------//
    // Create the Main NOMAD obj:
    NOMAD::MainStep opt;

    // Create the Params:
    // Problem Geometry:
    auto params = std::make_shared<NOMAD::AllParameters>();
    params->setAttributeValue("DIMENSION",         size_t(np));
    params->setAttributeValue("X0",                NOMAD::Point(initVals));
    params->setAttributeValue
      ("LOWER_BOUND", NOMAD::ArrayOfDouble(loBounds));
    params->setAttributeValue
      ("UPPER_BOUND", NOMAD::ArrayOfDouble(upBounds));

    // Stopping Criterion: XXX: Currently, only via the OptMaxEvals:
    params->setAttributeValue("MAX_BB_EVAL",       optMaxEvals);

    // XXX: The following must be compatible with "NOMADEvaluator::eval_x":
    NOMAD::BBOutputTypeList  bbTypes;
    bbTypes.push_back(NOMAD::BBOutputType::OBJ);

    // 2 constraints on the start consitions: h(-T) and V(-T),
    // to be satisfied at the solution point only:
    bbTypes.push_back(NOMAD::BBOutputType::PB);
    bbTypes.push_back(NOMAD::BBOutputType::PB);

    // Possible extra constraints: MaxQ, MaxSepQ, MaxLongG,
    // to be satisfied at the solution point only as well:
    if (IsFinite(QLimit))
      bbTypes.push_back(NOMAD::BBOutputType::PB);
    if (IsFinite(SepQLimit))
      bbTypes.push_back(NOMAD::BBOutputType::PB);
    if (IsFinite(LongGLimit))
      bbTypes.push_back(NOMAD::BBOutputType::PB);

    params->setAttributeValue("BB_OUTPUT_TYPE",    bbTypes );

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

    params->setAttributeValue("SEED", optSeed);
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
        &proto,
        actOpts,
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
      return std::make_pair(std::nullopt, std::nullopt);

    // If OK: Put the best solution back to "initVals":
    NOMAD::EvalPoint bestF = feasPts[0];

    for (int unsigned j = 0; j < unsigned(np); ++j)
      initVals[j] = bestF[j].todouble();

    // The "optMass" found (either the mininimised StartMass, or the maximised
    // PayLoadMass):
    double optVal = bestF.getF(NOMAD::defaultFHComputeType).todouble();
    optMass       = Mass(optVal);
  }
  // Otherwise (if np=0), no optimisation is performed, but the post-processing
  // still is:

  //-------------------------------------------------------------------------//
  // Post-Processing:                                                        //
  //-------------------------------------------------------------------------//
  // In order to get the final "OptRes",  put  the "initVals" (now containing
  // the opt args found) into the "proto" obj, and extract the "optRes"  from
  // the "proto":
  //
  for (int i = 0, j = 0; i < NP; ++i)
  {
    if (!actOpts[i])
      continue;
    assert(j < np);

    // XXX: The previously-declared variables are used here. They either hold
    // fixed values from Config.ini, or assigned the optimal vals found:
    switch (i)
    {
    case  0: thrustMult2 = initVals[unsigned(j)]; break;
    case  1: bHat2       = initVals[unsigned(j)]; break;
    case  2: muHat2      = initVals[unsigned(j)]; break;
    case  3: aAoAHat2    = initVals[unsigned(j)]; break;
    case  4: bAoAHat2    = initVals[unsigned(j)]; break;
    case  5: TGap        = initVals[unsigned(j)] * MaxTGap;              break;
    case  6: thrustMult1 = initVals[unsigned(j)]; break;
    case  7: bHat1       = initVals[unsigned(j)]; break;
    case  8: muHat1      = initVals[unsigned(j)]; break;
    case  9: aAoAHat1    = initVals[unsigned(j)]; break;
    case 10: bAoAHat1    = initVals[unsigned(j)]; break;
    case 11: alpha1      = initVals[unsigned(j)]; break;
    case 12: PayLoadMass = initVals[unsigned(j)] * proto.m_maxStartMass; break;
    default: assert(false);
    }
    ++j;
  }

  proto.ModifyLVParams
  (
    // Thrust and Mass Params: 
    thrustMult2, thrustMult1, alpha1,    PayLoadMass,
    // Ctl Params:      
    bHat2,       muHat2,      aAoAHat2,  bAoAHat2,  TGap,
    bHat1,       muHat1,      aAoAHat1,  bAoAHat1
  );

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

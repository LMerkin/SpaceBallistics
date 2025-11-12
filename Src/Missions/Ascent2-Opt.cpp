// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/Ascent2-OptGen.cpp":                    //
//     Ascent-to-Orbit for a "Model" 2-Stage LV: Parametric Optimisation     //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"
#include "SpaceBallistics/Missions/MkNOMADParams.hpp"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>
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

  // Stage1: XXX: Here we directly modify the "Base" flds which is undesirable:
  Base::m_fullMass1    =
    (m_maxStartMass   - m_fairingMass - m_payLoadMass) / (1.0 + m_alpha1) *
     m_alpha1;
  Base::m_emptyMass1   = Base::m_fullMass1   * (1.0 - Base::m_K1);
  Base::m_propMass1    = Base::m_fullMass1   * Base::m_K1;
  Base::m_unSpendable1 = Base::m_propMass1   * Base::m_propRem1;
  Base::m_spendable1   = Base::m_propMass1   - Base::m_unSpendable1;
  Base::m_thrustVacI1 *= a_thrustMult1       / m_thrustMult1;
  m_thrustMult1        = a_thrustMult1;
  Base::m_burnRateI1   = Base::m_thrustVacI1 / (Base::m_IspVac1 * g0K);
  Base::m_T1           = Base::m_spendable1  /  Base::m_burnRateI1;

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
  double           muHatLo2 = (2.0 * m_minThrtL2       + 1.0) / 3.0;
  double           muHatLo1 = (2.0 * Base::m_minThrtL1 + 1.0) / 3.0;
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
  double bHatUp1 =  2.0 * muHat1 * (3.0 * muHat1 - (2.0 + Base::m_minThrtL1));
  bHatUp1  = std::min(bHatUp1, 0.0);

  // NB: We must have bHatLo1 <= bHatUp1, up to rounding errors:
  m_bHat1  = a_bHat1;
  double bHat1   = (1.0 - a_bHat1) * bHatLo1 + a_bHat1 * bHatUp1;
  assert(bHat1  <=  0.0);

  // Then set the actual (dimensioned) coeffs for Stage1 (using the Spendable
  // PropMass1, ie one w/o the Remnants):
  Base::m_T1 = Base::m_spendable1 / (Base::m_burnRateI1 * muHat1);
  m_bMu1     = Sqr(m_burnRateI1)  /  Base::m_spendable1 * bHat1;
  m_aMu1     = 3.0 / Cube(m_T1)   *
               (Base::m_spendable1 - Base::m_burnRateI1 * Base::m_T1 -
                0.5 * m_bMu1 * Sqr(Base::m_T1));

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
  m_bAoA1        = 4.0 * m_maxAoA1 / Base::m_T1 * a_bAoAHat1;
  AngAcc aAoALo1 = - m_bAoA1       / Base::m_T1;
  AngAcc aAoAUp1 =
    (a_bAoAHat1  < 0.5)
    ? (m_maxAoA1     / Base::m_T1  - m_bAoA1) / Base::m_T1
    : - Sqr(m_bAoA1) / (4.0 * m_maxAoA1);
  m_aAoA1        =     (1.0 - a_aAoAHat1) * aAoALo1 + a_aAoAHat1 * aAoAUp1;
}

//===========================================================================//
// "OptRes" Non-Default Ctor:                                                //
//===========================================================================//
Ascent2::OptRes::OptRes(Ascent2 const& a_asc)
:
  // For Stage2:
  m_thrustMult2(a_asc.m_thrustMult2),
  m_thrustVacI2(a_asc.m_thrustVacI2),
  m_burnRateI2 (a_asc.m_burnRateI2),
  m_bHat2      (a_asc.m_bHat2),
  m_muHat2     (a_asc.m_muHat2),
  m_T2         (a_asc.m_T2),
  m_aMu2       (a_asc.m_aMu2),
  m_bMu2       (a_asc.m_bMu2),
  m_aAoAHat2   (a_asc.m_aAoAHat2),
  m_bAoAHat2   (a_asc.m_bAoAHat2),
  m_aAoA2      (a_asc.m_aAoA2),
  m_bAoA2      (a_asc.m_bAoA2),

  // Ballistic Gap:
  m_TGap       (a_asc.m_TGap),

  // For Stage1:
  m_thrustMult1(a_asc.m_thrustMult1),
  m_thrustVacI1(a_asc.m_thrustVacI1),
  m_burnRateI1 (a_asc.m_burnRateI1),
  m_bHat1      (a_asc.m_bHat1),
  m_muHat1     (a_asc.m_muHat1),
  m_T1         (a_asc.m_T1),
  m_aMu1       (a_asc.m_aMu1),
  m_bMu1       (a_asc.m_bMu1),
  m_aAoAHat1   (a_asc.m_aAoAHat1),
  m_bAoAHat1   (a_asc.m_bAoAHat1),
  m_aAoA1      (a_asc.m_aAoA1),
  m_bAoA1      (a_asc.m_bAoA1),

  // Over-All:
  m_alpha1     (a_asc.m_alpha1),
  m_payLoadMass(a_asc.m_payLoadMass)
{}

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
// "Ascent2::NOMADEvaluator": Helper Class used in NOMAD Optimisation:       //
//===========================================================================//
class Ascent2::NOMADEvaluator final: public NOMAD::Evaluator
{
private:
  //=========================================================================//
  // Data Flds:                                                              //
  //=========================================================================//
  Ascent2 const*      const  m_proto;
  std::array<bool,NP> const  m_actOpts;    // Flags for Active Opt Params
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
    // Done!
    assert(size_t(curr - buff) <= sizeof(buff));

    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 2)
    {
#     pragma omp critical(Output)
      *(m_proto->m_os) << buff << std::endl;
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
std::pair<std::optional<Ascent2::OptRes>,
          std::optional<Ascent2::Base::RunRes>> // If the FinalRun is performed
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
  double LoBounds[NP]
    { 0.7,         0.0,   0.0,    0.0,      0.0,      0.0,
      0.7,         0.0,   0.0,    0.0,      0.0,      1.0,    0.02 };

  double UpBounds[NP]
    { 1.5,         1.0,   1.0,    1.0,      1.0,      1.0,
      1.5,         1.0,   1.0,    1.0,      1.0,      7.0,    0.06 };

  double InitVals[NP]
    { 1.0,         0.5,   0.8,    0.5,      0.5,      0.1,
      1.0,         0.5,   0.8,    0.5,      0.5,      4.0,    0.04 };

  //-------------------------------------------------------------------------//
  // Open and Parse the Config.ini File:                                     //
  //-------------------------------------------------------------------------//
  boost::property_tree::ptree  pt;
  boost::property_tree::ini_parser::read_ini(a_config_ini, pt);

  // Get the LV and Mission params for the "prototype" "Ascent2" obj:
  // Stage2:
  double    k2          =      pt.get<double>("LV.K2");
  double    propRem2    =      pt.get<double>("LV.PropRem2");
  Time      IspVac2           (pt.get<double>("LV.IspVac2"));
  ForceK    thrustVacI2 = Mass(pt.get<double>("LV.ThrustVacI2")) * g0K;
  double    minThrtL2   =      pt.get<double>("LV.MinThrtL2");
  Angle_deg maxAoA2           (pt.get<double>("LV.MaxAoA2"));

  // Stage1:
  double    k1          =      pt.get<double>("LV.K1");
  double    propRem1    =      pt.get<double>("LV.PropRem1");
  Time      IspSL1            (pt.get<double>("LV.IspSL1" ));
  Time      IspVac1           (pt.get<double>("LV.IspVac1"));
  ForceK    thrustVacI1 = Mass(pt.get<double>("LV.ThrustVacI1")) * g0K;
  double    minThrtL1   =      pt.get<double>("LV.MinThrtL1");
  Angle_deg maxAoA1           (pt.get<double>("LV.MaxAoA1"));

  // Over-All:
  Mass      maxStartMass      (pt.get<double>("LV.MaxStartMass"));
  Mass      fairingMass       (pt.get<double>("LV.FairingMass" ));
  Len       diam              (pt.get<double>("LV.Diameter"    ));
  // NB: the initial "alpha1" and "payLoadMass" belong to the "Opt" section
  // (see below)...

  // Mission:
  LenK      perigee           (pt.get<double>("Mission.Perigee"));
  LenK      apogee            (pt.get<double>("Mission.Apogee" ));
  Angle_deg incl              (pt.get<double>("Mission.Inclination"));
  Angle_deg launchLat         (pt.get<double>("Mission.LaunchLat"));

  //-------------------------------------------------------------------------//
  // Now  get the "Opt" params:                                              //
  //-------------------------------------------------------------------------//
  // (*) If a param is give a numerical value in Config.ini, it is considered
  //     to be "fixed" (ie NOT an optimisation argument), and the corresp bit
  //     in "actOpts" is set to "false";
  // (*) if a param has the special value "OPTIMISE", is will serve as an opt-
  //     imisation argument, so the corresp "actOpts" bit is set to "true";
  // (*) "OPTIMISE" may be followed by 3 (or none) numerical vals which are
  //     [LoBound InitVal UpBound] of the corresp param; in that case, they
  //     will be used instead of built-in defaults.
  // (*) XXX: For the user convenience, in "Config.ini", the values of "TGap"
  //     and "payLoadMass" are given as actual magnitudes (in sec and kg, resp),
  //     whereas the Optimiser always uses relative vals normalised to [0..1]
  //     for the sake of uniformity:
  std::array<bool,NP> actOpts;
  std::vector<double> loBounds;
  std::vector<double> upBounds;
  std::vector<double> initVals;

# ifdef  GetOptParam
# undef  GetOptParam
# endif
# define GetOptParam(I, Type, Name, Scale) \
  Type   Name(NAN);       \
  { \
    std::string paramStr = pt.get<std::string>("Opt." #Name); \
    if (paramStr.substr(0,8) == "OPTIMISE") \
    { \
      /* Parse the params: */ \
      std::vector<std::string> tokens;  \
      boost::split(tokens, paramStr, boost::is_any_of(" \t")); \
      int nTokens = int(tokens.size()); \
      assert(nTokens >= 1); \
      if (nTokens == 2 || nTokens == 3 || nTokens > 4 || \
          tokens[0] != "OPTIMISE")  \
        throw std::invalid_argument \
          ("Ascent2::FindOptimalAscentCtls: Invalid Params: " + paramStr); \
      /* So yes, the "I"s variable is an optimisation arg: */ \
      actOpts [I] =  true;   \
      if (nTokens == 4)      \
      { \
        /* Set the Range and the InitVal: */ \
        LoBounds[I] = std::atof(tokens[1].data()); \
        InitVals[I] = std::atof(tokens[2].data()); \
        UpBounds[I] = std::atof(tokens[3].data()); \
      } \
      /* For the Dimensioned variable, need to apply the Scale: */ \
      Name  =  InitVals[I] * Scale;    \
      /* The following vectors will be passed to the optimiser: */ \
      loBounds.push_back(LoBounds[I]); \
      upBounds.push_back(UpBounds[I]); \
      initVals.push_back(InitVals[I]); \
    }    \
    else \
    { \
      /* NB: The value for "Name" already has the right magnitude: */ \
      Name       = Type(std::atof(paramStr.data())); \
      actOpts[I] = false; \
    } \
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
  GetOptParam(12, Mass,   payLoadMass, maxStartMass)

  // So: How many Optimisation Args have we got?
  int np = 0;
  for (int i = 0; i < NP; ++i)
    np += int(actOpts[size_t(i)]);

  assert(int(loBounds.size()) == np && int(upBounds.size()) == np &&
         int(initVals.size()) == np);

  // Limits for the Constained Variables:
  VelK     maxStartV (pt.get<double>("Opt.MaxStartV" ));
  Pressure QLimit    (pt.get<double>("Opt.QLimit"    ));
  Pressure sepQLimit (pt.get<double>("Opt.SepQLimit" ));
  double   longGLimit(pt.get<double>("Opt.LongGLimit"));

  //-------------------------------------------------------------------------//
  // Finally: Technical Params:                                              //
  //-------------------------------------------------------------------------//
  // Generic Params:
  Time   odeIntegrStep     (pt.get<double>("Technical.ODEIntegrStep"));
  bool   withFinalRun     = pt.get<bool>  ("Technical.WithFinalRun");
  int    finalRunLogLevel = pt.get<int>   ("Technical.FinalRunLogLevel");
  int    optMaxEvals      = pt.get<int>   ("Technical.OptMaxEvals");
  int    optLogLevel      = pt.get<int>   ("Technical.OptLogLevel");

  // NOMAD-Specific Params:
  int    optSeed          = pt.get<int>   ("Technical.NOMADSeed");
  bool   stopIfFeasible   = pt.get<bool>  ("Technical.NOMADStopIfFeasible");
  double useVNS           = pt.get<double>("Technical.NOMADUseVNS");
  if (useVNS < 0.0 || useVNS >= 1.0)      // 0: VNS not used
    throw std::invalid_argument ("NOMADUseVNS: The arg must be in [0..1)");
  bool   useMT            = pt.get<bool>  ("Technical.NOMADUseMT");

  //-------------------------------------------------------------------------//
  // Create the "prototype" "Ascent2" obj:                                   //
  //-------------------------------------------------------------------------//
  Ascent2 proto
  (
    k2,      propRem2,          IspVac2, thrustVacI2, minThrtL2, maxAoA2,
    k1,      propRem1,  IspSL1, IspVac1, thrustVacI1, minThrtL1, maxAoA1,
    alpha1,  maxStartMass,      fairingMass,    diam,            payLoadMass,
    perigee, apogee,    incl,   launchLat,      odeIntegrStep,   a_os,
    optLogLevel
  );

  proto.SetCtlParams
  (
    bHat2, muHat2, aAoAHat2, bAoAHat2, TGap,
    bHat1, muHat1, aAoAHat1, bAoAHat1
  );

  if (np != 0)
  {
    //-----------------------------------------------------------------------//
    // Create the NOMAD Optimiser:                                           //
    //-----------------------------------------------------------------------//
    // The Main NOMAD obj:
    NOMAD::MainStep opt;

    // Create and set the NOMAD params:
    // There are 4 Constraints: (StartV, MaxQ, SepQ, MaxLongG), but their actual
    // vals are not required yet:
    constexpr int NC = 4;
    std::shared_ptr<NOMAD::AllParameters> params =
      MkNOMADParams(initVals,    loBounds, upBounds,   NC,
                    optMaxEvals, optSeed,  stopIfFeasible, useVNS, useMT);

    opt.setAllParameters(params);

    std::unique_ptr<NOMADEvaluator> ev
      (new NOMADEvaluator
      (
        params->getEvalParams(),
        &proto,    actOpts,
        maxStartV, QLimit, sepQLimit, longGLimit
      ));
    opt.setEvaluator(std::move(ev));

    //-----------------------------------------------------------------------//
    // RUN the Optimiser:                                                    //
    //-----------------------------------------------------------------------//
    opt.start();
    opt.run();
    opt.end();

    //-----------------------------------------------------------------------//
    // Extract the Results:                                                  //
    //-----------------------------------------------------------------------//
    std::vector<NOMAD::EvalPoint> feasPts;
    (void) NOMAD::CacheBase::getInstance()->findBestFeas(feasPts);

    if (feasPts.empty())
      return std::make_pair(std::nullopt, std::nullopt);

    // If OK: Put the best solution back to "initVals":
    NOMAD::EvalPoint bestF = feasPts[0];

    for (int unsigned j = 0; j < unsigned(np); ++j)
      (initVals)[j] = bestF[j].todouble();

    // The "optMass" found (either the mininimised StartMass, or the maximised
    // PayLoadMass):
    // Mass optMass(bestF.getF(NOMAD::defaultFHComputeType).todouble());
  }
  // Otherwise (if np=0), no optimisation is performed, but the post-processing
  // still is...

  //-------------------------------------------------------------------------//
  // Post-Processing:                                                        //
  //-------------------------------------------------------------------------//
  // In order to get the final "OptRes",  put  the "initVals" (now containing
  // the opt args found) into the "proto" obj, and extract the "optRes"  from
  // the "proto":
  //
  for (int i = 0, j = 0; i < NP; ++i)
  {
    if (!actOpts[size_t(i)])
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
    case 12: payLoadMass = initVals[unsigned(j)] * proto.m_maxStartMass; break;
    default: assert(false);
    }
    ++j;
  }

  proto.ModifyLVParams
  (
    // Thrust and Mass Params:
    thrustMult2, thrustMult1, alpha1,    payLoadMass,
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
    proto.m_logLevel    = finalRunLogLevel;
    Base::RunRes runRes = proto.Run();
    return std::make_pair(optRes, runRes);
  }
  // Otherwise, there is no "runRes":
  return std::make_pair(optRes, std::nullopt);
}
}
// End namespace SpaceBallistics

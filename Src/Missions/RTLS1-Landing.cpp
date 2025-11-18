// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/RTLS1-Landing.cpp":                     //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include "SpaceBallistics/Missions/MkNOMADParams.hpp"
#include <boost/property_tree/ini_parser.hpp>
#include <vector>

namespace SpaceBallistics
{
//===========================================================================//
// "LandBurnApprox": Approxmilate Handling of LandBurn Conditions:           //
//===========================================================================//
RTLS1::Base::RunRes RTLS1::LandBurnApprox(LandBurnApproxExn const& a_lba) const
{
  // The "total equivalent velocity" using the Energy Integral:
  VelK V0 = SqRt(Sqr(a_lba.m_Vr) + Sqr(a_lba.m_Vhor) +
                 K * (2.0 / R - 2.0 / (R + MaxLandBurnH)));

  // The total LV mass @ "MaxLandBurnH":
  Mass   m         = Base::m_fullMass1  - a_lba.m_spentPropMass;
  // XXX: The folloiwng inequalities should hold; however, allow for various
  // rounding errors:
  assert(m        >= Base::m_emptyMass1 + Base::m_unSpendable1 - 5.0_kg);
  m    = std::max(m, Base::m_emptyMass1 + Base::m_unSpendable1);

  Mass   propMass  = Base::m_propMass1  - a_lba.m_spentPropMass;
  assert(propMass >= Base::m_unSpendable1 - 5.0_kg);
  propMass = std::max(propMass, Base::m_unSpendable1);

  // The Exhaust Velocity (assuming IspSL):
  VelK   W         = Base::m_IspSL1 * g0K;

  // Using the Tsiolkovsky formula, compute the PropMass required to cancel
  // "V0":
  Mass   propMass0 = m * (1.0 - Exp(- double(V0 / W)));

  // The landing velocity and mass:
  VelK   vL;
  Mass   mL;

  // Do we have sufficient propellant for soft landing?
  if (propMass0 <= propMass - Base::m_unSpendable1)
  {
    // Yes, "propMass" is sufficient, the landing velocity is 0:
    assert(IsZero(vL));
    mL   = m   - propMass0;
    assert(mL >=        Base::m_emptyMass1 + Base::m_unSpendable1 - 5.0_kg);
    mL   = std::max(mL, Base::m_emptyMass1 + Base::m_unSpendable1);
  }
  else
  {
    // No: we have spent the maximum spendable propellant mass and still hit
    // the ground:
    mL      = Base::m_emptyMass1 + Base::m_unSpendable1;
    VelK dV = W * Log(double(m / mL));
    assert(!IsNeg(dV) && dV < V0);
    vL      = V0 - dV;
  }

  // In any case, we reach "ZeroH", with the "vL" velocity. XXX:
  // (*) We do not compute the Fall Time, so the final time is NAN; this does
  //     not affect the optimisation process;
  // (*) More importantly, we still need to provide an approximation of the hor-
  //     izontal landing (or impact) co-ordinate. We assume that both "Vr"  and
  //     "Vhor" are cancelled in the same proportion,  so the landing point  is
  //     simply a projection of the curr "a_lba" point;
  // (*) the landing acceleration is also NAN in this case:
  //
  assert(IsNeg(a_lba.m_Vr));
  LenK   dL = - MaxLandBurnH * a_lba.m_Vhor / a_lba.m_Vr;
  LenK   L  = R * a_lba.m_phi / 1.0_rad;

  return Base::RunRes
        (Base::RunRC::ZeroH, Time(NAN), L + dL, vL, AccK(NAN), mL,
         m_maxQ,     m_sepQ, m_maxLongG);
}

//===========================================================================//
// "SetLandBurnParams":                                                      //
//===========================================================================//
// This function is constructed from the data obtained using
// "CalibrateLandBurnParams".
// It maps the observed velocities (just after crossing the "MaxLandBurnH"
// boundary) to the optimal "LandBurnH" and "LandBurnGamma"  which  should
// provide the soft landing (with Velocity and Acceleration constraints sa-
// tisfied):
//
//
void RTLS1::SetLandBurnParams
(
  Base::StateV const&  a_s,
  Time                 DEBUG_ONLY(  a_t),
  Base::DStateV const& UNUSED_PARAM(a_ds),
  Base::StateV  const& DEBUG_ONLY(  a_prev_s),
  Time                 DEBUG_ONLY(  a_prev_t),
  Base::DStateV const& UNUSED_PARAM(a_prev_ds)
)
{
  //-------------------------------------------------------------------------//
  // By default: No LandBurn: Hit the ground at the full end velocity:       //
  //-------------------------------------------------------------------------//
  m_landBurnH     = 0.0_km;
  m_landBurnGamma = 0.0;

  //-------------------------------------------------------------------------//
  // Curr and Prev State:                                                    //
  //-------------------------------------------------------------------------//
  LenK  r         = std::get<0>(a_s);
  LenK  h         = r      - R;
# ifndef NDEBUG
  LenK  prev_h    = std::get<0>(a_prev_s) - R;
# endif

  VelK  Vr        = std::get<1>(a_s);
  VelK  Vhor      = std::get<2>(a_s) * r / 1.0_rad;
  VelK  V         = SqRt(Sqr(Vr) + Sqr(Vhor));
  Angle psi       = Angle(ATan2(Vr, Vhor));
  Mass  propMass  = Base::m_propMass1 - std::get<3>(a_s);
  assert(IsPos(propMass));

  if (Base::m_os != nullptr && Base::m_logLevel >= 2)
#   pragma omp critical(Output)
    *Base::m_os
      << "# LandBurn Barrier: h=" << h.Magnitude() << " km, V="
      << V.Magnitude() << " km/sec, psi="
      << To_Angle_deg(psi).Magnitude() << " deg, propMass="
      << propMass.Magnitude()          << " kg"    << std::endl;

  // This methos is invoked at the point when we have just crossed the
  // "MaxLandBurnH" boundary, so:
  assert(IsPos(h) && IsPos(prev_h) && a_prev_t <  a_t  && IsNeg(Vhor) &&
         h <= MaxLandBurnH   &&   MaxLandBurnH <= prev_h);
}

//===========================================================================//
// "NOMADCalibrEvaluator" Class: Used in "CalibrateLandBurnParams":          //
//===========================================================================//
class RTLS1::NOMADCalibrEvaluator final: public NOMAD::Evaluator
{
private:
  //=========================================================================//
  // Data Flds:                                                              //
  //=========================================================================//
  // LV Params:
  RTLS1 const*  const m_proto;

  // Initial Conds (@ MaxLandBurnH):
  VelK          const m_VL;
  Angle         const m_psiL;
  Mass          const m_propMassL;

public:
  //=========================================================================//
  // Non-Default Ctor, Dtor:                                                 //
  //=========================================================================//
  NOMADCalibrEvaluator
  (
    std::shared_ptr<NOMAD::EvalParameters> const& a_params,
    RTLS1 const*                                  a_proto,
    VelK                                          a_VL,
    Angle                                         a_psiL,
    Mass                                          a_prop_massL
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_VL            (a_VL),
    m_psiL          (a_psiL),
    m_propMassL     (std::max(a_prop_massL, m_proto->m_unSpendable1))
  {
    assert(m_proto != nullptr);
    if (!IsPos(m_VL))
      throw std::invalid_argument
            ("RTLS1::NOMADCalibrEvaluator::Ctor: Invalid Param(s)");
  }

  ~NOMADCalibrEvaluator() override {}

  //=========================================================================//
  // "eval_x": The Actual Evaluation Method for MOMAD:                       //
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
    // For Thread-Safety, construct a new "RTLS1" obj from "Proto":          //
    //-----------------------------------------------------------------------//
    RTLS1 rtls
    (
      // Stage Params saved in the "Proto":
      m_proto->m_fplMass1,
      m_proto->m_fplK1,
      m_proto->m_fplPropRem1,
      m_proto->Base::m_IspSL1,
      m_proto->Base::m_IspVac1,
      m_proto->Base::m_thrustVacI1,
      m_proto->Base::m_minThrtL1,
      m_proto->m_diam,
      // IMPORTANT: The actual initial "PropMass":
      m_propMassL,

      // The Initial Conds:
      MaxLandBurnH,
      0.0_km,       // The horizontal distance does not matter here
      m_VL,
      m_psiL,

      // Estimates and Limits:
      m_proto->m_dVhorEst,
      m_proto->m_landDLLimit,
      m_proto->m_landVelLimit,
      m_proto->m_landAccLimit,
      m_proto->m_QLimit,
      m_proto->m_longGLimit,
      false,        // Do NOT use the LandBurn Approximation!

      // Optimisation Ranges:
      m_proto->m_propMassSRange,
      m_proto->m_bbBurnThetaMinPi,
      m_proto->m_bbBurnDurRange,
      m_proto->m_entryBurnDurRange,

      // Other Params:
      m_proto->Base::m_odeIntegrStep,
      m_proto->Base::m_os,
      m_proto->Base::m_logLevel
    );

    // Pre-set the LandBurn params:
    assert(a_x.size() == 2);
    double landBurnHN    = a_x[0].todouble();
    double landBurnGamma = a_x[1].todouble();

    assert(0.0 <= landBurnHN    && landBurnHN    <= 1.0 &&
           0.0 <= landBurnGamma && landBurnGamma <= 1.0);

    rtls.m_landBurnH     = landBurnHN * MaxLandBurnH;
    rtls.m_landBurnGamma = landBurnGamma;

    //-----------------------------------------------------------------------//
    // Run the Integrator:                                                   //
    //-----------------------------------------------------------------------//
    // NB: Exceptions are handled inside "Run":
    RunRes res = rtls.Run(RTLS1::FlightMode::EndoAtmDesc);

    // IMPORTANT: NOMAD allows us to indicate that evaluation has failed:
    if (res.m_rc == RunRC::Error)
      return false;

    //-----------------------------------------------------------------------//
    // Push the results back to NOMAD in string form:                        //
    //-----------------------------------------------------------------------//
    char  buff[512];
    char* curr = buff;    

    // The Objective Function Value (to me minimised): It is the PropMass spent:
    Mass propMassSpent = rtls.m_fullMass1 - res.m_mT;
    assert(!IsNeg(propMassSpent));
    curr += sprintf(curr, "%.16e",  propMassSpent.Magnitude());

    // Constraints:
    // XXX: At the moment, we only constrain the Landing Velocity and Landing
    // Absolute Acceleration:
    assert(!IsNeg(res.m_VT));
    curr += sprintf(curr, " %.16e",
                   (res.m_VT - m_proto->m_landVelLimit).Magnitude());

    assert(!IsNeg(res.m_aT));
    curr += sprintf(curr, " %.16e",
                   (res.m_aT - m_proto->m_landAccLimit).Magnitude());

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
// "CalibrateLandBurnParams":                                                //
//===========================================================================//
void RTLS1::CalibrateLandBurnParams(std::string const& a_config_ini)
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

  // "Technical" Params:
  Time   odeIntegrStep          (pt.get<double>("Technical.ODEIntegrStep"));
  int    optMaxEvals      =      pt.get<int>   ("Technical.OptMaxEvals");
  int    optLogLevel      =      pt.get<int>   ("Technical.OptLogLevel");
  int    optSeed          =      pt.get<int>   ("Technical.NOMADSeed");
  bool   stopIfFeasible   =      pt.get<bool>  ("Technical.NOMADStopIfFeasible");
  double useVNS           =      pt.get<double>("Technical.NOMADUseVNS");
  bool   useMT            =      pt.get<bool>  ("Technical.NOMADUseMT");

  if (useVNS < 0.0 || useVNS >= 1.0)      // 0: VNS not used
    throw std::invalid_argument
          ("RTLS1::CalibrateLandBurnParams: NOMADUseVNS: Must be in [0..1)");

  // The Landing Velocity and Acceleration Limits:
  VelK   landVelLimit           (pt.get<double>("Opt.LandVelLimit"));
  AccK   landAccLimit     =      pt.get<double>("Opt.LandAccLimitG") * g0K;

  // For completeness, we need the Q and LongG Limits as well, although they
  // are highly unlikely to be exceeded:
  Pressure QLimit               (pt.get<double>("Opt.QLimit"  ));
  double   longGLimit    =       pt.get<double>("Opt.LongGLimit");

  if (!(IsPos(landVelLimit) && IsPos(landAccLimit) && IsPos(QLimit) &&
        longGLimit > 0.0))
    throw std::invalid_argument
          ("RTLS1::CalibrateLandBurnParams: Invalid Limit(s)");

  //-------------------------------------------------------------------------//
  // State @ "MaxLandBurnH":                                                 //
  //-------------------------------------------------------------------------//
  constexpr VelK      VMin       (0.20);
  constexpr VelK      VMax       (0.20);
  constexpr VelK      dV         (0.01);

  constexpr Angle_deg psiMin     (-120.0);
  constexpr Angle_deg psiMax     ( -90.0);
  constexpr Angle_deg dPsi       (   3.0);

  constexpr Mass      propMMin =  2500.0_kg;
  constexpr Mass      propMMax = 10000.0_kg;
  constexpr Mass      dPropM   =   100.0_kg;

  //-------------------------------------------------------------------------//
  // Grid of State Vals:                                                     //
  //-------------------------------------------------------------------------//
  for (VelK      v     = VMin;     v     <= VMax;     v     += dV)
  for (Angle_deg psi   = psiMin;   psi   <= psiMax;   psi   += dPsi)
  for (Mass      propM = propMMin; propM <= propMMax; propM += dPropM)
  {
    // The "RTLS1" obj:
    // (*) the horizontal distance is not entirely irrelevant (because there are
    //     some built-in constraints on "dL" in "RTLS1"), but setting it to 0 is
    //     OK here;
    // (*) "dVhor" estimate and "landDLLimit" can be safely set to 0;
    // (*) obviously, approxLandBurn = false;
    // (*) optimisation Ranges are irrelevant for the LandBurn:
    //
    double propMassSRange   [2] {NAN, NAN};
    double bbBurnDurRange   [2] {NAN, NAN};
    Time   entryBurnDurRange[2] {Time(NAN), Time(NAN)};

    RTLS1 proto
    (
      maxFPLMass1,    fplK1, fplPropRem1,  IspSL1,   IspVac1,  thrustVacI1,
      minThrtL1,      diam,
      propM,          MaxLandBurnH,  0.0_km,    v,   To_Angle(psi),
      VelK(0.0),      0.0_km,        landVelLimit,   landAccLimit,
      QLimit,         longGLimit,    false,
      propMassSRange, NAN,           bbBurnDurRange, entryBurnDurRange,
      odeIntegrStep,  &std::cout,    optLogLevel
    );
  }
}
}

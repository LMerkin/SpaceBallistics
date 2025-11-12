// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/RTLS1-Landing.cpp":                     //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include "SpaceBallistics/Missions/MkNOMADParams.hpp"

namespace SpaceBallistics
{
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
  LenK UNUSED_PARAM(a_ref_h),
  VelK UNUSED_PARAM(a_ref_Vr),
  VelK UNUSED_PARAM(a_ref_Vhor),
  Mass UNUSED_PARAM(a_ref_prop_mass)
)
{
  // XXX: For the moment, it is a placeholder only: No LandBurn actually occurs,
  // we just hit the Earth surface:
  m_landBurnH     = 0.0_km;
  m_landBurnGamma = 0.0;
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

  // Optimisation Limits (Constraints):
  VelK          const m_landVelLimit;
  AccK          const m_landAccLimit;

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
    Mass                                          a_prop_massL,
    VelK                                          a_land_vel_limit,
    AccK                                          a_land_acc_limit
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_VL            (a_VL),
    m_psiL          (a_psiL),
    m_propMassL     (std::max(a_prop_massL, m_proto->m_unSpendable1)),
    m_landVelLimit  (a_land_vel_limit),
    m_landAccLimit  (a_land_acc_limit)
  {
    assert(m_proto != nullptr);
    if (!(IsPos(m_landVelLimit) && IsPos(m_landAccLimit) && IsPos(m_VL)))
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
      VelK(0.0),    // "dVhor" estimate does not matter here
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
    curr += sprintf(curr, " %.16e", (res.m_VT - m_landVelLimit).Magnitude());

    assert(!IsNeg(res.m_aT));
    curr += sprintf(curr, " %.16e", (res.m_aT - m_landAccLimit).Magnitude());

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
void RTLS1::CalibrateLandBurnParams()
{
}
}

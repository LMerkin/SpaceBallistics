// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/Missions/LVBase.h":                   //
//           Base Class for Ballistic Analysis of Launch Vehicles            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"
#include <utility>
#include <optional>
#include <ostream>

namespace SpaceBallistics
{
  namespace EAM = EarthAtmosphereModel;

  //=========================================================================//
  // "LVBase" Class:                                                         //
  //=========================================================================//
  // This class mostly contains re-usable helpers used by its Derived Classes,
  // and some Stage1 info.   The "Derived" must provide the following methods:
  // (*) "AoA"
  // (*) "PropBurnRate"
  // (*) "Thrust"
  // (*) "LVMass"
  //
  template<typename Derived>
  class LVBase
  {
  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    constexpr static GMK      K              = BodyData<Body::Earth>::K;
    constexpr static LenK     R /* Mean */   = BodyData<Body::Earth>::Rm;

    // ODE Integration Params: 1 msec step; it may only be reduced, never
    // increased beyond the original value:
    constexpr static double   ODERelPrec     = 1e-6;

    // Singular Point detection criteria: NB: "Vhor" approaches 0 much faster
    // than "Vr":
    constexpr static VelK     SingVr         = VelK(1e-2); // 10   m/sec
    constexpr static VelK     SingVhor       = VelK(1e-4); //  0.1 m/sec

  protected:
    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    using MassT2 =   decltype(MassRate(1.0) / 1.0_sec);
    using MassT3 =   decltype(MassT2  (1.0) / 1.0_sec);
    using AreaK  =   decltype(Sqr(1.0_km));

    //-----------------------------------------------------------------------//
    // "StateV":                                                             //
    //-----------------------------------------------------------------------//
    // State Vector (for some time t <= 0, where t=t0=0 corresponds to the Orbi-
    // tal Insertion):
    // Computations are performed in TopoCentric (Start) Planar Polar CoOrds, so
    // the first 3 components are r, rDot, omega = phiDot. The 4th  component is
    // the Total Spent Propellant Mass between the curr "t" and "t0" (which is a
    // continuous and DECREASING function of "t",  as opposed to the Total Mass
    // which is discontinuous when the Stage1 or Fairing are jettisones), and
    // the 5th component is the polar angle "phi" (integrated omega):
    //
    using StateV  = std::tuple<LenK, VelK, AngVel, Mass, Angle>;
    //                         r   rDot=Vr omega   spent phi

    //-----------------------------------------------------------------------//
    // "DStateV":                                                            //
    //-----------------------------------------------------------------------//
    // The Time Derivative of the "StateV". The "MassRate" components is the
    // negated (<= 0) BurnRate:
    //
    using DStateV = std::tuple<VelK, AccK, AngAcc, MassRate, AngVel>;

    //-----------------------------------------------------------------------//
    // "NearSingularityExn" Class:                                           //
    //-----------------------------------------------------------------------//
    // Exception thrown when the flight path is approaching the singular point:
    // the horizontal velocity is assumed to be 0 at that point,  so only "Vr"
    // matters (which is small as well):
    //
    struct NearSingularityExn
    {
    public:
      // Data Flds:
      LenK  const m_r;
      VelK  const m_Vr; // Because "Vhor" is assumed to be 0 near Singularity
      Mass  const m_spentPropMass;
      Angle const m_phi;
      Time  const m_t;

      // Non-Default Ctor:
      NearSingularityExn(LenK  a_r,   VelK a_Vr, Mass a_spent_prop_mass,
                         Angle a_phi, Time a_t)
      : m_r             (a_r),
        m_Vr            (a_Vr),
        m_spentPropMass (a_spent_prop_mass),
        m_phi           (a_phi),
        m_t             (a_t)
      {}
    };

  public:
    //-----------------------------------------------------------------------//
    // "RunRC":                                                              //
    //-----------------------------------------------------------------------//
    // Possible Return Code of Integration "Run":
    //
    enum class RunRC: int
    {
      // The "desirable" result: Singular Point reached: V = 0 at some h >= 0:
      Singularity = 0, 

      // Reached h = 0 at V > 0, still with unspent propellant (beyond the
      // minimum remnant):
      ZeroH       = 1,

      // Ran out of all available propellant (up to the minimum remnant), but
      // still h > 0 and V > 0:
      FlameOut    = 2,

      // If any exception occurred:
      Error       = 3
    };

    static char const* ToString(RunRC a_rc)
    {
      switch (a_rc)
      {
        case RunRC::Singularity: return "Singularity";
        case RunRC::ZeroH      : return "ZeroH";
        case RunRC::FlameOut   : return "FlameOut";
        default                : return "Error";
      }
    }

    //-----------------------------------------------------------------------//
    // "RunRes" Struct:                                                      //
    //-----------------------------------------------------------------------//
    // "RunRC" plus extra info. XXX: It does not contain the Final Altitude, as
    // it is always converted into the equivalent Velocity. However, the final
    // Down-Range distance is returned:
    //
    struct RunRes
    {
      RunRC    const m_rc;       // Return Code
      Time     const m_T;        // Final Flight Time (< 0 if Bwd integration)
      LenK     const m_LT;       // Final Down-Range distance
      VelK     const m_VT;       // Final Velocity
      Mass     const m_mT;       // Final LV Mass
      Pressure const m_maxQ;     // Max Dynamic Pressure (Q) encountered so far
      Pressure const m_sepQ;     // Q @ Stage1 Separation   (if encountered)
      double   const m_maxLongG; // Max Longitudinal G       encountered so far
    };

  protected:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // LV Params:                                                            //
    //-----------------------------------------------------------------------//
    // Stage1 Params: Always present (unlike Stage2 ones which are placed in
    // the Derived class) :
    // Const Stage1 Params:
    AreaK    const   m_crosS;  // XXX: Assume it is constant across all Stages
    double   const   m_K1;
    double   const   m_propRem1;
    Time     const   m_IspSL1;
    Time     const   m_IspVac1;
    double   const   m_minThrtL1;

    // Stage1 Params which may be subject to Optimisation:
    Mass             m_fullMass1;
    Mass             m_emptyMass1;
    Mass             m_propMass1;
    Mass             m_unSpendable1;
    Mass             m_spendable1;
    ForceK           m_thrustVacI1;    // Nominal (@ Max BurnRate)
    MassRate         m_burnRateI1;     // BurnRate @ Stage1 IgnTime
    Time             m_T1;             // Nominal Stage1 BurnTime

    // Integration and Output Params:
    Time             m_odeIntegrStep;  // Typically ~10 msec
    std::ostream*    m_os;
    int              m_logLevel;
    //
    // LogLevel:
    // 0: No output (even if "m_os" is non-NULL)
    // 1: Errors and warnings
    // 2: Results of each run (within the optimisation process)
    // 3: Important  events
    // 4: Full trajectory
    // 5: Debugging info

    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    LVBase
    (
      // Stage1:
      Mass           a_full_mass1,
      double         a_K1,                  // PropMass1 / FullMass1
      double         a_prop_rem1,
      Time           a_Isp_sl1,
      Time           a_Isp_vac1,
      ForceK         a_thrust_vac1,
      double         a_min_thrtl1,
      Len            a_diam,

      // Integration and Output Params:
      Time           a_ode_integr_step,
      std::ostream*  a_os,
      int            a_log_level
    );

    //-----------------------------------------------------------------------//
    // "Mach": The Mach Number:                                              //
    //-----------------------------------------------------------------------//
    static double Mach(EAM::AtmConds const& a_atm, VelK a_v);

    //-----------------------------------------------------------------------//
    // "NonGravForces":                                                      //
    //-----------------------------------------------------------------------//
    void NonGravForces
    (
      // Inputs:
      Time           a_t,
      LenK           a_r,
      VelK           a_Vr,
      VelK           a_Vhor,
      Mass           a_m,
      // Outputs:
      VelK*          a_V,
      Angle*         a_psi,
      Angle*         a_aoa,
      EAM::AtmConds* a_atm,
      MassRate*      a_burn_rate,
      ForceK*        a_thrust,
      double         a_lv_axis[2],  // In the (r, normal-to-r) frame
      AccK           a_ng_acc [2],  // ditto
      double*        a_long_g
    )
    const;

    //-----------------------------------------------------------------------//
    // "ODERHS":                                                             //
    //-----------------------------------------------------------------------//
    // NB: It is generic (can be put in the parent class) since it is implemen-
    // ted in terms of other methods:
    //
    DStateV ODERHS(StateV const& a_s, Time a_t) const;

    //-----------------------------------------------------------------------//
    // Integration Post-Processing:                                          //
    //-----------------------------------------------------------------------//
    // "LocateSingularPoint":
    // If we have arrived in a vicinity of the Singular Point:
    //
    RunRes LocateSingularPoint(NearSingularityExn const& a_nse) const;

    // "PostProcessRun":
    // If integration has come to completion but NOT to the Singular Point:
    //
    RunRes PostProcessRun(StateV const& a_sT, Time a_T) const;

  private:
    //-----------------------------------------------------------------------//
    // Cast into "Derived":                                                  //
    //-----------------------------------------------------------------------//
    Derived*       ToDer()       { return static_cast<Derived*>      (this); }
    Derived const* ToDer() const { return static_cast<Derived const*>(this); }
  };
}
// End namespave SpaceBallistics

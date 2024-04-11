// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CoOrds/CoOrds.h":                    //
//                   Co-Ordinate Systems and State Vectors                   //
//===========================================================================//
#pragma once
#include "SPaceBallistics/Types.hpp"

namespace SpaceBallistics
{
	//=========================================================================//
	// "StateVector" Struct:																									 //
	//=========================================================================//
  template<typename COS>
	struct StateVector
	{
  public:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Mandatory: For the Steady Motion of the CoM (in the specified COS):
    Time          m_t;      // TODO: TBD: May depend on COS
		PosV<COS>     m_r;      // (x,y,z) co-ords or the CoM
    Vel           m_v;      // (x_dot, y_dot, z_dot)  CoM velocity components

    // Mandatory: For the Rotational Motion around the CoM (in that COS):
    // Euler's Angles: (Pitch, Yaw, Roll). NB: They do not make a vector in the
    // mathematical sense; sometimes they are not well-defined  (eg for a snap-
    // shot of the EmbeddedCOS, the current Euler's angles would normally be 0):
    Angle         m_pitch;  // Aka "theta"
    Angle         m_yaw;    // Aka "psi"
    Angle         m_roll;   // Aka "phi"

    // Vector of Angular Velocity. NB: This is NOT the vector of derivatives
    // (Pitch_dot, Yaw_dot, Roll_dot); the latter are given separately below:
    AngVelV<COS>  m_angVel;

    // For info: Euler's Angles Derivatives:
    AngleVel      m_pitchDot;
    AngleVel      m_yawDot;
    AngleVel      m_rollDot;

    // For info: Linear and Angular Accelerations (in that COS). Similar
    // to "m_angVel", "m_angAcc
    AccV   <COS>  m_acc;    // LinearAccels: (x_ddot, y_ddot, z_ddot)
    AngAccV<COS>  m_angAcc; //

    // For info: 2nd Derivatives of of Euler's Angles:
    AngleAcc      m_pitchDDot;
    AngleAcc      m_yawDDot;
    AngleAcc      m_rollDDot;

  private:
    //-----------------------------------------------------------------------//
    // Conversion of Angular Velocities and Accelerations:                   //
    //-----------------------------------------------------------------------//
	};
}

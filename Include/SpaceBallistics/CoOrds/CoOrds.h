// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CoOrds/CoOrds.h":                    //
//                   Co-Ordinate Systems and State Vectors                   //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // Co-Ord Systems:                                                         //
  //=========================================================================//
  //=========================================================================//
  // "EmbeddedCOS" Class:                                                    //
  //=========================================================================//
  // Origin: A fixed physical point               in the LV or SC
  // Axes  : Fixed physical axes (eg of symmetry) of the LV or SC
  // NB    : This type is templatised by the LV/SC Type. It is intended to just
  //         stand for itself (no objs of this class are to be created):
  //
  enum class LVSC: int
  {
    Soyuz21b
    // Others TBD...
  };

  template<LVSC LVSCKind>
  class EmbeddedCOS
  {
    EmbeddedCOS() = delete;
  };

  //=========================================================================//
  // "TopoCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin: A point on the Earth surface given by "L"
  // Axes  : (X=East, Y=North, Z=Zenith)
  // NB    : Again, no objects of this class are to be created:
  //
  template<Location_WGS84 const* L>
  class TopoCentricCOS
  {
    TopoCentricCOS() = delete;
  };

  //=========================================================================//
  // "GeoCRotatingCOS" Class:                                                //
  //=========================================================================//
  class GeoCRotatingCOS
  {
  private:
    Time      m_epoch;  // E.g. 2024.0
  };

	//=========================================================================//
	// "StateVector" Struct:																									 //
	//=========================================================================//
  template<typename COS>
	struct StateVector
	{
    // Mandatory: For the Steady Motion of the CoM (in the specified COS):
    Time    m_t;          // TODO: TBD: May depend on COS
		Len     m_pos   [3];  // (x,y,z) co-ords or the CoM
    Vel     m_vel   [3];  // (x_dot, y_dot, z_dot)  CoM velocity components

    // Mandatory: For the Rotational Motion around the CoM (in that COS):
    double  m_angles[3];  // Euler's Angles:     (Pitch,    Yaw,     Roll)
    AngVel  m_angVel[3];  // Angular Velocities: (Pith_dot, Yaw_dot, Roll_dot)

    // Accelerations (for convenience) (in that COS):
    Acc     m_acc   [3];  // LinearAccels: (x_ddot,     y_ddot,   z_ddot)
    AngAcc  m_angAcc[3];  // AngAccels   : (Pitch_ddot, Yaw_ddot, Roll_ddot)

    // For convenience, we include the "MechElement" characterising this body;
    // Here CoM, MoIs and their "Dots" are converted from the
	};
}

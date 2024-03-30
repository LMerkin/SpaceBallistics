// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CoOrds/CoOrds.h":                    //
//                   Co-Ordinate Systems and State Vectors                   //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
/*
  //=========================================================================//
  // Co-Ord System Kinds:                                                    //
  //=========================================================================//
  enum class CoOrdSystemK: int
  {
    Embedded,     // Origin: Fixed Body Pt; Axes: Symmetry Axes
    TopoC,        // Origin: Earth Site;    Axes: (X=East, Y=North, Z=Zenith)
    GeoCRotating, // Origin: Earth Center;  Axes: Geographically-Fixed
    GeoCFixed,    // Origin: Earth Center;  Axes: ICRF 2000.0
    BaryC         // Origin: Solar System BaryCenter; Axes: ICRF 2000.0
  };
*/

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
  enum class LVCS: int
  {
    Soyuz21b
    // Others TBD...
  };

  template<LVCS Type>
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
    // Mandatory: For the Steady Motion of the CoM:
    Time    m_t;            // TODO: TBD: May depend on COS
		Len     m_pos   [3];    // (x,y,z) co-ords or the CoM
    Vel     m_vel   [3];    // (x_dot, y_dot, z_dot)  CoM velocity components

    // Mandatory: For the Rotational Motion around the CoM:
    double  m_angles[3];    // Euler's Angles: (Pitch, Yaw, Roll)
    AngVel  m_angVel[3];    // Euler's Angular Velocities

    // Auxiliary state components (for convenience):
    Acc     m_acc   [3];   // Accelertations wrt           (x,y,z)
    Force   m_force [3];   // Equiv.Force's components wrt (x,y,z)
    Mass    m_mass;        // Curr object's (LV/SpaceCraft) Mass

    AngAcc  m_angAcc[3];   // Angular Accelerations
    MoI     m_mois  [3];   // Moments of Inertia wrt       (x,y,z)
    // TODO: Add Moments of Forces wrt CoM (???)
	};
}

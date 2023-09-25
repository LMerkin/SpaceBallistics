// vim:ts=2:et
//===========================================================================//
//                                "Types.hpp":                               //
//        Declarations of Principal Types used in "SpaceBallistics"          //
//===========================================================================//
#pragma  once
#include <DimTypes/DimTypes.hpp>

namespace SpaceBallistics
{
  using namespace DimTypes;

	//=========================================================================//
	// Basic Dimension Types:                                                  //
	//=========================================================================//
	DECLARE_DIMS
	(
		double,
		(Len,  m, (km, 1000.0)),
		(Time, sec),
		(Mass, kg)
	)

  //=========================================================================//
  // Derived Dimension Types and Values:                                     //
  //=========================================================================//
  // Area and Volume:
  using Area       = decltype(Sqr    (Len_m_1));
  using Vol        = decltype(IPow<3>(Len_m_1));

  // Density (ie Volume Density) and Surface (Area) Density:
  using Density    = decltype(Mass_kg_1 / Vol (1.0));
  using SurfDens   = decltype(Mass_kg_1 / Area(1.0));

  // Moment of Inertia:
  using MoI        = decltype(Mass_kg_1 * Area(1.0));

  // Velocity and Acceleration (in SI units):
  using Vel        = decltype(Len_m_1   / Time_sec_1);
  using Acc        = decltype(Len_m_1   / Area(1.0));

  // Standard Gravity:
  constexpr Acc g0 = Acc(9.80665);
}
// End namespace SpaceBallistics

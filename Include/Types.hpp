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
  // Derived Dims and Units:                                                 //
  //=========================================================================//
  // Moment of Inertia:
  using MoI = decltype(Mass_kg_1 * Sqr(Len_m_1));
}
// End namespace SpaceBallistics

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
	// Mechanical Dims and Units:                                              //
	//=========================================================================//
	DECLARE_DIMS
	(
		double,
		(Len,  m, (km, 1000.0)),
		(Time, sec),
		(Mass, kg)
	)
}
// End namespace SpaceBallistics

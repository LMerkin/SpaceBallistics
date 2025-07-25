// vim:ts=2
//===========================================================================//
//                     "Tests/EarthAtmosphereTest.cpp":                      //
//===========================================================================//
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"
#include <iostream>

int main()
{
	namespace EAM = SpaceBallistics::EarthAtmosphereModel;
	using namespace SpaceBallistics;
	using namespace std;

	for (int i = 0; i <= 6; ++i)
	{
		EAM::LayerInfo const& l = EAM::Layers[i];

		// Output the info at the "base" of the Layer:
		cout << l.m_baseH << '\t' << l.m_baseP << '\t' << l.m_baseT << endl;

		if (i == 6)
			// This is the last Layer, so output the info at its upper boundary
			// as well:
			cout << l.m_endH << '\t' << l.P(l.m_endH) << '\t' << l.T(l.m_endH)
					 << endl;
	}
	return 0;
}

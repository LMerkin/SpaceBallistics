// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/PhysForces/DE440T.h":                   //
//                         JPL DE440T Ephemerides                            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/TimeScales.h"

namespace SpaceBallistics
{
	//=========================================================================//
	// "DE440T" Class:                                                         //
	//=========================================================================//
	class DE440T
	{
	private:
		//-----------------------------------------------------------------------//
		// Consts:																														   //
		//-----------------------------------------------------------------------//
		// The Temporal Range of Records (in TDB):
		constexpr static TDB From = TDB(Time_day(2'323'696.5)); // 1649-12-18.0
    constexpr static TDB To   = TDB(Time_day(2'506'320.5)); // 2149-12-21.0

    // The Number of Records:
    constexpr static int NR   = 5707;

    // Temporal Span of each Record:
    constexpr static Time_day   RecSpan =  Time_day(32.0);
    static_assert(To_Time(double(NR) * RecSpan) == To - From);

    // Size of each record in "double"s:
    constexpr static int ND  = 1122;

    //-----------------------------------------------------------------------//
    // Data Record Layout:                                                   //
    //-----------------------------------------------------------------------//
    struct Record
    {
      // Time Span of this Record, JD_TDB:
      Time_day    m_From;
      Time_day    m_To;

      // Chebyshev Coeffs for 3D Planetary Positions, in the "BaryCentricCOS";
      // the last dim is [X, Y, Z]:
      Len_km      m_Mercury       [4][14][3];
      Len_km      m_Venus         [2][10][3];
      Len_km      m_EMB           [2][13][3];  // Earth-Moon System BaryCenter
      Len_km      m_Mars             [11][3];
      Len_km      m_Jupiter           [8][3];
      Len_km      m_Saturn            [7][3];
      Len_km      m_Uranus            [6][3];
      Len_km      m_Neptune           [6][3];
      Len_km      m_Pluto             [6][3];

      // Chebyshev Coeffs for the 3D Moon Position, in the "GeoCentricFixedCOS";
      // the last dim is  [X, Y, Z]:
      Len_km      m_MoonGeoC      [8][13][3];

      // Chebyshev Coeffs for the Position of the Sun, in the "BaryCentricCOS":
      Len_km      m_Sun           [2][11][3];

      // Chebyshev Coeffs for Long-Period Earth Nutations:
      // the last dim is [d(psi), d(eps)]:
      Angle       m_EarthNutation [4][10][2];

      // Chebyshev Coeffs for Lunar Mantle (Exterior) Librations:
      // the last dim is [phi, theta, psi]:
      Angle       m_MoonLibration [4][10][3];

      // Chebyshev Coeffs for (TT-TDB):
      Time        m_TTmTDB        [8][13];
    };
    // The size of the above "Record" is the same as that of "ND" "double"s:
    static_assert(sizeof(Record) == size_t(ND) * sizeof(double));

    //-----------------------------------------------------------------------//
    // The Actual Data:                                                      //
    //-----------------------------------------------------------------------//
    static double const s_data[NR][ND];
	};
}
// End namespace SpaceBallistics

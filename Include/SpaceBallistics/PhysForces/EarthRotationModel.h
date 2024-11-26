// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/PhysForces/EarthRotationModel.h":            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/TimeScales.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "EarthRotationModel" Class:                                             //
  //=========================================================================//
  // An object of this class  describes  the orientation  of the Earth rotation
  // axis (CIP) and the RA origin (CEO), ie the CRS, for a given Time (Epoch).
  // The orientatation is represented as a certain rotation matrix  (expressing
  // GCRS co-ords via the CRS-of-Epoch co-ords) and changes very slowly with the
  // Epoch, due to Precession and Nutation. If the Epoch is J2000.0, the matrix
  // is extremely close to the Unit Matrix I (but not equal to I exactly because  // GCRS is slightly different from CRS_J2000.0). 
  // The "EarthRotationModel" object constructed then acts as a function trans-
  // forming GTRS co-ords (for a given TT) to GCRS.
  // XXX: For Nutations, an analytical model (integrated with the Precession mo-
  // del) is used rather than DE440T numerical data. This analytical model uses
  // power series expansions around J2000.0, and is only valid for a few centu-
  // ries  --  roughly for the same period as provided by our implementation of
  // DE440T:
  //
  class EarthRotationModel
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    double  m_M[3][3];      // r_GCRS = M * r_[CRS_Epoch]

    // Default Ctor is deleted: without the Epoch, this model makes no sense:
    EarthRotationModel() = delete;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Constructing the model for a given Year (may be fractional); the exact
    // TimeScale is presumably TT, but that does not matter here  because the
    // result depends on the Epoch VERY WEAKLY:
    //
    constexpr EarthRotationModel(Time_jyr a_epoch)
    {
      double T = (a_epoch - 2000.0_jyr).Magnitude() / 100.0;
    }
  };
}
// End namespace SpaceBallistics

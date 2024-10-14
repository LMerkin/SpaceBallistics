// vim:ts=2:et
//===========================================================================//
//                        "Src/ME/TanksPressrn.hpp":                         //
//            Modeling of Propellant Tanks Pressurisation Effects            //
//===========================================================================//
#include "SpaceBallistics/ME/MechElement.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TanksPressrn":                                                         //
  //=========================================================================//
  // Constructs an "ME" which is a result of filling in  the volumes  above the
  // (decreasing) Propellant Bulks in the Propellant Tanks, with Pressurisation
  // Gas, using the "method of complements":
  //
  template<LVSC LVSCKind>
  MechElement<LVSCKind> TanksPressrn
  (
    MechElement<LVSCKind> const& a_fuel_bulk_full,
    MechElement<LVSCKind> const& a_fuel_bulk_curr,
    MechElement<LVSCKind> const& a_oxid_bulk_full,
    MechElement<LVSCKind> const& a_oxid_bulk_curr,
    Mass                         a_gas_mass,
    MassRate                     a_gas_mass_dot   // Gas Generation/Supply Rate
  )
  {
    assert(IsPos(a_gas_mass) && !IsNeg(a_gas_mass_dot));

    // First, determine the volumes free of Fuel and Oxid -- they will be filled
    // with the Gas. The simplest way to do so is to compute the "Complementary"
    // MEs. Then obtain the volume, density and volume rate of "compl*ME"s -- as
    // yet, assuming that they are filled by the Propellants (of constant densi-
    // ty):
    // Fuel:
    auto compFuelME      = a_fuel_bulk_full - a_fuel_bulk_curr;
    auto compFuelVol     = compFuelME.GetEnclVol();
    auto compFuelMass    = compFuelME.GetMass   ();
    auto compFuelMassDot = compFuelME.GetMassDot();
    auto fuelDens        = compFuelMass    / compFuelVol;
    assert(IsPos(compFuelVol)     && IsPos(compFuelMass) &&
          !IsNeg(compFuelMassDot) && IsPos(fuelDens));
    auto compFuelVolDot  = compFuelMassDot / fuelDens;

    // Oxid (similar):
    auto compOxidME      = a_oxid_bulk_full - a_oxid_bulk_curr;
    auto compOxidVol     = compOxidME.GetEnclVol();
    auto compOxidMass    = compOxidME.GetMass   ();
    auto compOxidMassDot = compOxidME.GetMassDot();
    auto oxidDens        = compOxidMass    / compOxidVol;
    assert(IsPos(compOxidVol)     && IsPos(compOxidMass) &&
          !IsNeg(compOxidMassDot) && IsPos(oxidDens));
    auto compOxidVolDot  = compOxidMassDot / oxidDens;

    // The total volume to be filled by the Pressurisation Gas:
    auto gasVol          = compFuelVol    + compOxidVol;
    auto gasVolDot       = compFuelVolDot + compOxidVolDot;

    // The current gas density (XXX: however, in contrast to the Propellants,
    // it is variable over the time):
    auto gasDens         = a_gas_mass / gasVol;

    // Pro-Rating Coeffs (Propellants -> Gas):
    double compFuelPRC   = double(gasDens / fuelDens);
    double compOxidPRC   = double(gasDens / oxidDens);

    // However, in this case, "ProRateMass" must also take into account the
    // GasDensityDot = d/dt (GasMass / GasVol):
    auto gasDensDot      =
      a_gas_mass_dot / gasVol - gasDens / gasVol * gasVolDot;

    // We can now apply "ProRateMass":
    auto gasFME =
      MechElement<LVSCKind>::template ProRateMass<true>
        (compFuelME, compFuelPRC, gasDensDot);
    auto gasOME =
      MechElement<LVSCKind>::template ProRateMass<true>
        (compOxidME, compOxidPRC, gasDensDot);

    // The final result is the sum for both Tanks:
    return gasFME + gasOME;
  }
}
// End namespace SpaceBallistics

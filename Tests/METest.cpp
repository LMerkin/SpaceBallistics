// vim:ts=2:et
//===========================================================================//
//                            "Tests/METest.cpp":                            //
//                      Tests of the "MechElement" Class                     //
//===========================================================================//
#include "SpaceBallistics/ME/TrConeSpherSegm.hpp"
#include "SpaceBallistics/ME/ToricSegms.hpp"
#include <iostream>

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace SpaceBallistics;
  using namespace std;
  using ME  = MechElement   <LVSC::BoilerPlate0>;
  using TC  = TrCone        <LVSC::BoilerPlate0>;
  using SpS = SpherSegm     <LVSC::BoilerPlate0>;
  using TS  = ToricSegm     <LVSC::BoilerPlate0>;
  using DC  = DoubleCylinder<LVSC::BoilerPlate0>;

  // Used for both Cylinder and Spehere:
  constexpr Len       D   = 2.0_m;
  constexpr Len       HCyl= 5.0_m;
  constexpr Density   rho  (1.0);
  constexpr SurfDens  sigma(1.0);
  constexpr Len       R   = D / 2.0;
  constexpr Area      R2  = Sqr(R);
  constexpr Vol       R3  = R * R2;
  // Major Diameter used for Torus:
  constexpr Len       TD  = 3.0 * D;

  //=========================================================================//
  // MoIs and CoM of a Cylinder:                                             //
  //=========================================================================//
  constexpr Vol       vCyl        = Pi<double> * R2 * HCyl;
  constexpr Area      sCyl        = Pi<double> * D  * HCyl;
  constexpr Mass      propMassCyl = vCyl * rho;
  constexpr Mass      surfMassCyl = sCyl * sigma;

  // Upper Base is @ X=0; therefore, the center is @ (-HCyl/2):
  constexpr Len       XCCyl       = -HCyl / 2.0;
  constexpr TC        cyl(0.0_m, D,  HCyl, rho, surfMassCyl);

  // Computed CoM and MoIs for the Solid Cylinder (3D):
  constexpr   ME propCyl = cyl.GetPropBulkME(propMassCyl);
  auto const& comCyl     = propCyl.GetCoM ();
  auto const& moisCyl    = propCyl.GetMoIs();

  // Diffs with Expected CoM:
  Len dCCyl  [3] {Abs(comCyl[0] - XCCyl), Abs(comCyl[1]), Abs(comCyl[2])};

  // Expected MoIs for the Cylinder:
  constexpr MoI JxCyl  = propMassCyl *  Sqr(R) / 2.0;
  constexpr MoI JyzCyl = propMassCyl * (Sqr(R) / 4.0 + Sqr(HCyl) / 3.0);

  MoI dJCyl  [3]
    { Abs(moisCyl[0] - JxCyl ),
      Abs(moisCyl[1] - JyzCyl),
      Abs(moisCyl[2] - JyzCyl) };

  cout << "SOLID CYLINDER:"  << endl
       << "dxC=" << dCCyl[0] << "\ndyC=" << dCCyl[1] << "\ndzC=" << dCCyl[2]
       << endl
       << "dJx=" << dJCyl[0] << "\ndJy=" << dJCyl[1] << "\ndJz=" << dJCyl[2]
       << endl   << endl;

  // And now for the Empty (Hollow) Cylinder (2D):
  constexpr MoI JxECyl  = surfMassCyl * R2;
  constexpr MoI JyzECyl = surfMassCyl * (R2 / 2.0 + Sqr(HCyl) / 3.0);

  auto const& comECyl  = cyl.GetCoM ();
  auto const& moisECyl = cyl.GetMoIs();

  Len dMECyl [3]
    { Abs(comECyl[0] - XCCyl), Abs(comECyl[1]), Abs(comECyl[2]) };

  MoI dJECyl [3]
    { Abs(moisECyl[0] - JxECyl ),
      Abs(moisECyl[1] - JyzECyl),
      Abs(moisECyl[2] - JyzECyl) };

  cout << "HOLLOW CYLINDER (SHELL):"      << endl
       << "dxC=" << dMECyl[0] << "\ndyC=" << dMECyl[1] << "\ndzC=" << dMECyl[2]
       << endl
       << "dJx=" << dJECyl[0] << "\ndJy=" << dJECyl[1] << "\ndJz=" << dJECyl[2]
       << endl   << endl;

  //=========================================================================//
  // MoIs and CoM of a Sphere:                                               //
  //=========================================================================//
  // It is made of 2 HemiSpeheres: Facing Up and Low: Center at x=R:
  //
  constexpr Vol  vHS         = (2.0/3.0) * Pi<double> * R2 * R;
  constexpr Mass propMassHS  = vHS * rho;               // HemiSphere
  constexpr Mass propMassS   = 2.0 * propMassHS;        // FullSphere
  constexpr Mass emptyMassHS = TwoPi<double> * R2 * sigma;
  constexpr Mass emptyMassS  = 2.0 * emptyMassHS;

  constexpr SpS  upHS (true,  R, D, rho, emptyMassHS);
  constexpr SpS  lowHS(false, R, D, rho, emptyMassHS);

  constexpr ME   upPropME    = upHS .GetPropBulkME(propMassHS);
  constexpr ME   lowPropME   = lowHS.GetPropBulkME(propMassHS);

  auto const&    upCoM       = upPropME .GetCoM ();
  auto const&    upMoIs      = upPropME .GetMoIs();
  auto const&    lowCoM      = lowPropME.GetCoM ();
  auto const&    lowMoIs     = lowPropME.GetMoIs();

  // Over-All CoM and MoIs for the Solid Sphere (3D):
  Len comS[3]
  {
    0.5 * (upCoM[0] + lowCoM[0]),
    0.5 * (upCoM[1] + lowCoM[1]),
    0.5 * (upCoM[2] + lowCoM[2])
  };
  MoI moisS[3]
  {
    upMoIs[0] + lowMoIs[0],
    upMoIs[1] + lowMoIs[1],
    upMoIs[2] + lowMoIs[2]
  };

  // Diffs with theoretical vals:
  Len dCS [3] { Abs(comS[0] - R), Abs(comS[1]), Abs(comS[2]) };

  constexpr MoI JxS  = 0.4 * propMassS * R2;
  constexpr MoI JyzS = JxS + propMassS * R2;  // Steiner's Parallel Axis Theorem
  MoI dJS [3]
    { Abs(moisS[0] - JxS), Abs(moisS[1] - JyzS), Abs(moisS[2] - JyzS) };

  cout << "SOLID SPHERE:"    << endl
       << "dxC="   << dCS[0] << "\ndyC=" << dCS[1] << "\ndzC=" << dCS[2]
       << endl
       <<   "dJx=" << dJS[0] << "\ndJy=" << dJS[1] << "\ndJz=" << dJS[2]
       << endl     << endl;

  // Finally, the Empty (Hollow) Sphere (2D):
  auto const& upECoM   = upHS .GetCoM ();
  auto const& upEMoIs  = upHS .GetMoIs();
  auto const& lowECoM  = lowHS.GetCoM ();
  auto const& lowEMoIs = lowHS.GetMoIs();

  Len  comES [3]
  {
    0.5 * (upECoM[0] + lowECoM[0]),
    0.5 * (upECoM[1] + lowECoM[1]),
    0.5 * (upECoM[2] + lowECoM[2])
  };
  MoI  moisES[3]
  {
    upEMoIs[0] + lowEMoIs[0],
    upEMoIs[1] + lowEMoIs[1],
    upEMoIs[2] + lowEMoIs[2]
  };

  // Diffs with theoretical vals:
  Len dMES[3] { Abs(comES[0] - R), Abs(comES[1]), Abs(comES[2]) };

  constexpr MoI JxES  = 2.0/3.0 * emptyMassS * R2;
  constexpr MoI JyzES = JxES    + emptyMassS * R2;  // Steiner's Theorem...
  MoI dJES[3]
    { Abs(moisES[0] - JxES), Abs(moisES[1] - JyzES), Abs(moisES[2] - JyzES) };

  cout << "HOLLOW SPHERE (SHELL):"        << endl
       << "dxC="   << dMES[0] << "\ndyC=" << dMES[1] << "\ndzC=" << dMES[2]
       << endl
       <<   "dJx=" << dJES[0] << "\ndJy=" << dJES[1] << "\ndJz=" << dJES[2]
       << endl     << endl;

  //=========================================================================//
  // MoIs and CoM of a Torus:                                                //
  //=========================================================================//
  // It is made of 2 ToricSegments: Facing Up and Low: Center at x=R:
  //
  constexpr Len   Q           = (TD - D)  / 2.0;
  constexpr Vol   vHT         = Sqr(      Pi<double>  * R) * Q;
  constexpr Area  sHT         = 2.0 * Sqr(Pi<double>) * R  * Q;
  constexpr Mass  propMassHT  = vHT * rho;               // ToricSegm
  constexpr Mass  propMassT   = 2.0 * propMassHT;        // FullTorus
  constexpr Mass  emptyMassHT = sHT * sigma;
  constexpr Mass  emptyMassT  = 2.0 * emptyMassHT;

  constexpr TS upHT (true,  R, D, TD, rho, emptyMassHT);
  constexpr TS lowHT(false, R, D, TD, rho, emptyMassHT);

  constexpr ME upHTPropME  = upHT .GetPropBulkME(propMassHT);
  constexpr ME lowHTPropME = lowHT.GetPropBulkME(propMassHT);

  auto const&  upHTCoM     = upHTPropME .GetCoM ();
  auto const&  upHTMoIs    = upHTPropME .GetMoIs();
  auto const&  lowHTCoM    = lowHTPropME.GetCoM ();
  auto const&  lowHTMoIs   = lowHTPropME.GetMoIs();

  // Over-All CoM and MoIs for the Solid Torus (3D):
  Len comT[3]
  {
    0.5 * (upHTCoM[0] + lowHTCoM[0]),
    0.5 * (upHTCoM[1] + lowHTCoM[1]),
    0.5 * (upHTCoM[2] + lowHTCoM[2])
  };
  MoI moisT[3]
  {
    upHTMoIs[0] + lowHTMoIs[0],
    upHTMoIs[1] + lowHTMoIs[1],
    upHTMoIs[2] + lowHTMoIs[2]
  };

  // Diffs with theoretical vals:
  Len dCT [3] { Abs(comT[0] - R), Abs(comT[1]), Abs(comT[2]) };

  constexpr MoI JxT  =
    0.25  * propMassT * (4.0 * Sqr(Q) + 3.0 * Sqr(R));
  constexpr MoI JyzT =
    0.125 * propMassT * (4.0 * Sqr(Q) + 5.0 * Sqr(R)) +
    propMassT * R2;     // Steiner's Parallel Axis Theorem

  MoI dJT [3]
    { Abs(moisT[0] - JxT), Abs(moisT[1] - JyzT), Abs(moisT[2] - JyzT) };

  cout << "SOLID TORUS:"     << endl
       << "dxC="   << dCT[0] << "\ndyC=" << dCT[1] << "\ndzC=" << dCT[2]
       << endl
       <<   "dJx=" << dJT[0] << "\ndJy=" << dJT[1] << "\ndJz=" << dJT[2]
       << endl     << endl;

  // Finally, the Empty (Hollow) Torus (2D):
  auto const& upEHTCoM   = upHT .GetCoM ();
  auto const& upEHTMoIs  = upHT .GetMoIs();
  auto const& lowEHTCoM  = lowHT.GetCoM ();
  auto const& lowEHTMoIs = lowHT.GetMoIs();

  Len  comET [3]
  {
    0.5 * (upEHTCoM[0] + lowEHTCoM[0]),
    0.5 * (upEHTCoM[1] + lowEHTCoM[1]),
    0.5 * (upEHTCoM[2] + lowEHTCoM[2])
  };
  MoI  moisET[3]
  {
    upEHTMoIs[0] + lowEHTMoIs[0],
    upEHTMoIs[1] + lowEHTMoIs[1],
    upEHTMoIs[2] + lowEHTMoIs[2]
  };

  // Diffs with theoretical vals:
  Len dMET[3] { Abs(comET[0] - R), Abs(comET[1]), Abs(comET[2]) };

  // Theoretical vals obtained by manual integration:
  constexpr double q2    = Sqr(double(Q / R));
  constexpr MoI    JxET  =
    2.0 * Sqr(Pi<double>) * Q * R3 * (2.0 * q2 + 3.0)  * sigma;
  constexpr MoI    JyzET =
    JxET / 2.0 +
    2.0 * Sqr(Pi<double>) * Q * R3 * sigma  +
    emptyMassT * R2;      // Steiner's Theorem again...

  MoI dJET[3]
    { Abs(moisET[0] - JxET), Abs(moisET[1] - JyzET), Abs(moisET[2] - JyzET) };

  cout << "HOLLOW TORUS (SHELL):"         << endl
       << "dxC="   << dMET[0] << "\ndyC=" << dMET[1] << "\ndzC=" << dMET[2]
       << endl
       <<   "dJx=" << dJET[0] << "\ndJy=" << dJET[1] << "\ndJz=" << dJET[2]
       << endl     << endl;

  //=========================================================================//
  // MoIs and CoM of a Double Cylinder:                                      //
  //=========================================================================//
  // We will use "R" for the Inner Radius and the above-defined "Q" for the
  // outer one, with UpperBase @ x=R:
  //
  constexpr Vol   vDC          = Pi<double> * (Sqr(Q) - Sqr(R)) * HCyl;
  constexpr Area  sDC          = (2.0 * Pi<double>)   * (Q + R) * HCyl;
  constexpr Mass  propMassDC   = vDC  * rho;
  constexpr Mass  emptyMassDC  = sDC  * sigma;

  // Upper Base is @ X=R; therefore, the center is @ (R-HCyl/2):
  constexpr Len   XCDC         = R - HCyl / 2.0;
  constexpr DC    dc(R, 2.0*Q, 2.0*R, HCyl, rho, emptyMassDC);

  constexpr ME    dcPropME     = dc.GetPropBulkME(propMassDC);
  auto const&     dcCoM        = dcPropME.GetCoM ();
  auto const&     dcMoIs       = dcPropME.GetMoIs();

  // Diffs with theoretical vals:
  Len dcDC [3] { Abs(dcCoM[0] - XCDC), Abs(dcCoM[1]), Abs(dcCoM[2]) };

  constexpr MoI JxDC   =
    0.5   * propMassDC * (Sqr(Q) + Sqr(R));
  constexpr MoI JyzDC  =
    propMassDC * (3.0  * (Sqr(Q) + Sqr(R)) + Sqr(HCyl)) / 12.0 +
    propMassDC * Sqr(XCDC);        // Steiner's Parallel Axis Theorem

  MoI dJDC [3]
    { Abs(dcMoIs[0] - JxDC), Abs(dcMoIs[1] - JyzDC), Abs(dcMoIs[2] - JyzDC) };

  cout << "SOLID DOUBLE CYLINDER (THICK TUBE):"      << endl
       << "dxC="   << dcDC[0] << "\ndyC=" << dcDC[1] << "\ndzC=" << dcDC[2]
       << endl
       <<   "dJx=" << dJDC[0] << "\ndJy=" << dJDC[1] << "\ndJz=" << dJDC[2]
       << endl     << endl;

  // Finally, the Empty (Hollow) Double Cylinder (2D):
  auto const& dcECoM    = dc.GetCoM ();
  auto const& dcEMoIs   = dc.GetMoIs();

  // Diffs with theoretical vals:
  Len dcEDC[3] { Abs(dcECoM[0] - XCDC), Abs(dcECoM[1]), Abs(dcECoM[2]) };

  constexpr auto Q3     = Sqr(Q) * Q;
  constexpr MoI  JxEDC  = sigma  * (TwoPi<double> * (R3 + Q3) * HCyl);
  constexpr MoI  JyzEDC =
    TwoPi<double> * HCyl * sigma *
    (0.5 * (R3 + Q3) + (R + Q) * (Sqr(HCyl) / 12.0 + Sqr(XCDC)));

  MoI dJEDC[3]
    { Abs(dcEMoIs[0] - JxEDC),
      Abs(dcEMoIs[1] - JyzEDC),
      Abs(dcEMoIs[2] - JyzEDC) };

  cout << "HOLLOW DOUBLE CYLINDER (2-WALLED SHELL):" << endl
       << "dxC="   << dcEDC[0] << "\ndyC=" << dcEDC[1] << "\ndzC=" << dcEDC[2]
       << endl
       <<   "dJx=" << dJEDC[0] << "\ndJy=" << dJEDC[1] << "\ndJz=" << dJEDC[2]
       << endl     << endl;
  return 0;
}

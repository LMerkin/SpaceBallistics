// vim:ts=2:et
//===========================================================================//
//                            "Tests/CETest.cpp":                            //
//                   Tests of the "ConstrElement" Class                      //
//===========================================================================//
#include "SpaceBallistics/ConstrElement.hpp"
#include <iostream>

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  // Used for both Cylinder amd Spehere:
  constexpr Len       D   = 2.0_m;
  constexpr Density   rho  (1.0);
  constexpr SurfDens  sigma(1.0);
  constexpr Len       R   = D / 2.0;
  constexpr Area      R2  = Sqr(R);

  //-------------------------------------------------------------------------//
  // MoIs and CoM of a Cylinder:                                             //
  //-------------------------------------------------------------------------//
  constexpr Len       hCyl        = 5.0_m;
  constexpr Vol       vCyl        = Pi<double> * R2 * hCyl;
  constexpr Area      sCyl        = Pi<double> * D  * hCyl;
  constexpr Mass      propMassCyl = vCyl * rho;
  constexpr Mass      surfMassCyl = sCyl * sigma;

  constexpr TrCone cyl(0.0_m, D, D, hCyl, rho, surfMassCyl);

  // Computed CoM and MoIs for the Prop Cylinder (3D, w/o the Shell):
  constexpr   ConstrElement propCyl = cyl.GetPropCE(propMassCyl);
  auto const& comCyl      = propCyl.GetCoM ();
  auto const& moisCyl     = propCyl.GetMoIs();

  // Diffs with Expected CoM:
  Len dCCyl  [3] {Abs(comCyl[0] - hCyl/2.0), Abs(comCyl[1]), Abs(comCyl[2])};

  // Expected MoIs for the Cylinder:
  constexpr MoI JxCyl  = propMassCyl *  Sqr(R) / 2.0;
  constexpr MoI JyzCyl = propMassCyl * (Sqr(R) / 4.0 + Sqr(hCyl) / 3.0);

  MoI dJCyl  [3]
    { Abs(moisCyl[0] - JxCyl ),
      Abs(moisCyl[1] - JyzCyl),
      Abs(moisCyl[2] - JyzCyl) };

  cout << "FILLED CYLINDER (w/o Shell):" << endl
       << "dxC=" << dCCyl[0] << "\ndyC=" << dCCyl[1] << "\ndzC=" << dCCyl[2]
       << endl
       << "dJx=" << dJCyl[0] << "\ndJy=" << dJCyl[1] << "\ndJz=" << dJCyl[2]
       << endl   << endl;

  // And now for the Empty Cylinder (2D):
  constexpr MoI JxECyl  = surfMassCyl * R2;
  constexpr MoI JyzECyl = surfMassCyl * (R2 / 2.0 + Sqr(hCyl) / 3.0);

  auto const& comECyl  = cyl.GetCoM ();
  auto const& moisECyl = cyl.GetMoIs();

  Len dCECyl [3]
    { Abs(comECyl[0] - hCyl/2.0), Abs(comECyl[1]), Abs(comECyl[2]) };

  MoI dJECyl [3]
    { Abs(moisECyl[0] - JxECyl ),
      Abs(moisECyl[1] - JyzECyl),
      Abs(moisECyl[2] - JyzECyl) };

  cout << "EMPTY CYLINDER (SHELL):"       << endl
       << "dxC=" << dCECyl[0] << "\ndyC=" << dCECyl[1] << "\ndzC=" << dCECyl[2]
       << endl
       << "dJx=" << dJECyl[0] << "\ndJy=" << dJECyl[1] << "\ndJz=" << dJECyl[2]
       << endl   << endl;

  //-------------------------------------------------------------------------//
  // MoIs and CoM of a Sphere:                                               //
  //-------------------------------------------------------------------------//
  // It is made of 2 HemiSpeheres: Facing left and right: Center at x=D/2:
  //
  constexpr Vol       vHS         = (2.0/3.0) * Pi<double> * R2 * R;
  constexpr Mass      propMassHS  = vHS * rho;               // HemiSphere
  constexpr Mass      propMassS   = 2.0 * propMassHS;        // FullSphere
  constexpr Mass      emptyMassHS = TwoPi<double> * R2 * sigma;
  constexpr Mass      emptyMassS  = 2.0 * emptyMassHS;

  constexpr SpherSegm leftHS (true,  R, D, rho, emptyMassHS);
  constexpr SpherSegm rightHS(false, R, D, rho, emptyMassHS);

  constexpr ConstrElement leftPropCE  = leftHS .GetPropCE(propMassHS);
  constexpr ConstrElement rightPropCE = rightHS.GetPropCE(propMassHS);

  auto const& leftCoM   = leftPropCE .GetCoM();
  auto const& rightCoM  = rightPropCE.GetCoM();
  auto const& leftMoIs  = leftPropCE .GetMoIs();
  auto const& rightMoIs = rightPropCE.GetMoIs();

  // Over-All CoM and MoIs for the Spehere (3D, w/o the Shell):
  Len comS[3]
  {
    0.5 * (leftCoM[0] + rightCoM[0]),
    0.5 * (leftCoM[1] + rightCoM[1]),
    0.5 * (leftCoM[2] + rightCoM[2])
  };
  MoI moisS[3]
  {
    leftMoIs[0] + rightMoIs[0],
    leftMoIs[1] + rightMoIs[1],
    leftMoIs[2] + rightMoIs[2]
  };

  // Diffs with theoretical vals:
  Len dCS [3] { Abs(comS[0] - R), Abs(comS[1]), Abs(comS[2]) };

  constexpr MoI JxS  = 0.4 * propMassS * R2;
  constexpr MoI JyzS = JxS + propMassS * R2;  // Steiner's Parallel Axis Theorem
  MoI dJS [3]
    { Abs(moisS[0] - JxS), Abs(moisS[1] - JyzS), Abs(moisS[2] - JyzS) };

  cout << "FILLED SPHERE (w/o Shell):"   << endl
       << "dxC="   << dCS[0] << "\ndyC=" << dCS[1] << "\ndzC=" << dCS[2]
       << endl
       <<   "dJx=" << dJS[0] << "\ndJy=" << dJS[1] << "\ndJz=" << dJS[2]
       << endl     << endl;

  // Finally, the Empty Sphere:
  auto const& leftECoM   = leftHS .GetCoM ();
  auto const& rightECoM  = rightHS.GetCoM ();
  auto const& leftEMoIs  = leftHS .GetMoIs();
  auto const& rightEMoIs = rightHS.GetMoIs();

  Len  comES [3]
  {
    0.5*(leftECoM[0] + rightECoM[0]),
    0.5*(leftECoM[1] + rightECoM[1]),
    0.5*(leftECoM[2] + rightECoM[2])
  };
  MoI  moisES[3]
  {
    leftEMoIs[0] + rightEMoIs[0],
    leftEMoIs[1] + rightEMoIs[1],
    leftEMoIs[2] + rightEMoIs[2]
  };

  // Diffs with theoretical vals:
  Len dCES[3] { Abs(comES[0] - R), Abs(comES[1]), Abs(comES[2]) };

  constexpr MoI JxES  = 2.0/3.0 * emptyMassS * R2;
  constexpr MoI JyzES = JxES    + emptyMassS * R2;  // Steiner's Theorem...
  MoI dJES[3]
    { Abs(moisES[0] - JxES), Abs(moisES[1] - JyzES), Abs(moisES[2] - JyzES) };

  cout << "EMPTY SPHERE (SHELL)"          << endl
       << "dxC="   << dCES[0] << "\ndyC=" << dCES[1] << "\ndzC=" << dCES[2]
       << endl
       <<   "dJx=" << dJES[0] << "\ndJy=" << dJES[1] << "\ndJz=" << dJES[2]
       << endl     << endl;
  return 0;
}

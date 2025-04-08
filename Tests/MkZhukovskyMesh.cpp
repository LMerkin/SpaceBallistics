// vim:ts=2:et
//===========================================================================//
//                          "Tests/MkZhukovskyMesh.cpp":                     //
//              Constructing a 2D Conformal Mesh  for the SU2 Solver         //
//              (Model: Compressible, Inviscid or Viscous Turbulent)         //
//                        around a 2D Thin Elliptical Wing                   //
//===========================================================================//
#include "SpaceBallistics/Types.hpp"
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <bit>
#include <vector>
#include <utility>
#include <unistd.h>

using namespace SpaceBallistics;
using namespace std;

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main(int argc, char* argv[])
{
  //-------------------------------------------------------------------------//
  // Params:                                                                 //
  //-------------------------------------------------------------------------//
  // Ratio of the Ellipse Axes: K = a/b:
  double K       = 10.0;

  // Number of Points along each Ellipse (assumed to constant for all Layers):
  int    N       = 400;

  // The over-all size of the mesh:
  Len  aMax      = 20.0_m;

  // The output file:
  string outFile = "mesh.su2";

  // The FreeStream velocity. It is only required if the mesh is to be generated
  // for a viscous (BoundaryLayer) problem:
  double Moo     = NaN<double>;

  // The default value of y^+ for the thickness (NOT half-thickness!) of the
  // inner-most Mesh layer. Not used unless "Moo" is configured:
  double yPlus   = 10.0;

  // Radial Expansion Ratio: Again, not used unless "Moo" is configured:
  double RER     = 1.1;

  // Possibly modify the above params from the command-line:
  while (true)
  {
    int c = getopt(argc, argv, "K:N:A:M:y:R:o:h");
    if (c < 0)
      break;

    switch (c)
    {
    case 'K':
      K       = atof(optarg);
      break;
    case 'N':
      N       = atoi(optarg);
      break;
    case 'A':
      aMax    = Len(atof(optarg));
      break;
    case 'o':
      outFile = string(optarg);
      break;
    case 'M':
      Moo     = atof(optarg);
      break;
    case 'y':
      yPlus   = atof(optarg);
      break;
    case 'R':
      RER     = atof(optarg);
      break;
    case 'h':
    default:
      cerr << "PARAMS: "                      << endl;
      cerr << "\t-K {A/B EllipticWing Ratio; default: " << K << '}'   << endl;
      cerr << "\t-N {Number of Tangential Intervals; default: "       << N
           << '}' << endl;
      cerr << "\t-A {Semi-Major Axis of the Whole Mesh, m; default: "
           << aMax.Magnitude()  << '}'        << endl;
      cerr << "\t-o {Output File; default: "  << outFile << '}'       << endl;
      cerr << "\t-M {FreeStream Mach Number; default: UnDefined => "
              " No BoundaryLayer}"            << endl;
      cerr << "\t-y {The y^+ value for the Inner-Most Mesh Layer; default: "
           << yPlus << '}' << endl;
      cerr << "\t-R {Radial Expansion Ratio; default: " << RER << '}' << endl;
      cerr << "\t-h: display this help"       << endl;
      return 1;
    }
  }
  // Verify the params:
  if (K    <= 1.0 || N   < 50  || aMax < 1.5_m || outFile.empty() ||
      yPlus < 1.0 || RER < 1.0 || RER  > 2.0   ||
      (IsFinite(Moo) && !IsPos(Moo)))
  {
    cerr << "ERROR: Invalid Param(s)" << endl;
    return 1;
  }
  if (!outFile.ends_with(".su2"))
    outFile.append(".su2");

  // If "Moo" is set, the BoundaryLayer-aware mesh will be constructed:
  bool const withBoundLayer = IsFinite(Moo);

  // The Radius of the circle which is the Zhukovsky Conformal Inverse-Image of
  // the Thin Elliptical AirFoil secton (a Circle).    XXX: This param is NOT a
  // "Len":
  double const r0 = SqRt((K-1.0)/(K+1.0));
	assert(r0 < 1.0);

  // The Semi-Major and Semi-Minor Axes of the Thin Elliptic AirFoil:
  Len const a0(0.5 * (1.0 / r0 + r0));
  Len const b0(0.5 * (1.0 / r0 - r0));
  assert(IsPos(b0) && b0 < a0 && a0 > 1.0_m);

  // Kinematic Viscosity of the Air (at 15 C = 288.15 K):
  // TODO: make it variable:
  constexpr auto nu   = 1.47e-5 * Area(1.0) / 1.0_sec;

  // Speed of Sound again (at 15 C = 288.15 K); TODO: make it variable:
  constexpr Vel  Va   = Vel(340.294);

  // The FreeStream Air Velocity:
  Vel const Uoo       = Moo * Va;

  // XXX: Assuming the Chord of 1 m. Since it appears with the degree (1/14),
  // its actual value has VERY LITTLE EFFECT:
  constexpr Len Chord = 1.0_m;

  // Physical thickness of the Turbulent Boundary Layer:
  Len    hB  =
    withBoundLayer
    ? 8.7706 * yPlus * ((nu / Uoo).IPow<13>() * Chord).RPow<1,14>()
    : Len(NaN<double>);

  //-------------------------------------------------------------------------//
  // Open and initialise the output file:                                    //
  //-------------------------------------------------------------------------//
  FILE* f = fopen(outFile.data(), "w");
  if (f == nullptr)
  {
    cerr << "ERROR: Cannot open \"" << outFile << "\" for writing" << endl;
    return 1;
  }
  // XXX: Currently,  up to 99999999 Points are possible; the actual number will
  // be updated after they are generated:
  constexpr int MaxPoints = 99999999;
  fputs  ("NDIME= 2\nNPOIN= XXXXXXXX\n", f);
  long NPointsOff = ftell(f) - 9;

  //-------------------------------------------------------------------------//
  // Points Generation Loop:                                                 //
  //-------------------------------------------------------------------------//
  Len a             = a0;
  Len b             = b0;
  int NPoints       = 0;
  double const dPhi = TwoPi<double> / double(N);

  // The semi-major axis increment ("da"):
  // Then the maximum "height" (normal side size) of the mesh cell  occurs @
  // phi=Pi/2, and is equal to (a/b) * da; for the inner-most layer, it must
  // be <= hB  if the BondaryLayer is used.
  // Otherwise, we select "da" to make the inner-most layer cells square; for
  // that, consider the point phi=0:
  Len da =
    withBoundLayer
    ? double(b0 / a0) * hB
    : b0 * dPhi;

  // The number of Ellipses constructed so far:
  int NEllipses = 0;
 
  while (LIKELY(a <= aMax))
	{
    // Generate the ellipse with the given "a" and "b":
    for (int i = 0; i < N; ++i)
    {
      double phi    = dPhi * double(i);
      double cosPhi = Cos(phi);
      double sinPhi = Sin(phi);

      Len x = a * cosPhi;
      Len y = b * sinPhi;

      fprintf(f, "%.12lf %.12lf   %d\n", x.Magnitude(), y.Magnitude(), NPoints);
      ++NPoints;
    }
    ++NEllipses;

    // Move to the next Ellipse.
    // If there is no Boundary Layer at all, change nothing. Otherwise, increase
    // "ar", but be careful not to do so prematurely:
    if (withBoundLayer && NEllipses >= 2)
      da *= RER;

    // Adjust "a", then re-calculate "b":
    a       += da;
    double r = (a - SqRt(Sqr(a) - Area(1.0))).Magnitude();
    assert(r > 0.0);
    b      = Len(0.5 * (1.0 / r - r));
    assert(IsPos(b) && b < a);
	}
  // All Points and Ellipses Done!

  if (UNLIKELY(NEllipses < 2))
  {
    cerr << "ERROR: Too few Ellipses generated (" << NEllipses << ')' << endl;
    fclose(f);
    return 1;
  }
  if (UNLIKELY(NPoints > MaxPoints))
  {
    cerr << "ERROR: Too many Points generated (" << NPoints << " > "
         << MaxPoints << ')' << endl;
    fclose(f);
    return 1;
  }

  // Store the number of Points in the file (up to 8 digits):
  fseek  (f, NPointsOff, SEEK_SET);
  fprintf(f, "%8d", NPoints);
  fseek  (f, 0, SEEK_END);

  //-------------------------------------------------------------------------//
  // Generate the Rectangular Mesh Elements:                                 //
  //-------------------------------------------------------------------------//
  int  NElems = (NEllipses - 1) * N;
  fprintf(f, "NELEM= %d\n", NElems);

  for (int e = 0; e <= NEllipses-2; ++e)
  {
    // NB: point indices used in Elems generation are 0-based:
    int pFrom = e * N;
    int pTo   = pFrom + N - 1;

    // Generate Rectangles (Code=9):
    for (int i = pFrom; i <= pTo; ++i)
    {
      int i1 = (i <= pTo-1) ? (i+1) : pFrom;

      fprintf(f, "9 %d %d %d %d   %d\n",
              i, i1, i1+N, i+N, NElems++);
    }
  }
  //-------------------------------------------------------------------------//
  // Generate the Boundaries ("Markers"):                                    //
  //-------------------------------------------------------------------------//
  // The Inner-Most Ellipse ("Wall"):
  int wFrom = 0;
  int wTo   = N - 1;
  fprintf(f, "NMARK= 2\nMARKER_TAG= Wall\nMARKER_ELEMS= %d\n",
          wTo - wFrom + 1);

  for (int i = wFrom; i <= wTo; ++i)
    fprintf(f, "3 %d %d\n", i, (i != wTo) ? (i+1) : wFrom);

  // The Outer-Most Ellipse ("FarAway"):
  int fFrom = (NEllipses - 1) * N;
  int fTo   = fFrom  + N - 1;
  assert(fTo == NPoints  - 1);

  fprintf(f, "MARKER_TAG= FarAway\nMARKER_ELEMS= %d\n",
          fTo - fFrom + 1);

  for (int i = fFrom; i <= fTo; ++i)
    fprintf(f, "3 %d %d\n", i, (i != fTo) ? (i+1) : fFrom);

  //-------------------------------------------------------------------------//
  // Store the generation params at the end:                                 //
  //-------------------------------------------------------------------------//
  if (withBoundLayer)
    fprintf(f, "%% M = %lf, y+ = %lf, RER = %lf\n", Moo, yPlus, RER);

  fclose(f);
  return 0;
}

// vim:ts=2:et
//===========================================================================//
//                        "Tests/MkZhukovskyMesh.cpp":                       //
//            Constructing a 2D Conformal Mesh for the SU2 Solver            //
//                       around a Thin Elliptical Plate                      //
//===========================================================================//
#include "SpaceBallistics/Types.hpp"
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <unistd.h>

int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  //-------------------------------------------------------------------------//
  // Params:                                                                 //
  //-------------------------------------------------------------------------//
  // Ratio of the Ellipse Axes: K = a/b:
  double K       = 10.0;

  // Number of Points along each Ellipse:
  int    NP      = 100;

  // The over-all size of the mesh:
  Len  aMax      = 7.0_m;

  // The output file:
  string outFile = "mesh.su2";

  // Possibly modify them from the command-line params:
  while (true)
  {
    int c = getopt(argc, argv, "K:N:A:o:");
    if (c < 0)
      break;

    switch (c)
    {
    case 'K':
      K       = atof(optarg);
      break;
    case 'N':
      // NB: "NP" must be even:
      NP      = 2 * ((atoi(optarg) + 1) / 2);
      break;
    case 'A':
      aMax    = Len(atof(optarg));
      break;
    case 'o':
      outFile = string(optarg);
      break;
    default:
      assert(false);
    }
  }
  // Verify the params:
  if (K <= 1.0 || NP < 10 || NP % 2 != 0 || aMax < 1.5_m || outFile.empty())
  {
    cerr << "ERROR: Invalid Param(s)" << endl;
    return 1;
  }
  if (!outFile.ends_with(".su2"))
    outFile.append(".su2");

  // The Radius of the circle which is the Zhukovsky Conformal Inverse-Image of
  // the Thin Elliptical AirFoil secton (a Circle).    XXX: This param is NOT a
  // "Len":
  double const r0 = SqRt((K-1.0)/(K+1.0));
	assert(r0 < 1.0);

  // The Semi-Major and Semi-Minor Axes of the Thin Elliptic AirFoil:
  Len const a0(0.5 * (1.0 / r0 + r0));
  Len const b0(0.5 * (1.0 / r0 - r0));
  assert(IsPos(b0) && b0 < a0 && a0 > 1.0_m);

  // The transversal mesh size (in m) is approximately the same as the tangent-
  // ial mesh size, to get a near-square mesh:
  Len const h0S = Pi<double> * (a0 + b0) / double(NP);

  //-------------------------------------------------------------------------//
  // Open and initialise the output file:                                    //
  //-------------------------------------------------------------------------//
  FILE* f = fopen(outFile.data(), "w");
  if (f == nullptr)
  {
    cerr << "ERROR: Cannot open \"" << outFile << "\" for writing" << endl;
    return 1;
  }
  // XXX: Currently,  up to 99999 points are possible; the actual number will
  // be updated after they are generated:
  constexpr int MaxPoints = 99999;
  fputs("NDIME= 2\nNPOIN= XXXXX\n", f);

  // Generate the Sin and Cos values of the Polar Angles to be used:
  double cosPhi[NP];
  double sinPhi[NP];

  for (int i = 0; i < NP; ++i)
  {
    double phi = TwoPi<double> * double(i) / double(NP);
    cosPhi[i]  = Cos(phi);
    sinPhi[i]  = Sin(phi);
  }
  //-------------------------------------------------------------------------//
  // Points Generation Loop:                                                 //
  //-------------------------------------------------------------------------//
  Len    a = a0;
  Len    b = b0;
  double r = r0;
  int    TotPts    = 0;
  int    NEllipses = 0;

  while (LIKELY(a <= aMax))
	{
    // Generate the ellipse with the given "a" and "b":
    for (int i = 0; i < NP; ++i)
    {
      Len x = a * cosPhi[i];
      Len y = b * sinPhi[i];
      fprintf(f, "%.12lf %.12lf   %d\n", x.Magnitude(), y.Magnitude(), TotPts);
      ++TotPts;
    }
    ++NEllipses;

    // Move to the next elliptic layer; "h0S" will be the step wrt "a":
    a += h0S;
    r  = (a - SqRt(Sqr(a) - Area(1.0))).Magnitude();
    b  = Len(0.5 * (1.0 / r - r));
	}
  // All Points Done!

  if (UNLIKELY(NEllipses < 2))
  {
    cerr << "ERROR: Too few Ellipses generated (" << NEllipses << ')' << endl;
    fclose(f);
    return 1;
  }
  if (UNLIKELY(TotPts > MaxPoints))
  {
    cerr << "ERROR: Too many Points generated (" << TotPts << ')' << endl;
    fclose(f);
    return 1;
  }
  assert(TotPts == NEllipses * NP);

  // Store the number of points in the file:
  fseek  (f, 16, SEEK_SET);
  fprintf(f, "%5d", TotPts);
  fseek  (f,  0, SEEK_END);

  //-------------------------------------------------------------------------//
  // Generate the Squre Mesh Elements:                                       //
  //-------------------------------------------------------------------------//
  // Consider "layers of near-squares", each layer located between two ellipses.
  // Thus, the number of layers in (NEllipses-1). Each layer contains NP squar-
  // es, so:
  int NSquares = (NEllipses - 1) * NP;
  fprintf(f, "NELEM= %d\n", NSquares);

  int s = 0;
  for (int e = 0; e < NEllipses-1; ++e)  // For all Ellipses but the last
  {
    // First, process all points on this Ellipse in [0 .. NP-2]:
    for (int i = 0; i < NP-1; ++i, ++s)
      // NB: Code=9 is for the Square;
      // the "generic" square is made of points (s, s+1, s+NP+1, s+NP), where
      // the first two lie on this ellipse, and the last two  -- on the  next
      // ellipse:
      fprintf(f, "9 %d %d %d %d   %d\n", s, s+1,    s+1+NP, s+NP, s);

    // The last square (i==NP-1):
    fprintf  (f, "9 %d %d %d %d   %d\n", s, s-NP+1, s+1,    s+NP, s);
    ++s;
  }
  assert(s == TotPts - NP && s == NSquares);   // Except the outer-most ellipse

  //-------------------------------------------------------------------------//
  // Generate the Boundaries ("Markers"):                                    //
  //-------------------------------------------------------------------------//
  // "Wall":
  fprintf(f, "NMARK= 2\nMARKER_TAG= Wall\nMARKER_ELEMS= %d\n", NP);
  for (int p = 0; p < NP; ++p)
    fprintf(f, "3 %d %d\n", p, (p != NP-1)     ? (p+1) : 0);

  fprintf(f, "MARKER_TAG= FarAway\nMARKER_ELEMS= %d\n", NP);
  for (int p = TotPts-NP; p < TotPts; ++p)
    fprintf(f, "3 %d %d\n", p, (p != TotPts-1) ? (p+1) : (TotPts-NP));

  // All Done:
  fclose(f);
  return 0;
}

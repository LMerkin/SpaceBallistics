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

  // Number of Points along each Ellipse
  // (in the Outer, ie Non-Boundary, Region); must be a power of 2:
  int    NO      = 512;

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
  double RER     = 1.2;

  // Number of inner-most Mesh layers of constant thickness (before it starts
  // increasing); reasonable values are 1..3:
  int    nInner  = 1;

  // Possibly modify the above params from the command-line:
  while (true)
  {
    int c = getopt(argc, argv, "K:N:A:M:y:R:n:o:h");
    if (c < 0)
      break;

    switch (c)
    {
    case 'K':
      K       = atof(optarg);
      break;
    case 'N':
      // NB: "NO" must be a power of 2:
      NO      = int(bit_ceil(unsigned(atoi(optarg))));
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
    case 'n':
      nInner  = atoi(optarg);
      break;
    case 'h':
    default:
      cerr << "PARAMS: "                      << endl;
      cerr << "\t-K {A/B EllipticWing Ratio; default: " << K << '}'   << endl;
      cerr << "\t-N {Number of Elliptic Arcs in the Far-Away Region; "
              "default: " << NO << '}'        << endl;
      cerr << "\t-A {Semi-Major Axis of the Whole Mesh, m; default: "
           << aMax.Magnitude()  << '}'        << endl;
      cerr << "\t-o {Output File; default: "  << outFile << '}'       << endl;
      cerr << "\t-M {FreeStream Mach Number; default: UnDefined => "
              " No BoundaryLayer}"            << endl;
      cerr << "\t-y {The y^+ value for the Inner-Most Mesh Layer; default: "
           << yPlus << '}' << endl;
      cerr << "\t-R {Radial Expansion Ratio; default: " << RER << '}' << endl;
      cerr << "\t-n {Number of Inner-Most Mesh Layers of Const Thickness; "
              "default: "  << nInner << '}'   << endl;
      cerr << "\t-h: display this help"       << endl;
      return 1;
    }
  }
  // Verify the params:
  if (K    <= 1.0 || NO  < 16  || aMax < 1.5_m || outFile.empty() ||
      yPlus < 1.0 || RER < 1.0 || RER  > 2.0   || nInner < 1      ||
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
  // its actual value has very little effect:
  constexpr Len Chord = 1.0_m;

  // Then the params of the inner-most mesh layer (corrsep to the Viscous
  // BoundaryLayer) can be approximated as follows:
  int    NB  = NO;
  Len    hB(NaN<double>);   // Not used unless "withBoundLayer" is set
  double AR0 = 1.0;

  if (withBoundLayer)
  {
    // Physical thickness of the Turbulent Boundary Layer:
    hB = 8.7706 * yPlus * ((nu / Uoo).IPow<13>() * Chord).RPow<1,14>();

    // We begin with the assumption that the inner-most layer elements are
    // roughly Squares. Then the max size of the square side is around phi=Pi/2,
    // and is equal to (a0 * dPhi).
    // We require that HALF of this size (ie the distance from the Wall to the
    // Square center) should not exceed "hB":
    // a0 * dPhi / 2 <= hB,     and thus dPhi <= 2 * hB / a0;
    // if we have  "NB" intervals, then  dPhi  = 2 * Pi / NB, and thus
    // NB >= Pi * a0 / hB;
    // furthermore, we round NB up (in most cases) to a Power-of-2:
    //
    NB = int(bit_ceil(unsigned(Round(Pi<double> * double(a0 / hB)))));
    if (NB < NO)
    {
      cerr << "ERROR: Inconsistency: Got NB=" << NB << ", NO=" << NO << endl;
      return 1;
    }
    // The Aspect Ratio (Height / Width) of the Square @ phi=Pi/2:
    AR0 = double(2.0 * hB / (TwoPi<double> / double(NB) * a0));
  }
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
  Len a        = a0;
  Len b        = b0;
  int NPoints  = 0;

  // The number of Points to be generated for the curr Ellipse: variable for the
  // ("fine", "inner") mesh, or constant  for the ("outer", "coarse") mesh:
  int const NP = withBoundLayer ? NB : NO;
  int np       = NP;

  // For each Ellipse constructed, we memoise its range of Points:
  vector<pair<int,int>> ellipses;

  // The number of ellipses generated with the current "np":
  int ne       = 0;
 
  // The Aspect Ratio of the rectangles to be generated (Height / Width):
  // variable:
  double ar    = 1.0;

  // The "ar" limit after which we double the mesh size (must be > 2):
  constexpr double MaxAR = 2.75;

  while (LIKELY(a <= aMax))
	{
    // Generate the ellipse with the given "a", "b" and "np":
    double dPhi = TwoPi<double> / double(np);
    for (int i = 0; i < np; ++i)
    {
      double phi    = dPhi * double(i);
      double cosPhi = Cos(phi);
      double sinPhi = Sin(phi);

      Len x = a * cosPhi;
      Len y = b * sinPhi;

      fprintf(f, "%.12lf %.12lf   %d\n", x.Magnitude(), y.Magnitude(), NPoints);
      ++NPoints;
    }
    // Save the range of indices of points making this ellipse:
    ellipses.push_back   (make_pair(NPoints-np, NPoints-1));
    ++ne;

    // The semi-major axis increment  is the orthogonal step at phi=0,  which
    // is computed using the tangential step (b * dPhi), incl the "AR0" factor;
    // it may be overridden in the following:
    Len da = b * dPhi * AR0;

    // Move to the next Ellipse.
    // If we are still in the "nInner" Boundary Mesh Layers, or there is no
    // Boundary Layer at all, change nothing.
    // Otherwise, increase "ar", but be careful not to do so prematurely:
    //
    if (withBoundLayer && int(ellipses.size()) >= nInner+1 && ne >= 2 &&
        np > NO)
    {
      ar *= RER;

      // If we got a mesh which is "too much stretched" in the radial direction,
      // it's time to double the cell size (incl the base), until we reach the
      // "far away" region:
      if (ar > MaxAR)
      {
        // "da" will be for the "interleaving" layer:
        da *= (ar - 2.0);

        // New "np" and "ar" will be for the next ellipse containing ~2 times
        // less arcs.  At the beginning of the next iteration, "dPhi" will be
        // re-calculated:
        np  = max(np / 2, NO);
        ar  = 1.0;
        ne  = 0;
      }
      else
        // Use the new "ar" to stretch "da" from its generic value:
        da *= ar;
    }
    // In all cases:
    assert(NP % np == 0);

    // For the inner-most layers, make sure that "da" (NOT da/2 !) does not
    // exceed "hB":
    if (UNLIKELY(withBoundLayer && int(ellipses.size()) <= nInner))
      da = min(da,  hB);

    // Adjust "a", then re-calculate "b":
    a     += da;
    Len r  = a - SqRt(Sqr(a) - Area(1.0));
    assert(IsPos(r));
    b      = 0.5 * (Area(1.0) / r - r);
    assert(IsPos(b) && b < a);
	}
  // All Points and Ellipses Done!

  int NEllipses = int(ellipses.size());
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
  // Generate the Squre and Triangular Mesh Elements:                        //
  //-------------------------------------------------------------------------//
  // Consider "layers of near-squares", each layer located between two ellipses.
  // The total number of squares is determined dynamically:
  int MaxElems   = 99999999;
  fputs   ("NELEM= XXXXXXXX\n", f);
  long NElemsOff = ftell(f)  -  9;
  int  NElems    = 0;

  for (int e = 0; e <= NEllipses-2; ++e)
  {
    // First, process all Points  on this Ellipse:
    auto [pFrom, pTo] = ellipses[size_t(e)];

    // The number of Points on this Ellipse:
    int  ne           = pTo - pFrom + 1;
    assert(ne >= NO && ne % 2 == 0);

    // The next Ellipse may have either the same number of Points, or 2 times
    // less:
    auto [nFrom, nTo] = ellipses[size_t(e+1)];
    int  nn           = nTo - nFrom + 1;
    assert(nn >= 2 && nn % 2 == 0);

    bool   same = (ne == nn);
    assert(same   ||  ne == nn * 2);
    // In particular, if e <= nInner, we must have the "same" flag set, other-
    // wise we would not be able to complete the inner squares:
    assert(e > nInner || same);

    // If the "same" flag is set, then we use all Points of the Ellipse "e" as
    // the Square vertices;  otherwise, use each 2nd point as Triangle vertices:
    int step = same ? 1 : 2;

    // "Ort": the function which, for a given Point idx ("i"), returns the idx
    // of the corresp point on the next Ellipse, ie in the orthogonal (wrt the
    // curr Ellipse) direction from "i":
    auto Ort =
      [pFrom,  nFrom, step](int i) -> int
      {
        assert(i >= pFrom && (i - pFrom) % step == 0);
        return      nFrom +  (i - pFrom) / step;
      };
    assert(pTo+1 - ne == pFrom);

    if (same)
    {
      // Generate Squares (Code=9), except the last one (i==pTo):
      for (int i = pFrom; i <= pTo; ++i)
      {
        int i1 = (i <= pTo-1) ? (i+1) : pFrom;

        fprintf(f, "9 %d %d %d %d   %d\n",
                i, i1, Ort(i1), Ort(i), NElems++);
      }
    }
    else
    {
      // Generate Triangles (Code=5):
      assert(pFrom % 2 == 0 && (pTo-1) % 2 == 0);

      for (int i = pFrom; i <= pTo-1; i += 2)
      {
        // The right-most point (i+2): again, (pTo+1) is wrapped to "pFrom":
        int i2  = (i <= pTo-2) ? (i+2) : pFrom;
        int Oi  = Ort(i);
        int Oi2 = Ort(i2);

        // The left  triangle ("pointing outwards"):
        fprintf(f, "5 %d %d %d   %d\n", i,   i+1, Oi,  NElems++);

        // The mid   triangle ("pointing inwards"):
        fprintf(f, "5 %d %d %d   %d\n", Oi,  i+1, Oi2, NElems++);

        // The right triangle ("pointing outwards"):
        fprintf(f, "5 %d %d %d   %d\n", i+1, i2,  Oi2, NElems++);
      }
    }
  }
  if (UNLIKELY(NElems > MaxElems))
  {
    cerr << "ERROR: Too many Elems generated (" << NElems << " > "
         << MaxElems << ')' << endl;
    fclose(f);
    return 1;
  }
  // Store the number of Elems in the file (up to 8 digits):
  fseek  (f, NElemsOff, SEEK_SET);
  fprintf(f, "%8d",     NElems);
  fseek  (f, 0,         SEEK_END);

  //-------------------------------------------------------------------------//
  // Generate the Boundaries ("Markers"):                                    //
  //-------------------------------------------------------------------------//
  // The Inner-Most Ellipse ("Wall"):
  auto  [wFrom, wTo] = ellipses.front();
  assert(wFrom == 0 && wTo == NP - 1);

  fprintf(f, "NMARK= 2\nMARKER_TAG= Wall\nMARKER_ELEMS= %d\n",
          wTo - wFrom + 1);

  for (int i = wFrom; i <= wTo; ++i)
    fprintf(f, "3 %d %d\n", i, (i != wTo) ? (i+1) : wFrom);

  // The Outer-Most Ellipse ("FarAway"):
  auto  [fFrom, fTo] = ellipses.back();

  fprintf(f, "MARKER_TAG= FarAway\nMARKER_ELEMS= %d\n",
          fTo - fFrom + 1);

  for (int i = fFrom; i <= fTo; ++i)
    fprintf(f, "3 %d %d\n", i, (i != fTo) ? (i+1) : fFrom);

  //-------------------------------------------------------------------------//
  // Store the generation params at the end:                                 //
  //-------------------------------------------------------------------------//
  if (withBoundLayer)
    fprintf(f, "%% M = %lf, y+ = %lf, nInner = %d, RER = %lf\n",
            Moo, yPlus, nInner, RER);

  fclose(f);
  return 0;
}

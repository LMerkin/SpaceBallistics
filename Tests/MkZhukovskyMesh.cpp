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

namespace
{
  //=========================================================================//
  // Approximation for the Perimeter of the Ellipse:                         //
  //=========================================================================//
  Len PerEll(Len a, Len b)
  {
    assert(a > b && IsPos(b));
    double h = Sqr(double((a-b)/(a+b)));
    assert(0 < h && h < 1);

    constexpr double c1 =    1.0 /       4.0;
    constexpr double c2 =    1.0 /      64.0;
    constexpr double c3 =    1.0 /     256.0;
    constexpr double c4 =   25.0 /   16384.0;
    constexpr double c5 =   49.0 /   65536.0;
    constexpr double c6 =  441.0 / 1048576.0;
    constexpr double c7 = 1089.0 / 4194304.0;
    double corr =
      ((((((c7 * h + c6) * h + c5) * h + c4) * h + c3) * h + c2) * h + c1) * h;

    return Pi<double> * (a + b) * (1.0 + corr);
  }
}

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
  int    NO      = 256;

  // The over-all size of the mesh:
  Len  aMax      = 7.5_m;

  // The output file:
  string outFile = "mesh.su2";

  // The FreeStream velocity. It is only required if the mesh is to be generated
  // for a viscous (Boundary Layer) problem:
  Vel    Uoo(NaN<double>);

  // Whether the Mesh is to be generated fotr use with Wall Functions; only eff-
  // ective if "Uoo" is set:
  bool   withWallFuncs  = false;

  // The default value of y^+ to be used with Wall Functions:
  double yPlusWall      = 50.0;

  // Possibly modify them from the command-line params:
  while (true)
  {
    int c = getopt(argc, argv, "K:N:A:U:Y:o:Wh");
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
    case 'U':
      Uoo     = Vel(atof(optarg));
      break;
    case 'W':
      withWallFuncs = true;
      break;
    case 'Y':
      yPlusWall     = atof(optarg);
      break;
    case 'h':
    default:
      cerr << "PARAMS: "                      << endl;
      cerr << "\t-K {A/B EllipticWing Ratio; default: " << K << '}'   << endl;
      cerr << "\t-N {Number of Elliptic Arcs in the Non-Boundary Region; "
              "default: " << NO << '}'        << endl;
      cerr << "\t-A {Semi-Major Axis of the Whole Mesh, m; default: "
           << aMax.Magnitude()  << '}'        << endl;
      cerr << "\t-o {Output File; default: "  << outFile << '}'       << endl;
      cerr << "\t-U {FreeStream Velocity, m/sec; default: UnDefined}" << endl;
      cerr << "\t-W: set withWallFuncs=true"  << endl;
      cerr << "\t-Y {y^+ for use with WallFuncs; default: " << yPlusWall
           << '}' << endl;
      cerr << "\t-h: display this help"       << endl;
      return 1;
    }
  }
  // Verify the params:
  if (K <= 1.0 || NO < 16 || aMax < 1.5_m || outFile.empty() ||
      yPlusWall < 1.0     || (IsFinite(Uoo) && !IsPos(Uoo)))
  {
    cerr << "ERROR: Invalid Param(s)" << endl;
    return 1;
  }
  if (!outFile.ends_with(".su2"))
    outFile.append(".su2");

  // If "Uoo" is set, the BoundaryLayer-aware mesh will be constructed:
  bool const withBoundLayer = IsFinite(Uoo);

  // The Radius of the circle which is the Zhukovsky Conformal Inverse-Image of
  // the Thin Elliptical AirFoil secton (a Circle).    XXX: This param is NOT a
  // "Len":
  double const r0 = SqRt((K-1.0)/(K+1.0));
	assert(r0 < 1.0);

  // The Semi-Major and Semi-Minor Axes of the Thin Elliptic AirFoil:
  Len const a0(0.5 * (1.0 / r0 + r0));
  Len const b0(0.5 * (1.0 / r0 - r0));
  assert(IsPos(b0) && b0 < a0 && a0 > 1.0_m);

  // The minimal transversal mesh size (in m) in the BoundaryLayer depends on
  // the value of "y+":
  double const yPlus = withWallFuncs ? yPlusWall : 1.0;

  // Kinematic Viscosity of the Air (at 15 C = 288.15 K):
  constexpr auto nu  = 1.47e-5 * Area(1.0) / 1.0_sec;

  // Then the thickness of the inner-most mesh layer (corrsep to the Viscous
  // Boundary Layer) can be approximated as:
  int NB = NO;
  if (withBoundLayer)
  {
    Len perim0 = PerEll(a0, b0);
    Len hB     = 7.5 * yPlus * ((nu / Uoo).IPow<13>() * 1.0_m).RPow<1,14>();
    NB         = int(bit_ceil(unsigned(Round(double(perim0 / hB)))));
    if (NB <= NO)
    {
      cerr << "ERROR: Inconsistency: Got NB=" << NB << ", NO=" << NO << endl;
      return 1;
    }
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
  // XXX: Currently,  up to 9999999 Points are possible; the actual number will
  // be updated after they are generated:
  constexpr int MaxPoints = 9999999;
  fputs  ("NDIME= 2\nNPOIN= XXXXXXX\n", f);
  long NPointsOff = ftell(f) - 8;

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

  while (LIKELY(a <= aMax))
	{
    // Perimeter of this ellipse and the "small arc" length:
    Len dl = PerEll(a, b) / double(np);

    // Generate the ellipse with the given "a" and "b":
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
    ellipses.push_back   (make_pair(NPoints-np, NPoints-1));

    // Move to the next Ellipse:
    // "fine" mesh is enlarged by the factor of 2 at each increase (UNLESS mov-
    // ing from the 0th to the 1st Ellipse); "outer" mesh is unchanged:
    np =
      withBoundLayer
      ? ((ellipses.size() == 1) ? np : max(np / 2, NO))
      : NO;
    assert(NP % np == 0);

    // The semi-major axis increment should be equal to the orthogonal step,
    // which is in turn approx equal to the "dl" (if we want the mesh elements
    // to be ~ squares):
    a        += dl;
    double r  = (a - SqRt(Sqr(a) - Area(1.0))).Magnitude();
    b         = Len(0.5 * (1.0 / r - r));
    assert(b < a);
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

  // Store the number of Points in the file (up to 6 digits):
  fseek  (f, NPointsOff, SEEK_SET);
  fprintf(f, "%7d", NPoints);
  fseek  (f, 0, SEEK_END);

  //-------------------------------------------------------------------------//
  // Generate the Squre Mesh Elements:                                       //
  //-------------------------------------------------------------------------//
  // Consider "layers of near-squares", each layer located between two ellipses.
  // The total number of squares is determined dynamically:
  int MaxSquares = 9999999;
  fputs   ("NELEM= XXXXXXX\n", f);
  long NSqrsOff  = ftell(f) - 8;
  int  NSquares  = 0;

  for (int e = 0; e <= NEllipses-2; ++e)
  {
    // First, process all Points  on this Ellipse except the last one:
    auto [pFrom, pTo] = ellipses[size_t(e)];
    // The number of Points on this Ellipse:
    int  ne           = pTo - pFrom + 1;
    assert(ne >= 2 && ne % 2 == 0);

    // The next Ellipse may have either the same number of Points, or 2 times
    // less:
    auto [nFrom, nTo] = ellipses[size_t(e+1)];
    int  nn           = nTo - nFrom + 1;
    assert(nn >= 2 && nn % 2 == 0);

    bool   same = (ne == nn);
    assert(same   ||  ne == nn * 2);
    // In particular, if e==0, we must have the "same" flag set, otherwise we
    // would not be able to complete the squares:
    assert(e != 0 || same);

    // If the "same" flag is set, then we use all Points of the Ellipse "e" as
    // the Square vertices;  otherwise, use each 2nd point:
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

    // "Generic" maximum "i":
    int iGMax  = pTo - (same ? 1 : 3);

    for (int i = pFrom; i <= iGMax; i += step)
      // NB: Code=9 is for the Square:
      fprintf(f, "9 %d %d %d %d   %d\n",
              i, i + step, Ort(i + step), Ort(i), NSquares++);

    // The last square (i == iGMax + step):
    int iMax   = pTo - (same ? 0 : 1);
    assert(iMax == iGMax + step);
    assert(iMax  + step  - ne == pFrom);

    fprintf(f, "9 %d %d %d %d   %d\n",
            iMax, pFrom, Ort(pFrom), Ort(iMax), NSquares++);
  }
  if (UNLIKELY(NSquares > MaxSquares))
  {
    cerr << "ERROR: Too many Squares generated (" << NSquares << " > "
         << MaxSquares << ')' << endl;
    fclose(f);
    return 1;
  }
  // Store the number of Squares in the file (up to 6 digits):
  fseek  (f, NSqrsOff, SEEK_SET);
  fprintf(f, "%7d",    NSquares);
  fseek  (f, 0,        SEEK_END);

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

  // All Done:
  fclose(f);
  return 0;
}

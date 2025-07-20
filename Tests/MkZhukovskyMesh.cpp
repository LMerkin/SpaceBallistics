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

  // Number of Points along each Ellipse (in the Far-Away Region); must be a
  // power of 2:
  int    NF      = 512;

  // The over-all size of the mesh:
  Len    aMax    = 20.0_m;

  // The output file:
  string outFile = "mesh.su2";

  // THE FOLLOWING PARAMS ARE APPLICABLE TO THE BOUNDARY-LAYER CASE ONLY:

  // The FreeStream velocity. It is only required if the mesh is to be generated
  // for a viscous (BoundaryLayer) problem:
  double Moo     = NaN<double>;

  // The default value of y^+ for the thickness (NOT half-thickness!) of the
  // inner-most Mesh layer. Not used unless "Moo" is configured:
  double yPlus   = 1.0;

  // Radial Expansion Ratio: Again, not used unless "Moo" is configured:
  double RER     = 1.1;

  // Number of inner-most Mesh layers of constant thickness (before it starts
  // increasing); reasonable values are 1..3:
  int    nInner  = 1;

  // Aspect Ratio (TangentialLength / NormalThickness) of rectangles in the
  // inner-most layer:
  double ar0     = 1000.0;

  // Whether we generate "fixed-angular" mesh (ie the mesh cell size changes in
  // the radial direction only); by default this is False, so that we apply the
  // "doubling" stratregy in the angular domain; 
  bool fixedAng = false;

  // Verbose mode?
  bool verbose  = false;

  // CONSTANTS (TODO: May need to make them configurable as well):
  // Kinematic Viscosity of the Air (at 15 C = 288.15 K):
  constexpr auto nu   = 1.47e-5 * Area(1.0) / 1.0_sec;

  // Speed of Sound again (at 15 C = 288.15 K):
  constexpr Vel  Va   = Vel(340.294);

  //-------------------------------------------------------------------------//
  // Possibly modify the above params from the command-line:                 //
  //-------------------------------------------------------------------------//
  while (true)
  {
    int c = getopt(argc, argv, "K:N:A:M:y:R:n:s:o:fvh");
    if (c < 0)
      break;

    switch (c)
    {
    case 'K':
      K       = atof(optarg);
      break;
    case 'N':
      // NB: "NF" must be a power of 2:
      NF       = int(bit_ceil(unsigned(atoi(optarg))));
      break;
    case 'A':
      aMax     = Len(atof(optarg));
      break;
    case 'o':
      outFile  = string(optarg);
      break;
    case 'M':
      Moo      = atof(optarg);
      break;
    case 'y':
      yPlus    = atof(optarg);
      break;
    case 'R':
      RER      = atof(optarg);
      break;
    case 'n':
      nInner   = atoi(optarg);
      break;
    case 's':
      ar0      = atof(optarg);
      break;
    case 'f':
      fixedAng = true;
      break;
    case 'v':
      verbose  = true;
      break;
    case 'h':
    default:
      cerr << "PARAMS: "                      << endl;
      cerr << "\t-h: display this help"       << endl;
      cerr << "\t-K {A/B EllipticWing Ratio; default: " << K << '}'   << endl;
      cerr << "\t-N {Number of Elliptic Arcs in the Far-Away Region; "
              "default: " << NF << '}'        << endl;
      cerr << "\t-A {Semi-Major Axis of the Whole Mesh, m; default: "
           << aMax.Magnitude()  << '}'        << endl;
      cerr << "\t-o {Output File; default: "  << outFile << '}'       << endl;
      cerr << "The following params are for the BoundaryLayer case only:"
           << endl;
      cerr << "\t-M {FreeStream Mach Number; default: UnDefined => "
              " No BoundaryLayer}"            << endl;
      cerr << "\t-y {The y^+ value for the Inner-Most Mesh Layer; default: "
           << yPlus << '}' << endl;
      cerr << "\t-R {Radial Expansion Ratio; default: " << RER << '}' << endl;
      cerr << "\t-n {Number of Inner-Most Mesh Layers of Const Thickness; "
              "default: "  << nInner << '}'   << endl;
      cerr << "\t-s {Max Aspect Ratio of Inner-Most Mesh Rectangles; default: "
           << ar0 << '}' << endl;
      cerr << "\t-f : Construct a \"Fixed-Angular\" Mesh; default: "
           << (fixedAng ? "true" : "false")   << endl;
      cerr << "\t-v : Verbose mode; default: "
           << (verbose  ? "true" : "false")   << endl;
      return 1;
    }
  }
  // Verify the params:
  constexpr int MinN = 32;
  if (K     <= 1.0 || NF  < MinN || aMax < 1.5_m || outFile.empty() ||
      yPlus <= 0.0 || RER < 1.0  || RER  > 2.0   || nInner < 1      ||
      ar0   <= 0.0 || (IsFinite(Moo) && !IsPos(Moo)))
  {
    cerr << "ERROR: Invalid Param(s)" << endl;
    return 1;
  }
  if (!outFile.ends_with(".su2"))
    outFile.append(".su2");

  //-------------------------------------------------------------------------//
  // Derived Params:                                                         //
  //-------------------------------------------------------------------------//
  // The Radius of the circle which is the Zhukovsky Conformal Inverse-Image of
  // the Thin Elliptical AirFoil secton (a Circle).    XXX: This param is NOT a
  // "Len":
  double const r0 = SqRt((K-1.0)/(K+1.0));
  assert(r0 < 1.0);

  // The Semi-Major and Semi-Minor Axes of the Thin Elliptic AirFoil:
  Len const a0(0.5 * (1.0 / r0 + r0));
  Len const b0(0.5 * (1.0 / r0 - r0));
  assert(IsPos(b0) && b0 < a0 && a0 > 1.0_m);

  // Then the params of the inner-most mesh layer (corrsep to the Viscous
  // BoundaryLayer) can be approximated as follows:
  int    NB   = NF;

  // "hB" is not used unless "withBoundLayer" is set:
  Len    hB(NaN<double>);

  // "da" is the initial "a" increment; the default value is for the case w/o
  // BoundaryLayer (dPhi * b0); otherwise, it is overriden below:
  double dPhi = TwoPi<double> / double(NF);
  Len    da   = dPhi * b0;

  // If "Moo" is set, the BoundaryLayer-aware mesh will be constructed:
  bool const withBoundLayer = IsFinite(Moo);

  //-------------------------------------------------------------------------//
  // The BoundaryLayer Setup:                                                //
  //-------------------------------------------------------------------------//
  if (withBoundLayer)
  {
    // XXX: Assuming the Chord of 1 m. Since it appears with the degree (1/14),
    // its actual value has very little effect:
    constexpr Len Chord = 1.0_m;

    // The FreeStream Air Velocity:
    Vel const Uoo = Moo * Va;

    // Physical thickness of the Turbulent Boundary Layer (as computed with the
    // given "yPlus"):
    hB = 8.7706 * yPlus * ((nu / Uoo).IPow<13>() * Chord).RPow<1,14>();

    // Then the max thickness of the inner-most mesh layer is @ phi=Pi/2, and
    // is equal to (a0/b0) * da, where  "da" is the thickness @ phi=0.  Thus,
    // we need     (a0/b0) * da  <= hB, where (a0/b0) = K:
    da = hB / K;

    // Then apply the "ar0" to the inner-most cells @ phi=0 (and since the map
    // is conformal, the same ratio will apply automatically to all other "phi"
    // didirections):
    // their tangential size    is (b0 * dPhi),
    // the normal (radial) size is da,
    // and we must have
    // b0 * dPhi <= da * ar0, so the (yet-unadjusted) "dPhi" will be:
    dPhi = double(da / b0) * ar0;

    // The number of points (and arcs) at the "boundary wall": must be a power
    // of 2 as well, so the actual "dPhi" may be smaller than the above value:
    NB   = int(bit_ceil(unsigned(Round(TwoPi<double> / dPhi))));
    if (NB < MinN)
    {
      cerr << "ERROR: Got too low NB=" << NB << endl;
      return 1;
    }
    // Adjust "dPhi" for the sake of consistency:
    dPhi = TwoPi<double> / double(NB);

    // If a "fixed-angular" mesh is requested, then "NB" will be used for the
    // whole mesh, overriding "NF"  (and the "NF"-based "da" has already been
    // overridden):
    if (fixedAng)
      NF = NB;
    else
    if (NB < NF)
    {
      cerr << "ERROR: Inconsistency: Got NB=" << NB << ", NF=" << NF << endl;
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
  // ("fine", "inner") mesh, or constant  for the ("far-away", "coarse") mesh:
  int const NP = withBoundLayer ? NB : NF;
  int       np = NP;

  // For each Ellipse constructed, we memoise its range of Points:
  vector<pair<int,int>> ellipses;

  // The number of ellipses generated with the current "np":
  int       nc = 0;
 
  while (LIKELY(a <= aMax))
  {

    // Generate the ellipse with the given "a", "b" and "np":
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
    ++nc;

    // For the next Ellipse:
    // If we are still in the "nInner" Boundary Mesh Layers,  or there is no
    // Boundary Layer at all, or we are in the "far away" region,   keep the
    // current "da".
    // Otherwise, may increase "da", but be careful not to do so prematurely:
    //
    if (withBoundLayer && int(ellipses.size()) >= nInner+1 && nc >= 2)
    {
      // Tangential Width @ phi=0:
      Len dw = b * dPhi;

      if (np > NF)
      {
        // We are NOT in the "far-away" region yet, so expand "da":
        da *= RER;

        // If we got a mesh which is "too much stretched" in the radial direct-
        // ion, it's time to double the cell size (incl the base). So check the
        // (Height / Width) @ phi=0:
        double rs = double(da / dw);

        if (rs > 2.0)
        {
          // Time to double the Elements Width:
          // New "np" will be for the next ellipse containing ~2 times less
          // arcs; re-calcuilate "dPhi" with the new "np":
          np   = max(np / 2, NF);
          nc   = 0;
          dPhi = TwoPi<double> / double(np);

          // Set "da" to make the next layer "square"-like:
          da = b * dPhi;
        }
        // Otherwise, keep the curr expanded "da"
      }
      else
        // We are in the "far-away" region; make square-like cells there:
        da = dw;
    }
    // In all cases:
    assert(NP % np == 0);

    // Increment "a", then re-calculate "b":
    a       += da;
    double r = (a - SqRt(Sqr(a) - Area(1.0))).Magnitude();
    assert(r > 0.0);
    b      = Len(0.5 * (1.0 / r - r));
    assert(IsPos(b) && b < a);

    if (verbose)
      cerr << ellipses.size() << " --> np=" << np << ", da="
           << da.Magnitude()  << ", a="     << a.Magnitude() << endl;
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
  // Generate the Rectangular and Triangular Mesh Elements:                  //
  //-------------------------------------------------------------------------//
  // The total number of Elements is determined dynamically:
  int MaxElems   = 99999999;
  fputs   ("NELEM= XXXXXXXX\n", f);
  long NElemsOff = ftell(f)  -  9;
  int  NElems    = 0;

  for (int e = 0; e <= NEllipses-2; ++e)
  {
    // First, process all Points  on this Ellipse:
    auto [pFrom, pTo] = ellipses[size_t(e)];

    // The number of Points on this Ellipse:
    int  npe          = pTo - pFrom + 1;
    assert(npe >= NF && npe % 2 == 0);

    // The next Ellipse may have either the same number of Points, or 2 times
    // less:
    auto [nFrom, nTo] = ellipses[size_t(e+1)];
    int  nn           = nTo - nFrom + 1;
    assert(nn >= 2 && nn % 2 == 0);

    bool   same = (npe == nn);
    assert(same   ||  npe == nn * 2);
    // In particular, if e <= nInner, we must have the "same" flag set, other-
    // wise we would not be able to complete the inner squares:
    assert(e > nInner || same);

    // If the "same" flag is set, then we use all Points of the Ellipse "e" as
    // the Rect vertices;  otherwise, use each 2nd point as Triangle vertices:
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
    assert(pTo+1 - npe == pFrom);

    if (same)
    {
      // Generate the Rectangles (Code=9):
      for (int i = pFrom; i <= pTo; ++i)
      {
        int i1 = (i <= pTo-1) ? (i+1) : pFrom;

        fprintf(f, "9 %d %d %d %d   %d\n",
                i, i1, Ort(i1), Ort(i), NElems++);
      }
    }
    else
    {
      // Generate the Triangles (Code=5):
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
    fprintf(f, "%% M = %.3lf, y^+ = %.3lf, nInner = %d, RER = %.3lf, NF = %d, "
               "ar0 = %.3lf\n",
            Moo, yPlus, nInner, RER, NF, ar0);

  fclose(f);
  return 0;
}

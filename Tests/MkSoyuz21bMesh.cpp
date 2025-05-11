// vim:ts=2:et
//===========================================================================//
//                        "Tests/MkSoyuz21bMesh.cpp":                        //
//  Construction of a 3D Mesh for AeroDynamic Simulation of "Soyuz-2.1b" LV  //
//                           using the "gmsh" API                            //
//===========================================================================//
#include <gmsh.h>

int main(/*int argc, char* argv[]*/)
{
	using namespace gmsh;
	initialize();

  model::add("Soyuz21b");
  model::occ::addCylinder(0,0,0, 0,0,25.0, 2.5);

	finalize();
	return 0;
}

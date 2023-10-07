// vim:ts=2:et
//===========================================================================//
//                              "Soyuz_21b.h":                               //
//          Mathematical Model of the "Soyuz-2.1b" Space Launcher            //
//===========================================================================//
#pragma  once
#include "Types.hpp"

namespace SpaceBallistics
{
	//=========================================================================//
	// "Soyuz_21b" Class:																											 //
	//=========================================================================//
	// "Soyuz_21b" is considered to be a 3-Stage rocket. This class is parameter-
	// ised by the "Payload"  which may in particular  include the 4th Stage, eg
	// "Fregat":
	//
	template<typename Payload>
  class Soyuz_21b
  {
  public:
    //=======================================================================//
    // Consts and Data Flds (in the "public" part):                          //
    //=======================================================================//
    Payload   const m_payload;

    // Fairing:

    //-----------------------------------------------------------------------//
    // 3rd stage:                                                            //
    //-----------------------------------------------------------------------//
    // NB: The following consts are actually static, but since they are comput-
    // ed by the Ctor, it is more convenient to make them member flds:
    Mass_kg   const m_stage3AftMass;      // Jetisonable Aft Mass
    MoI       const m_stage3EmptyMoI;     // With Aft Secrion

    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    Soyuz_21b(Payload const& a_payload);

    // Default and Copy Ctors are deleted, Dtor is trivial and auto-generated:
    Soyuz_21b()                 = delete;
    Soyuz_21b(Soyuz_21b const&) = delete;

  private:
    //-----------------------------------------------------------------------//
    // Moment of Inertia Computation for the Empty Rocket:                   //
    //-----------------------------------------------------------------------//
    MoI EmptyMoI() const;
  };
}

// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CoOrds/Bodies.h":                    //
// Nomenclature of Bodies endowed with Gravitational Field Models and COSes  //
//===========================================================================//
#pragma  once
#include <string_view>
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Body" Enum Class:                                                      //
  //=========================================================================//
  enum class Body: int
  {
    Sun     = 0,
    Mercury = 1,
    Venus   = 2,
    Earth   = 3,   // Earth alone, w/o the Moon!
    Mars    = 4,   // Mars    and its moons
    Jupiter = 5,   // Jupiter and its moons
    Saturn  = 6,   // Saturn  and its moons
    Uranus  = 7,   // Uranus  and its moons
    Neptune = 8,   // Neptune and its moons
    PlChB   = 9,   // Pluto-Charon System BaryCenter

    // We place the Moon after PlChB to preserve the classical numbering sequ-
    // ence for the Sun and Major Planets:
    Moon    = 10,

    // We also provide a "virtual Body"  EMB, the Earth-Moon System BaryCenter,
    // for compatibility with DE440T and for convenience of interplanetray traj-
    // ectores integration:
    EMB     = 11
  };

  //-------------------------------------------------------------------------//
  // Conversion To/From Strings:                                             //
  //-------------------------------------------------------------------------//
  // XXX: This could be done automatically, eg using "MagicEnum", but we do not
  // want to introduce extra external dependencies:
  //
  constexpr char const* ToString(Body a_body)
  {
    switch (a_body)
    {
      case Body::Sun     : return "Sun";
      case Body::Mercury : return "Mercury";
      case Body::Venus   : return "Venus";
      case Body::Earth   : return "Earth";
      case Body::Mars    : return "Mars";
      case Body::Jupiter : return "Jupiter";
      case Body::Saturn  : return "Saturn";
      case Body::Uranus  : return "Uranus";
      case Body::Neptune : return "Neptune";
      case Body::PlChB   : return "PlChB";
      case Body::Moon    : return "Moon";
      case Body::EMB     : return "EMB";
      default            : assert(false); return nullptr;
    }
  }

  constexpr Body ToBody(char const* a_body_str)
  {
    assert(a_body_str != nullptr);
    std::string_view bsv(a_body_str);
    return
      (bsv == "Sun")     ? Body::Sun     :
      (bsv == "Mercury") ? Body::Mercury :
      (bsv == "Venus")   ? Body::Venus   :
      (bsv == "Earth")   ? Body::Earth   :
      (bsv == "Mars")    ? Body::Mars    :
      (bsv == "Jupiter") ? Body::Jupiter :
      (bsv == "Saturn")  ? Body::Saturn  :
      (bsv == "Uranus")  ? Body::Uranus  :
      (bsv == "Neptune") ? Body::Neptune :
      (bsv == "PlChB")   ? Body::PlChB   :
      (bsv == "Moon")    ? Body::Moon    :
      (bsv == "EMB")     ? Body::EMB     :
      throw "Invalid BodyName";
  }
}
// End namespace SpaceBallistics

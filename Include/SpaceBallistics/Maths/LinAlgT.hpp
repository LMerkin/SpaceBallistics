// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/Maths/LinAlgT.hpp":                    //
//   Linear Operations on Heterogeneous Vectors (Represented by Tuples):     //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include <tuple>
#include <utility>

namespace SpaceBallistics
{
  //=========================================================================//
  // Implementation Functions:                                               //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "AddImpl":                                                              //
  //-------------------------------------------------------------------------//
  // Compute (x + c * y), where "x" and "y" are tuples (or any other container
  // types for which std::get<> is implemented):
  //
  template<typename TupleX, typename TupleY, typename C, size_t... Is>
  constexpr TupleX  AddImpl
  (
    TupleX const&  a_x,
    C              a_c,
    TupleY const&  a_y,
    std::index_sequence<Is...>
  )
  { return std::tuple{(std::get<Is>(a_x) + a_c * std::get<Is>(a_y)) ...}; }

  //-------------------------------------------------------------------------//
  // "MultImpl":                                                             //
  //-------------------------------------------------------------------------//
  // As above, but w/o the free term: returns (c * y):
  template<typename TupleY, typename C, size_t... Is>
  constexpr auto    MultImpl
  (
    C              a_c,
    TupleY const&  a_y,
    std::index_sequence<Is...>
  )
  { return std::tuple{(a_c * std::get<Is>(a_y)) ...}; }

  //-------------------------------------------------------------------------//
  // "AbsRelErrImpl":                                                        //
  //-------------------------------------------------------------------------//
  // For the tuples "e" and "y", computes the tuple of "double"s:
  //   |e| / (max(1.0_Y, |y|) :
  //
  template<typename TupleY, size_t... Is>
  constexpr auto AbsRelErrImpl
  (
    TupleY const&   a_e,
    TupleY const&   a_y,
    std::index_sequence<Is...>
  )
  {
    return std::tuple
    {
      double(Abs(std::get<Is>(a_e)) /
             Max(decltype(std::get<Is>(a_y))(1.0), Abs(std::get<Is>(a_y))))
      ...
    };
  }

  //-------------------------------------------------------------------------//
  // "FoldImpl":                                                             //
  //-------------------------------------------------------------------------//
  // Generic Fold on a Vector given by a tuple:        
  //
  template<typename T, typename F, typename... Ts>
  constexpr T FoldImpl
  (
    std::tuple<Ts...> const& a_tuple,
    T const&                 a_init,
    F const&                 a_f
  )
  {
    return std::apply
    (
      [&a_init, &a_f](auto const&... elems) -> T
      {
        T acc = a_init;
        (..., (acc = a_f(acc, elems))); // Left Fold
        return acc;
      },
      a_tuple
    );
  }

  //=========================================================================//
  // User-Level API:                                                         //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "Add": Computes (x + c * y):                                            //
  //-------------------------------------------------------------------------//
  template<typename TupleX, typename C, typename... TYs>
  constexpr TupleX  Add
  (
    TupleX const&             a_x,
    C                         a_c,
    std::tuple<TYs...> const& a_y
  )
  { return AddImpl(a_x, a_c, a_y, std::index_sequence_for<TYs...>{}); }

  //-------------------------------------------------------------------------//
  // "Mult": Computes (c * y):                                               //
  //-------------------------------------------------------------------------//
  template<typename C, typename... Ts>
  constexpr auto Mult
  (
    C                        a_c,
    std::tuple<Ts...> const& a_y
  )
  { return MultImpl(a_c, a_y, std::index_sequence_for<Ts...>{}); }

  //-------------------------------------------------------------------------//
  // "LinfNorm" (for tuples of "DimQ"s):                                     //
  //-------------------------------------------------------------------------//
  template<typename... Ts>
  constexpr double LinfNorm(std::tuple<Ts...> const& a_x)
  {
    return FoldImpl
    (
      a_x,
      0.0,
      [](double a_max, auto const& a_val) -> double
        { return std::max(a_max,
                          std::fabs(double(a_val / decltype(a_val)(1.0)))); }
    );
  }

  //-------------------------------------------------------------------------//
  // "MaxAbsRelError":                                                       //
  //-------------------------------------------------------------------------//
  template<typename... Ts>
  constexpr double MaxAbsRelError
  (
    std::tuple<Ts...> const& a_e,
    std::tuple<Ts...> const& a_y
  )
  {
    // Compute the tuple of AbsRelErrs (as "double"s):
    auto absRelErrs = AbsRelErrImpl(a_e, a_y, std::index_sequence_for<Ts...>{});

    // And the Linf norm of that tuple:
    return LinfNorm(absRelErrs);
  }
}
// End namespace SpaceBallistics

// vim:ts=2:et
//===========================================================================//
//                       "Src/PhysForces/DE440T.cpp":                        //
//===========================================================================//
#include "SpaceBallistics/PhysForces/DE440T.h"
#include <exception>

namespace SpaceBallistics
{
  //=========================================================================//
  // "SelfTest":                                                             //
  //=========================================================================//
  void DE440T::SelfTest()
  {
      // Verify the Temporal Continuity of DE440T Data:
      TDB expFrom = From;
      TDB expTo;    // Empty as yet

      for (int r = 0; r < NR; ++r)
      {
        Record const* rec =
          reinterpret_cast<Record const*>(&(s_data[r][0]));

        TDB from(rec->m_From);
        TDB to  (rec->m_To);
        expTo  = from + To_Time(32.0_day);

        if (UNLIKELY(from != expFrom || to != expTo))
          throw std::logic_error
                ("DE440T: Temporal Inconsistency in r=" + std::to_string(r));

        expFrom = to;
      }
      if (UNLIKELY(expTo != To))
        throw std::logic_error("DE440T: Temporal inconsistency at the end");
  }
}
// End namespace SpaceBallistics

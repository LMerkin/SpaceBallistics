// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/EmbeddedCOS.h":                 //
//       Embedded Co-Ords System (Bound to the Physical Axes of LV/SC)       //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/LVSC/LVSC.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "EmbeddedCOS" Class:                                                    //
  //=========================================================================//
  // Origin: A fixed physical point               in the LV or SC
  // Axes  : Fixed physical axes (eg of symmetry) of the LV or SC
  // NB    : This type is templatised by the LV/SC Type. It is intended to just
  //         stand for itself (no objs of this class are to be created):
  //
  template<LVSC LVSCKind>
  class EmbeddedCOS
  {
    EmbeddedCOS() = delete;	  // No objects construction at all!
  };
}
//End namespace SpaceBallistics

// vim:ts=2:et
//===========================================================================//
//                    "Src/LVSC/Soyuz-2.1b/Stage1.cpp":                      //
//                         Actual Code Generation                            //
//===========================================================================//
#define  STAGE1
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stages1-2.hpp"

namespace SpaceBallistics
{
  // StrapOn Boosters of Stage1: "Blocks B, V, G, D":
  template class Soyuz21b_Stage1_Booster<'B'>;
  template class Soyuz21b_Stage1_Booster<'V'>;
  template class Soyuz21b_Stage1_Booster<'G'>;
  template class Soyuz21b_Stage1_Booster<'D'>;
}
// End namespace SpaceBallistics

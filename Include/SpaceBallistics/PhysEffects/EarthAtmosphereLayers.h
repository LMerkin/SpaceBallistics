// vim:ts=2:et
//===========================================================================//
//          "SpaceBallistics/PhysEffects/EarthAtmosphereLayers.h":           //
//===========================================================================//
// XXX: The following fragment can be included into either
// "EarthAtmosphereModel.hpp" or "EarthAtmosphereModel.cpp";
// the latter is required with CLang which as yet (major version <= 20) does
// not support some C++23-mandated "constexpr" functions:
//
    LayerInfo
         (0.0_km, P0, T0, TempGrad(-6.5), 11.0_km),  // 0: TropoSpehere
    Layers[0].MkNextLayer(TempGrad( 0.0), 20.0_km),  // 1: TropoPause
    Layers[1].MkNextLayer(TempGrad( 1.0), 32.0_km),  // 2: LowerStratoSphere
    Layers[2].MkNextLayer(TempGrad( 2.8), 47.0_km),  // 3: UpperStratoSphere
    Layers[3].MkNextLayer(TempGrad( 0.0), 51.0_km),  // 4: StratoPause
    Layers[4].MkNextLayer(TempGrad(-2.8), 71.0_km),  // 5: LowerMesoSphere
    Layers[5].MkNextLayer(TempGrad(-2.0), 84.852_km) // 6: UpperMesoSphere


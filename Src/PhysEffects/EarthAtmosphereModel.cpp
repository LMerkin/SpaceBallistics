// vim:ts=2:et
//===========================================================================//
//               "Src/PhysEffects/EarthAtmosphereModel.cpp":                 //
//===========================================================================//
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"

namespace SpaceBallistics::EarthAtmosphereModel
{
  //=========================================================================//
  // "LowLayers" (z = 0 .. 93 km, but the arg is h < z):                     //
  //=========================================================================//
# ifdef   __clang__
# pragma  clang diagnostic push
# pragma  clang diagnostic ignored "-Wglobal-constructors"
# endif
  LowLayerInfo const LowLayers[NLowLayers]
  {
    LowLayerInfo
            (0.0_km, P0, T0, TempGrad(-6.5), 11.0_km), // 0: TropoSpehere
    LowLayers[0].MkNextLayer(TempGrad( 0.0), 20.0_km), // 1: TropoPause
    LowLayers[1].MkNextLayer(TempGrad( 1.0), 32.0_km), // 2: LowerStratoSphere
    LowLayers[2].MkNextLayer(TempGrad( 2.8), 47.0_km), // 3: UpperStratoSphere
    LowLayers[3].MkNextLayer(TempGrad( 0.0), 51.0_km), // 4: StratoPause
    LowLayers[4].MkNextLayer(TempGrad(-2.8), 71.0_km), // 5: LowerMesoSphere
    LowLayers[5].MkNextLayer(TempGrad(-2.0), 85.0_km), // 6: UpperMesoSphere
    LowLayers[6].MkNextLayer(TempGrad( 0.0), 93.0_km), // 7: MesoPause
  };
# ifdef   __clang__
# pragma  clang diagnostic pop
# endif

  //=========================================================================//
  // "UpperLayers1" (z = 93+ .. 300 km, step = 1 km):                        //
  //=========================================================================//
  UpperLayerInfo const UpperLayers1[NUpperLayers1]
  {
    {   93.0_km,  186.650_K, Pressure(1.07416e-1), Density(2.00484e-6) },
    {   94.0_km,  186.650_K, Pressure(8.99258e-2), Density(1.67840e-6) },
    {   95.0_km,  186.622_K, Pressure(7.52834e-2), Density(1.40510e-6) },
    {   96.0_km,  188.240_K, Pressure(6.30558e-2), Density(1.16617e-6) },
    {   97.0_km,  190.917_K, Pressure(5.29545e-2), Density(9.64450e-7) },
    {   98.0_km,  192.912_K, Pressure(4.45890e-2), Density(7.99924e-7) },
    {   99.0_km,  194.775_K, Pressure(3.76413e-2), Density(6.65312e-7) },
    {  100.0_km,  196.605_K, Pressure(3.18606e-2), Density(5.54951e-7) },
    {  101.0_km,  198.404_K, Pressure(2.70326e-2), Density(4.64111e-7) },
    {  102.0_km,  200.172_K, Pressure(2.29898e-2), Density(3.89128e-7) },
    {  103.0_km,  201.911_K, Pressure(1.95962e-2), Density(3.27067e-7) },
    {  104.0_km,  203.618_K, Pressure(1.67431e-2), Density(2.75614e-7) },
    {  105.0_km,  211.752_K, Pressure(1.43661e-2), Density(2.26166e-7) },
    {  106.0_km,  220.721_K, Pressure(1.24169e-2), Density(1.86514e-7) },
    {  107.0_km,  229.589_K, Pressure(1.08011e-2), Density(1.55121e-7) },
    {  108.0_km,  238.336_K, Pressure(9.45321e-3), Density(1.30060e-7) },
    {  109.0_km,  246.962_K, Pressure(8.32015e-3), Density(1.09860e-7) },
    {  110.0_km,  255.487_K, Pressure(7.35900e-3), Density(9.34035e-8) },
    {  111.0_km,  263.891_K, Pressure(6.54011e-3), Density(7.99116e-8) },
    {  112.0_km,  272.184_K, Pressure(5.83750e-3), Density(6.87683e-8) },
    {  113.0_km,  280.356_K, Pressure(5.23190e-3), Density(5.94575e-8) },
    {  114.0_km,  288.427_K, Pressure(4.70615e-3), Density(5.17245e-8) },
    {  115.0_km,  296.377_K, Pressure(4.24849e-3), Density(4.51814e-8) },
    {  116.0_km,  304.207_K, Pressure(3.84826e-3), Density(3.96418e-8) },
    {  117.0_km,  311.926_K, Pressure(3.49640e-3), Density(3.49222e-8) },
    {  118.0_km,  319.535_K, Pressure(3.18584e-3), Density(3.08814e-8) },
    {  119.0_km,  327.033_K, Pressure(2.91073e-3), Density(2.74061e-8) },
    {  120.0_km,  334.417_K, Pressure(2.66618e-3), Density(2.44041e-8) },
    {  121.0_km,  345.68_K,  Pressure(2.50215e-3), Density(2.20814e-8) },
    {  122.0_km,  356.94_K,  Pressure(2.36291e-3), Density(2.01267e-8) },
    {  123.0_km,  368.20_K,  Pressure(2.22898e-3), Density(1.83438e-8) },
    {  124.0_km,  379.46_K,  Pressure(2.10075e-3), Density(1.67199e-8) },
    {  125.0_km,  390.72_K,  Pressure(1.97851e-3), Density(1.52430e-8) },
    {  126.0_km,  401.97_K,  Pressure(1.86245e-3), Density(1.39015e-8) },
    {  127.0_km,  413.23_K,  Pressure(1.75266e-3), Density(1.26846e-8) },
    {  128.0_km,  424.49_K,  Pressure(1.64917e-3), Density(1.15819e-8) },
    {  129.0_km,  435.75_K,  Pressure(1.55191e-3), Density(1.05836e-8) },
    {  130.0_km,  447.01_K,  Pressure(1.46076e-3), Density(9.68064e-9) },
    {  131.0_km,  458.27_K,  Pressure(1.37554e-3), Density(8.86428e-9) },
    {  132.0_km,  469.53_K,  Pressure(1.29601e-3), Density(8.12642e-9) },
    {  133.0_km,  480.79_K,  Pressure(1.22190e-3), Density(7.45948e-9) },
    {  134.0_km,  492.05_K,  Pressure(1.15289e-3), Density(6.85639e-9) },
    {  135.0_km,  503.31_K,  Pressure(1.08865e-3), Density(6.31060e-9) },
    {  136.0_km,  514.56_K,  Pressure(1.02883e-3), Density(5.81605e-9) },
    {  137.0_km,  525.82_K,  Pressure(9.73053e-4), Density(5.36717e-9) },
    {  138.0_km,  537.08_K,  Pressure(9.20958e-4), Density(4.95889e-9) },
    {  139.0_km,  548.34_K,  Pressure(8.72181e-4), Density(4.58659e-9) },
    {  140.0_km,  559.60_K,  Pressure(8.26375e-4), Density(4.24614e-9) },
    {  141.0_km,  566.40_K,  Pressure(7.77095e-4), Density(3.93386e-9) },
    {  142.0_km,  573.20_K,  Pressure(7.31027e-4), Density(3.64654e-9) },
    {  143.0_km,  580.00_K,  Pressure(6.87814e-4), Density(3.38138e-9) },
    {  144.0_km,  586.80_K,  Pressure(6.47162e-4), Density(3.13605e-9) },
    {  145.0_km,  593.60_K,  Pressure(6.08841e-4), Density(2.90864e-9) },
    {  146.0_km,  600.40_K,  Pressure(5.72689e-4), Density(2.69768e-9) },
    {  147.0_km,  607.20_K,  Pressure(5.38618e-4), Density(2.50210e-9) },
    {  148.0_km,  614.00_K,  Pressure(5.06620e-4), Density(2.32127e-9) },
    {  149.0_km,  620.80_K,  Pressure(4.76771e-4), Density(2.15494e-9) },
    {  150.0_km,  627.60_K,  Pressure(4.49233e-4), Density(2.00329e-9) },
    {  151.0_km,  634.40_K,  Pressure(4.31065e-4), Density(1.89681e-9) },
    {  152.0_km,  641.20_K,  Pressure(4.14031e-4), Density(1.79797e-9) },
    {  153.0_km,  648.00_K,  Pressure(3.97769e-4), Density(1.70494e-9) },
    {  154.0_km,  654.80_K,  Pressure(3.82257e-4), Density(1.61741e-9) },
    {  155.0_km,  661.60_K,  Pressure(3.67473e-4), Density(1.53510e-9) },
    {  156.0_km,  668.40_K,  Pressure(3.53393e-4), Density(1.45770e-9) },
    {  157.0_km,  675.20_K,  Pressure(3.39992e-4), Density(1.38496e-9) },
    {  158.0_km,  682.00_K,  Pressure(3.27247e-4), Density(1.31660e-9) },
    {  159.0_km,  688.80_K,  Pressure(3.15131e-4), Density(1.25237e-9) },
    {  160.0_km,  695.60_K,  Pressure(3.03620e-4), Density(1.19204e-9) },
    {  161.0_km,  699.57_K,  Pressure(2.91509e-4), Density(1.13535e-9) },
    {  162.0_km,  703.54_K,  Pressure(2.80055e-4), Density(1.08210e-9) },
    {  163.0_km,  707.51_K,  Pressure(2.69225e-4), Density(1.03206e-9) },
    {  164.0_km,  711.48_K,  Pressure(2.58982e-4), Density(9.85037e-10) },
    {  165.0_km,  715.49_K,  Pressure(2.49295e-4), Density(9.40825e-10) },
    {  166.0_km,  719.42_K,  Pressure(2.40129e-4), Density(8.99240e-10) },
    {  167.0_km,  723.39_K,  Pressure(2.31452e-4), Density(8.60103e-10) },
    {  168.0_km,  727.36_K,  Pressure(2.23234e-4), Density(8.23245e-10) },
    {  169.0_km,  731.33_K,  Pressure(2.15444e-4), Density(7.88506e-10) },
    {  170.0_km,  735.30_K,  Pressure(2.08053e-4), Density(7.55732e-10) },
    {  171.0_km,  739.27_K,  Pressure(2.01034e-4), Density(7.24781e-10) },
    {  172.0_km,  743.24_K,  Pressure(1.94359e-4), Density(6.95515e-10) },
    {  173.0_km,  747.21_K,  Pressure(1.88004e-4), Density(6.67809e-10) },
    {  174.0_km,  751.18_K,  Pressure(1.81943e-4), Density(6.41542e-10) },
    {  175.0_km,  755.15_K,  Pressure(1.76153e-4), Density(6.16603e-10) },
    {  176.0_km,  759.12_K,  Pressure(1.70614e-4), Density(5.92889e-10) },
    {  177.0_km,  763.09_K,  Pressure(1.65304e-4), Density(5.70303e-10) },
    {  178.0_km,  767.06_K,  Pressure(1.60206e-4), Density(5.48759e-10) },
    {  179.0_km,  771.03_K,  Pressure(1.55301e-4), Density(5.28175e-10) },
    {  180.0_km,  775.00_K,  Pressure(1.50574e-4), Density(5.08478e-10) },
    {  181.0_km,  778.97_K,  Pressure(1.46011e-4), Density(4.89604e-10) },
    {  182.0_km,  782.94_K,  Pressure(1.41600e-4), Density(4.71493e-10) },
    {  183.0_km,  786.91_K,  Pressure(1.37329e-4), Density(4.54096e-10) },
    {  184.0_km,  790.88_K,  Pressure(1.33190e-4), Density(4.37368e-10) },
    {  185.0_km,  794.85_K,  Pressure(1.29177e-4), Density(4.21273e-10) },
    {  186.0_km,  798.82_K,  Pressure(1.25282e-4), Density(4.05781e-10) },
    {  187.0_km,  802.79_K,  Pressure(1.21505e-4), Density(3.90870e-10) },
    {  188.0_km,  806.76_K,  Pressure(1.17842e-4), Density(3.76524e-10) },
    {  189.0_km,  810.73_K,  Pressure(1.14294e-4), Density(3.62734e-10) },
    {  190.0_km,  814.70_K,  Pressure(1.10866e-4), Density(3.49499e-10) },
    {  191.0_km,  818.67_K,  Pressure(1.07561e-4), Density(3.36823e-10) },
    {  192.0_km,  822.64_K,  Pressure(1.04386e-4), Density(3.24717e-10) },
    {  193.0_km,  826.61_K,  Pressure(1.01351e-4), Density(3.13199e-10) },
    {  194.0_km,  830.58_K,  Pressure(9.84678e-5), Density(3.02294e-10) },
    {  195.0_km,  834.55_K,  Pressure(9.57500e-5), Density(2.92032e-10) },
    {  196.0_km,  838.52_K,  Pressure(9.32141e-5), Density(2.82451e-10) },
    {  197.0_km,  842.49_K,  Pressure(9.08788e-5), Density(2.73595e-10) },
    {  198.0_km,  846.46_K,  Pressure(8.87655e-5), Density(2.65513e-10) },
    {  199.0_km,  850.43_K,  Pressure(8.68979e-5), Density(2.58262e-10) },
    {  200.0_km,  854.40_K,  Pressure(8.53026e-5), Density(2.51904e-10) },
    {  201.0_km,  856.15_K,  Pressure(8.23175e-5), Density(2.42171e-10) },
    {  202.0_km,  857.90_K,  Pressure(8.00431e-5), Density(2.34594e-10) },
    {  203.0_km,  859.65_K,  Pressure(7.78409e-5), Density(2.27283e-10) },
    {  204.0_km,  861.40_K,  Pressure(7.57087e-5), Density(2.20229e-10) },
    {  205.0_km,  863.15_K,  Pressure(7.36446e-5), Density(2.13424e-10) },
    {  206.0_km,  864.90_K,  Pressure(7.16466e-5), Density(2.06859e-10) },
    {  207.0_km,  866.65_K,  Pressure(6.97127e-5), Density(2.00526e-10) },
    {  208.0_km,  868.40_K,  Pressure(6.78409e-5), Density(1.94416e-10) },
    {  209.0_km,  870.15_K,  Pressure(6.60295e-5), Density(1.88523e-10) },
    {  210.0_km,  871.90_K,  Pressure(6.42764e-5), Density(1.82838e-10) },
    {  211.0_km,  873.65_K,  Pressure(6.25800e-5), Density(1.77353e-10) },
    {  212.0_km,  875.40_K,  Pressure(6.09383e-5), Density(1.72062e-10) },
    {  213.0_km,  877.15_K,  Pressure(5.93496e-5), Density(1.66958e-10) },
    {  214.0_km,  878.90_K,  Pressure(5.78122e-5), Density(1.62033e-10) },
    {  215.0_km,  880.65_K,  Pressure(5.63242e-5), Density(1.57280e-10) },
    {  216.0_km,  882.40_K,  Pressure(5.48842e-5), Density(1.52694e-10) },
    {  217.0_km,  884.15_K,  Pressure(5.34903e-5), Density(1.48268e-10) },
    {  218.0_km,  885.90_K,  Pressure(5.21409e-5), Density(1.43996e-10) },
    {  219.0_km,  887.65_K,  Pressure(5.08346e-5), Density(1.39871e-10) },
    {  220.0_km,  889.40_K,  Pressure(4.95697e-5), Density(1.35889e-10) },
    {  221.0_km,  891.15_K,  Pressure(4.83448e-5), Density(1.32043e-10) },
    {  222.0_km,  892.90_K,  Pressure(4.71583e-5), Density(1.28328e-10) },
    {  223.0_km,  894.65_K,  Pressure(4.60088e-5), Density(1.24738e-10) },
    {  224.0_km,  896.40_K,  Pressure(4.48950e-5), Density(1.21269e-10) },
    {  225.0_km,  898.15_K,  Pressure(4.38153e-5), Density(1.17916e-10) },
    {  226.0_km,  899.90_K,  Pressure(4.27686e-5), Density(1.14674e-10) },
    {  227.0_km,  901.65_K,  Pressure(4.17535e-5), Density(1.11538e-10) },
    {  228.0_km,  903.40_K,  Pressure(4.07687e-5), Density(1.08504e-10) },
    {  229.0_km,  905.15_K,  Pressure(3.98131e-5), Density(1.05567e-10) },
    {  230.0_km,  906.90_K,  Pressure(3.88854e-5), Density(1.02724e-10) },
    {  231.0_km,  908.65_K,  Pressure(3.79845e-5), Density(9.99709e-11) },
    {  232.0_km,  910.40_K,  Pressure(3.71094e-5), Density(9.73033e-11) },
    {  233.0_km,  912.15_K,  Pressure(3.62589e-5), Density(9.47179e-11) },
    {  234.0_km,  913.90_K,  Pressure(3.54321e-5), Density(9.22112e-11) },
    {  235.0_km,  915.65_K,  Pressure(3.46279e-5), Density(8.97799e-11) },
    {  236.0_km,  917.40_K,  Pressure(3.38455e-5), Density(8.74209e-11) },
    {  237.0_km,  919.15_K,  Pressure(3.30839e-5), Density(8.51310e-11) },
    {  238.0_km,  920.90_K,  Pressure(3.23424e-5), Density(8.29076e-11) },
    {  239.0_km,  922.65_K,  Pressure(3.16200e-5), Density(8.07479e-11) },
    {  240.0_km,  924.40_K,  Pressure(3.09161e-5), Density(7.86493e-11) },
    {  241.0_km,  926.15_K,  Pressure(3.02299e-5), Density(7.66094e-11) },
    {  242.0_km,  927.90_K,  Pressure(2.95607e-5), Density(7.46260e-11) },
    {  243.0_km,  929.65_K,  Pressure(2.89080e-5), Density(7.26970e-11) },
    {  244.0_km,  931.40_K,  Pressure(2.82711e-5), Density(7.08203e-11) },
    {  245.0_km,  933.15_K,  Pressure(2.76496e-5), Density(6.89942e-11) },
    {  246.0_km,  934.90_K,  Pressure(2.70429e-5), Density(6.72170e-11) },
    {  247.0_km,  936.65_K,  Pressure(2.64506e-5), Density(6.54871e-11) },
    {  248.0_km,  938.40_K,  Pressure(2.58723e-5), Density(6.38031e-11) },
    {  249.0_km,  940.15_K,  Pressure(2.53076e-5), Density(6.21637e-11) },
    {  250.0_km,  941.90_K,  Pressure(2.47564e-5), Density(6.05679e-11) },
    {  251.0_km,  942.47_K,  Pressure(2.42279e-5), Density(5.91295e-11) },
    {  252.0_km,  943.04_K,  Pressure(2.36864e-5), Density(5.76710e-11) },
    {  253.0_km,  943.61_K,  Pressure(2.31584e-5), Density(5.62529e-11) },
    {  254.0_km,  944.18_K,  Pressure(2.26439e-5), Density(5.48740e-11) },
    {  255.0_km,  944.75_K,  Pressure(2.21424e-5), Density(5.35333e-11) },
    {  256.0_km,  945.32_K,  Pressure(2.16535e-5), Density(5.22299e-11) },
    {  257.0_km,  945.89_K,  Pressure(2.11771e-5), Density(5.09625e-11) },
    {  258.0_km,  946.46_K,  Pressure(2.07127e-5), Density(4.97302e-11) },
    {  259.0_km,  947.03_K,  Pressure(2.02602e-5), Density(4.85321e-11) },
    {  260.0_km,  947.60_K,  Pressure(1.98191e-5), Density(4.73671e-11) },
    {  261.0_km,  948.17_K,  Pressure(1.93892e-5), Density(4.62343e-11) },
    {  262.0_km,  948.74_K,  Pressure(1.89701e-5), Density(4.51328e-11) },
    {  263.0_km,  949.31_K,  Pressure(1.85617e-5), Density(4.40617e-11) },
    {  264.0_km,  949.88_K,  Pressure(1.81636e-5), Density(4.30201e-11) },
    {  265.0_km,  950.45_K,  Pressure(1.77756e-5), Density(4.20071e-11) },
    {  266.0_km,  951.02_K,  Pressure(1.73973e-5), Density(4.10218e-11) },
    {  267.0_km,  951.59_K,  Pressure(1.70286e-5), Density(4.00635e-11) },
    {  268.0_km,  952.16_K,  Pressure(1.66691e-5), Density(3.91313e-11) },
    {  269.0_km,  952.73_K,  Pressure(1.63185e-5), Density(3.82245e-11) },
    {  270.0_km,  953.30_K,  Pressure(1.59767e-5), Density(3.73422e-11) },
    {  271.0_km,  953.87_K,  Pressure(1.56434e-5), Density(3.64837e-11) },
    {  272.0_km,  954.44_K,  Pressure(1.53183e-5), Density(3.56482e-11) },
    {  273.0_km,  955.01_K,  Pressure(1.50012e-5), Density(3.48351e-11) },
    {  274.0_km,  955.58_K,  Pressure(1.46919e-5), Density(3.40437e-11) },
    {  275.0_km,  956.15_K,  Pressure(1.43901e-5), Density(3.32732e-11) },
    {  276.0_km,  956.72_K,  Pressure(1.40956e-5), Density(3.25230e-11) },
    {  277.0_km,  957.29_K,  Pressure(1.38081e-5), Density(3.17925e-11) },
    {  278.0_km,  957.86_K,  Pressure(1.35276e-5), Density(3.10810e-11) },
    {  279.0_km,  958.43_K,  Pressure(1.32537e-5), Density(3.03879e-11) },
    {  280.0_km,  959.00_K,  Pressure(1.29862e-5), Density(2.97125e-11) },
    {  281.0_km,  959.57_K,  Pressure(1.27250e-5), Density(2.90544e-11) },
    {  282.0_km,  960.14_K,  Pressure(1.24699e-5), Density(2.84130e-11) },
    {  283.0_km,  960.71_K,  Pressure(1.22205e-5), Density(2.77876e-11) },
    {  284.0_km,  961.28_K,  Pressure(1.19769e-5), Density(2.71778e-11) },
    {  285.0_km,  961.85_K,  Pressure(1.17388e-5), Density(2.65831e-11) },
    {  286.0_km,  962.42_K,  Pressure(1.15060e-5), Density(2.60029e-11) },
    {  287.0_km,  962.99_K,  Pressure(1.12783e-5), Density(2.54368e-11) },
    {  288.0_km,  963.56_K,  Pressure(1.10556e-5), Density(2.48843e-11) },
    {  289.0_km,  964.13_K,  Pressure(1.08377e-5), Density(2.43449e-11) },
    {  290.0_km,  964.70_K,  Pressure(1.06244e-5), Density(2.38182e-11) },
    {  291.0_km,  965.27_K,  Pressure(1.04157e-5), Density(2.33038e-11) },
    {  292.0_km,  965.84_K,  Pressure(1.02112e-5), Density(2.28013e-11) },
    {  293.0_km,  966.41_K,  Pressure(1.00110e-5), Density(2.23102e-11) },
    {  294.0_km,  966.98_K,  Pressure(9.81493e-6), Density(2.18302e-11) },
    {  295.0_km,  967.55_K,  Pressure(9.62274e-6), Density(2.13610e-11) },
    {  296.0_km,  968.12_K,  Pressure(9.43437e-6), Density(2.09021e-11) },
    {  297.0_km,  968.69_K,  Pressure(9.24969e-6), Density(2.04532e-11) },
    {  298.0_km,  969.26_K,  Pressure(9.06858e-6), Density(2.00141e-11) },
    {  299.0_km,  969.83_K,  Pressure(8.89094e-6), Density(1.95843e-11) },
    {  300.0_km,  970.40_K,  Pressure(8.71665e-6), Density(1.91637e-11) }
  };

  //=========================================================================//
  // "UpperLayers2" (z = 300+ .. 500 km, step = 2 km):                       //
  //=========================================================================//
  UpperLayerInfo const UpperLayers2[NUpperLayers2]
  {
    {  300.0_km,  970.40_K,  Pressure(8.71665e-6), Density(1.91637e-11) },
    {  302.0_km,  971.54_K,  Pressure(8.37773e-6), Density(1.83486e-11) },
    {  304.0_km,  972.68_K,  Pressure(8.05108e-6), Density(1.75668e-11) },
    {  306.0_km,  973.82_K,  Pressure(7.73604e-6), Density(1.68164e-11) },
    {  308.0_km,  974.96_K,  Pressure(7.43202e-6), Density(1.60957e-11) },
    {  310.0_km,  976.10_K,  Pressure(7.13852e-6), Density(1.54032e-11) },
    {  312.0_km,  977.24_K,  Pressure(6.85513e-6), Density(1.47378e-11) },
    {  314.0_km,  978.38_K,  Pressure(6.58153e-6), Density(1.40984e-11) },
    {  316.0_km,  979.52_K,  Pressure(6.31746e-6), Density(1.34842e-11) },
    {  318.0_km,  980.66_K,  Pressure(6.06278e-6), Density(1.28945e-11) },
    {  320.0_km,  981.80_K,  Pressure(5.81741e-6), Density(1.23289e-11) },
    {  322.0_km,  982.94_K,  Pressure(5.58438e-6), Density(1.17871e-11) },
    {  324.0_km,  984.08_K,  Pressure(5.35478e-6), Density(1.12692e-11) },
    {  326.0_km,  984.80_K,  Pressure(5.13564e-6), Density(1.07752e-11) },
    {  328.0_km,  985.10_K,  Pressure(4.92450e-6), Density(1.03055e-11) },
    {  330.0_km,  985.40_K,  Pressure(4.72399e-6), Density(9.86051e-12) },
    {  332.0_km,  985.70_K,  Pressure(4.53455e-6), Density(9.44101e-12) },
    {  334.0_km,  986.00_K,  Pressure(4.35670e-6), Density(9.04785e-12) },
    {  336.0_km,  986.30_K,  Pressure(4.19106e-6), Density(8.68207e-12) },
    {  338.0_km,  986.60_K,  Pressure(4.03833e-6), Density(8.34490e-12) },
    {  340.0_km,  986.90_K,  Pressure(3.89929e-6), Density(8.03773e-12) },
    {  342.0_km,  987.20_K,  Pressure(3.77483e-6), Density(7.76214e-12) },
    {  344.0_km,  987.50_K,  Pressure(3.66592e-6), Density(7.51986e-12) },
    {  346.0_km,  987.80_K,  Pressure(3.57361e-6), Density(7.31282e-12) },
    {  348.0_km,  988.10_K,  Pressure(3.49906e-6), Density(7.14309e-12) },
    {  350.0_km,  988.40_K,  Pressure(3.44350e-6), Density(7.01292e-12) },
    {  352.0_km,  988.70_K,  Pressure(3.32323e-6), Density(6.75194e-12) },
    {  354.0_km,  989.00_K,  Pressure(3.20680e-6), Density(6.50003e-12) },
    {  356.0_km,  989.30_K,  Pressure(3.09469e-6), Density(6.25809e-12) },
    {  358.0_km,  989.60_K,  Pressure(2.98675e-6), Density(6.02573e-12) },
    {  360.0_km,  989.90_K,  Pressure(2.88283e-6), Density(5.80258e-12) },
    {  362.0_km,  990.20_K,  Pressure(2.78280e-6), Density(5.58829e-12) },
    {  364.0_km,  990.50_K,  Pressure(2.68650e-6), Density(5.38249e-12) },
    {  366.0_km,  990.80_K,  Pressure(2.59381e-6), Density(5.18486e-12) },
    {  368.0_km,  991.10_K,  Pressure(2.50458e-6), Density(4.99506e-12) },
    {  370.0_km,  991.40_K,  Pressure(2.41869e-6), Density(4.81277e-12) },
    {  372.0_km,  991.70_K,  Pressure(2.33601e-6), Density(4.63769e-12) },
    {  374.0_km,  992.00_K,  Pressure(2.25642e-6), Density(4.46952e-12) },
    {  376.0_km,  992.30_K,  Pressure(2.17979e-6), Density(4.30797e-12) },
    {  378.0_km,  992.60_K,  Pressure(2.10602e-6), Density(4.15277e-12) },
    {  380.0_km,  992.90_K,  Pressure(2.03498e-6), Density(4.00363e-12) },
    {  382.0_km,  993.20_K,  Pressure(1.96657e-6), Density(3.86031e-12) },
    {  384.0_km,  993.50_K,  Pressure(1.90068e-6), Density(3.72256e-12) },
    {  386.0_km,  993.80_K,  Pressure(1.83720e-6), Density(3.59013e-12) },
    {  388.0_km,  994.10_K,  Pressure(1.77605e-6), Density(3.46280e-12) },
    {  390.0_km,  994.40_K,  Pressure(1.71711e-6), Density(3.34034e-12) },
    {  392.0_km,  994.70_K,  Pressure(1.66031e-6), Density(3.22253e-12) },
    {  394.0_km,  995.00_K,  Pressure(1.60555e-6), Density(3.10918e-12) },
    {  396.0_km,  995.30_K,  Pressure(1.55274e-6), Density(3.00010e-12) },
    {  398.0_km,  995.60_K,  Pressure(1.50180e-6), Density(2.89508e-12) },
    {  400.0_km,  995.90_K,  Pressure(1.45265e-6), Density(2.79396e-12) },
    {  402.0_km,  995.94_K,  Pressure(1.40486e-6), Density(2.69942e-12) },
    {  404.0_km,  995.98_K,  Pressure(1.35874e-6), Density(2.60748e-12) },
    {  406.0_km,  996.02_K,  Pressure(1.31422e-6), Density(2.51874e-12) },
    {  408.0_km,  996.06_K,  Pressure(1.27124e-6), Density(2.43308e-12) },
    {  410.0_km,  996.10_K,  Pressure(1.22973e-6), Density(2.35038e-12) },
    {  412.0_km,  996.14_K,  Pressure(1.18964e-6), Density(2.27052e-12) },
    {  414.0_km,  996.18_K,  Pressure(1.15090e-6), Density(2.19340e-12) },
    {  416.0_km,  996.22_K,  Pressure(1.11348e-6), Density(2.11891e-12) },
    {  418.0_km,  996.26_K,  Pressure(1.07731e-6), Density(2.04695e-12) },
    {  420.0_km,  996.30_K,  Pressure(1.04236e-6), Density(1.97745e-12) },
    {  422.0_km,  996.34_K,  Pressure(1.00859e-6), Density(1.91031e-12) },
    {  424.0_km,  996.38_K,  Pressure(9.75957e-7), Density(1.84547e-12) },
    {  426.0_km,  996.42_K,  Pressure(9.44428e-7), Density(1.78285e-12) },
    {  428.0_km,  996.46_K,  Pressure(9.13974e-7), Density(1.72241e-12) },
    {  430.0_km,  996.50_K,  Pressure(8.84568e-7), Density(1.66408e-12) },
    {  432.0_km,  996.54_K,  Pressure(8.56189e-7), Density(1.60781e-12) },
    {  434.0_km,  996.58_K,  Pressure(8.28816e-7), Density(1.55357e-12) },
    {  436.0_km,  996.62_K,  Pressure(8.02435e-7), Density(1.50131e-12) },
    {  438.0_km,  996.66_K,  Pressure(7.77033e-7), Density(1.45103e-12) },
    {  440.0_km,  996.70_K,  Pressure(7.52601e-7), Density(1.40268e-12) },
    {  442.0_km,  996.74_K,  Pressure(7.29136e-7), Density(1.35626e-12) },
    {  444.0_km,  996.78_K,  Pressure(7.06635e-7), Density(1.31176e-12) },
    {  446.0_km,  996.82_K,  Pressure(6.85102e-7), Density(1.26918e-12) },
    {  448.0_km,  996.86_K,  Pressure(6.64542e-7), Density(1.22852e-12) },
    {  450.0_km,  996.90_K,  Pressure(6.44964e-7), Density(1.18980e-12) },
    {  452.0_km,  996.94_K,  Pressure(6.24047e-7), Density(1.14872e-12) },
    {  454.0_km,  996.98_K,  Pressure(6.05019e-7), Density(1.11124e-12) },
    {  456.0_km,  997.02_K,  Pressure(5.86581e-7), Density(1.07496e-12) },
    {  458.0_km,  997.06_K,  Pressure(5.68718e-7), Density(1.03985e-12) },
    {  460.0_km,  997.10_K,  Pressure(5.51415e-7), Density(1.00587e-12) },
    {  462.0_km,  997.14_K,  Pressure(5.34656e-7), Density(9.73002e-13) },
    {  464.0_km,  997.18_K,  Pressure(5.18428e-7), Density(9.41207e-13) },
    {  466.0_km,  997.22_K,  Pressure(5.02715e-7), Density(9.10457e-13) },
    {  468.0_km,  997.26_K,  Pressure(4.87503e-7), Density(8.80722e-13) },
    {  470.0_km,  997.30_K,  Pressure(4.72778e-7), Density(8.51973e-13) },
    {  472.0_km,  997.34_K,  Pressure(4.58527e-7), Density(8.24181e-13) },
    {  474.0_km,  997.38_K,  Pressure(4.44736e-7), Density(7.97319e-13) },
    {  476.0_km,  997.42_K,  Pressure(4.31391e-7), Density(7.71357e-13) },
    {  478.0_km,  997.46_K,  Pressure(4.18480e-7), Density(7.46270e-13) },
    {  480.0_km,  997.50_K,  Pressure(4.05990e-7), Density(7.22030e-13) },
    {  482.0_km,  997.54_K,  Pressure(3.93907e-7), Density(6.98612e-13) },
    {  484.0_km,  997.58_K,  Pressure(3.82220e-7), Density(6.75989e-13) },
    {  486.0_km,  997.62_K,  Pressure(3.70917e-7), Density(6.54137e-13) },
    {  488.0_km,  997.66_K,  Pressure(3.59985e-7), Density(6.33031e-13) },
    {  490.0_km,  997.70_K,  Pressure(3.49413e-7), Density(6.12647e-13) },
    {  492.0_km,  997.74_K,  Pressure(3.39190e-7), Density(5.92960e-13) },
    {  494.0_km,  997.78_K,  Pressure(3.29304e-7), Density(5.73949e-13) },
    {  496.0_km,  997.82_K,  Pressure(3.19744e-7), Density(5.55590e-13) },
    {  498.0_km,  997.86_K,  Pressure(3.10501e-7), Density(5.37862e-13) },
    {  500.0_km,  997.90_K,  Pressure(3.01562e-7), Density(5.20743e-13) }
  };

  //=========================================================================//
  // "UpperLayers5" (z = 500+ .. 1200 km, step = 5 km):                      //
  //=========================================================================//
  UpperLayerInfo const UpperLayers5[NUpperLayers5]
  {
    {  500.0_km,  997.90_K,  Pressure(3.01562e-7), Density(5.20743e-13) },
    {  505.0_km,  998.00_K,  Pressure(2.80484e-7), Density(4.80474e-13) },
    {  510.0_km,  998.10_K,  Pressure(2.61098e-7), Density(4.43567e-13) },
    {  515.0_km,  998.20_K,  Pressure(2.43258e-7), Density(4.09731e-13) },
    {  520.0_km,  998.30_K,  Pressure(2.26830e-7), Density(3.78691e-13) },
    {  525.0_km,  998.40_K,  Pressure(2.11689e-7), Density(3.50196e-13) },
    {  530.0_km,  998.50_K,  Pressure(1.97719e-7), Density(3.24011e-13) },
    {  535.0_km,  998.60_K,  Pressure(1.84812e-7), Density(2.99926e-13) },
    {  540.0_km,  998.70_K,  Pressure(1.72871e-7), Density(2.77744e-13) },
    {  545.0_km,  998.80_K,  Pressure(1.61807e-7), Density(2.57292e-13) },
    {  550.0_km,  998.90_K,  Pressure(1.51543e-7), Density(2.38414e-13) },
    {  555.0_km,  999.00_K,  Pressure(1.42008e-7), Density(2.20971e-13) },
    {  560.0_km,  999.10_K,  Pressure(1.33142e-7), Density(2.04842e-13) },
    {  565.0_km,  999.20_K,  Pressure(1.24894e-7), Density(1.89925e-13) },
    {  570.0_km,  999.30_K,  Pressure(1.17222e-7), Density(1.76131e-13) },
    {  575.0_km,  999.40_K,  Pressure(1.10095e-7), Density(1.63389e-13) },
    {  580.0_km,  999.50_K,  Pressure(1.03490e-7), Density(1.51643e-13) },
    {  585.0_km,  999.60_K,  Pressure(9.73924e-8), Density(1.40850e-13) },
    {  590.0_km,  999.70_K,  Pressure(9.17988e-8), Density(1.30981e-13) },
    {  595.0_km,  999.80_K,  Pressure(8.67143e-8), Density(1.22020e-13) },
    {  600.0_km,  999.90_K,  Pressure(8.21535e-8), Density(1.13960e-13) },
    {  605.0_km,  999.90_K,  Pressure(7.76907e-8), Density(1.06205e-13) },
    {  610.0_km,  999.91_K,  Pressure(7.35576e-8), Density(9.90517e-14) },
    {  615.0_km,  999.91_K,  Pressure(6.96888e-8), Density(9.23975e-14) },
    {  620.0_km,  999.91_K,  Pressure(6.60695e-8), Density(8.62103e-14) },
    {  625.0_km,  999.91_K,  Pressure(6.26855e-8), Density(8.04593e-14) },
    {  630.0_km,  999.92_K,  Pressure(5.95228e-8), Density(7.51152e-14) },
    {  635.0_km,  999.92_K,  Pressure(5.65682e-8), Density(7.01498e-14) },
    {  640.0_km,  999.92_K,  Pressure(5.38089e-8), Density(6.55367e-14) },
    {  645.0_km,  999.92_K,  Pressure(5.12324e-8), Density(6.12504e-14) },
    {  650.0_km,  999.93_K,  Pressure(4.88269e-8), Density(5.72670e-14) },
    {  655.0_km,  999.93_K,  Pressure(4.65811e-8), Density(5.35490e-14) },
    {  660.0_km,  999.93_K,  Pressure(4.44839e-8), Density(5.01275e-14) },
    {  665.0_km,  999.93_K,  Pressure(4.25250e-8), Density(4.69735e-14) },
    {  670.0_km,  999.94_K,  Pressure(4.06945e-8), Density(4.40643e-14) },
    {  675.0_km,  999.94_K,  Pressure(3.89828e-8), Density(4.13788e-14) },
    {  680.0_km,  999.94_K,  Pressure(3.73810e-8), Density(3.88976e-14) },
    {  685.0_km,  999.94_K,  Pressure(3.58806e-8), Density(3.66028e-14) },
    {  690.0_km,  999.95_K,  Pressure(3.44737e-8), Density(3.44783e-14) },
    {  695.0_km,  999.95_K,  Pressure(3.31526e-8), Density(3.25089e-14) },
    {  700.0_km,  999.95_K,  Pressure(3.19103e-8), Density(3.06811e-14) },
    {  705.0_km,  999.95_K,  Pressure(3.07403e-8), Density(2.89823e-14) },
    {  710.0_km,  999.95_K,  Pressure(2.96365e-8), Density(2.74014e-14) },
    {  715.0_km,  999.96_K,  Pressure(2.85933e-8), Density(2.59281e-14) },
    {  720.0_km,  999.96_K,  Pressure(2.76057e-8), Density(2.45532e-14) },
    {  725.0_km,  999.96_K,  Pressure(2.66689e-8), Density(2.32685e-14) },
    {  730.0_km,  999.97_K,  Pressure(2.57788e-8), Density(2.20664e-14) },
    {  735.0_km,  999.97_K,  Pressure(2.49319e-8), Density(2.09405e-14) },
    {  740.0_km,  999.97_K,  Pressure(2.41248e-8), Density(1.98849e-14) },
    {  745.0_km,  999.97_K,  Pressure(2.33550e-8), Density(1.88944e-14) },
    {  750.0_km,  999.98_K,  Pressure(2.26202e-8), Density(1.79647e-14) },
    {  755.0_km,  999.98_K,  Pressure(2.19187e-8), Density(1.70917e-14) },
    {  760.0_km,  999.98_K,  Pressure(2.12493e-8), Density(1.62723e-14) },
    {  765.0_km,  999.98_K,  Pressure(2.06112e-8), Density(1.55035e-14) },
    {  770.0_km,  999.99_K,  Pressure(2.00041e-8), Density(1.47832e-14) },
    {  775.0_km,  999.99_K,  Pressure(1.94283e-8), Density(1.41093e-14) },
    {  780.0_km,  999.99_K,  Pressure(1.88844e-8), Density(1.34806e-14) },
    {  785.0_km,  999.99_K,  Pressure(1.83736e-8), Density(1.28959e-14) },
    {  790.0_km, 1000.00_K,  Pressure(1.78977e-8), Density(1.23545e-14) },
    {  795.0_km, 1000.00_K,  Pressure(1.74587e-8), Density(1.18561e-14) },
    {  800.0_km, 1000.00_K,  Pressure(1.70593e-8), Density(1.14006e-14) },
    {  805.0_km, 1000.00_K,  Pressure(1.66060e-8), Density(1.09248e-14) },
    {  810.0_km, 1000.00_K,  Pressure(1.61971e-8), Density(1.04935e-14) },
    {  815.0_km, 1000.00_K,  Pressure(1.58032e-8), Density(1.00859e-14) },
    {  820.0_km, 1000.00_K,  Pressure(1.54236e-8), Density(9.70075e-15) },
    {  825.0_km, 1000.00_K,  Pressure(1.50578e-8), Density(9.33686e-15) },
    {  830.0_km, 1000.00_K,  Pressure(1.47051e-8), Density(8.99305e-15) },
    {  835.0_km, 1000.00_K,  Pressure(1.43650e-8), Density(8.66819e-15) },
    {  840.0_km, 1000.00_K,  Pressure(1.40371e-8), Density(8.36124e-15) },
    {  845.0_km, 1000.00_K,  Pressure(1.37206e-8), Density(8.07122e-15) },
    {  850.0_km, 1000.00_K,  Pressure(1.34153e-8), Density(7.79720e-15) },
    {  855.0_km, 1000.00_K,  Pressure(1.31204e-8), Density(7.53828e-15) },
    {  860.0_km, 1000.00_K,  Pressure(1.28357e-8), Density(7.29363e-15) },
    {  865.0_km, 1000.00_K,  Pressure(1.25606e-8), Density(7.06247e-15) },
    {  870.0_km, 1000.00_K,  Pressure(1.22946e-8), Density(6.84404e-15) },
    {  875.0_km, 1000.00_K,  Pressure(1.20374e-8), Density(6.63764e-15) },
    {  880.0_km, 1000.00_K,  Pressure(1.17886e-8), Density(6.44261e-15) },
    {  885.0_km, 1000.00_K,  Pressure(1.15476e-8), Density(6.25832e-15) },
    {  890.0_km, 1000.00_K,  Pressure(1.13142e-8), Density(6.08416e-15) },
    {  895.0_km, 1000.00_K,  Pressure(1.10880e-8), Density(5.91959e-15) },
    {  900.0_km, 1000.00_K,  Pressure(1.08686e-8), Density(5.76406e-15) },
    {  905.0_km, 1000.00_K,  Pressure(1.06557e-8), Density(5.60198e-15) },
    {  910.0_km, 1000.00_K,  Pressure(1.04490e-8), Density(5.44733e-15) },
    {  915.0_km, 1000.00_K,  Pressure(1.02481e-8), Density(5.29993e-15) },
    {  920.0_km, 1000.00_K,  Pressure(1.00529e-8), Density(5.15935e-15) },
    {  925.0_km, 1000.00_K,  Pressure(9.86287e-9), Density(5.02514e-15) },
    {  930.0_km, 1000.00_K,  Pressure(9.67793e-9), Density(4.89694e-15) },
    {  935.0_km, 1000.00_K,  Pressure(9.49779e-9), Density(4.77435e-15) },
    {  940.0_km, 1000.00_K,  Pressure(9.32222e-9), Density(4.65704e-15) },
    {  945.0_km, 1000.00_K,  Pressure(9.15102e-9), Density(4.54468e-15) },
    {  950.0_km, 1000.00_K,  Pressure(8.98397e-9), Density(4.43696e-15) },
    {  955.0_km, 1000.00_K,  Pressure(8.82091e-9), Density(4.33361e-15) },
    {  960.0_km, 1000.00_K,  Pressure(8.66167e-9), Density(4.23434e-15) },
    {  965.0_km, 1000.00_K,  Pressure(8.50609e-9), Density(4.13892e-15) },
    {  970.0_km, 1000.00_K,  Pressure(8.35405e-9), Density(4.04712e-15) },
    {  975.0_km, 1000.00_K,  Pressure(8.20542e-9), Density(3.95872e-15) },
    {  980.0_km, 1000.00_K,  Pressure(8.06010e-9), Density(3.87353e-15) },
    {  985.0_km, 1000.00_K,  Pressure(7.91801e-9), Density(3.79136e-15) },
    {  990.0_km, 1000.00_K,  Pressure(7.77907e-9), Density(3.71205e-15) },
    {  995.0_km, 1000.00_K,  Pressure(7.64322e-9), Density(3.63544e-15) },
    { 1000.0_km, 1000.00_K,  Pressure(7.51043e-9), Density(3.56140e-15) },
    { 1005.0_km, 1000.00_K,  Pressure(7.41143e-9), Density(3.50435e-15) },
    { 1010.0_km, 1000.00_K,  Pressure(7.31128e-9), Density(3.44759e-15) },
    { 1015.0_km, 1000.00_K,  Pressure(7.21239e-9), Density(3.39217e-15) },
    { 1020.0_km, 1000.00_K,  Pressure(7.11474e-9), Density(3.33801e-15) },
    { 1025.0_km, 1000.00_K,  Pressure(7.01834e-9), Density(3.28502e-15) },
    { 1030.0_km, 1000.00_K,  Pressure(6.92319e-9), Density(3.23313e-15) },
    { 1035.0_km, 1000.00_K,  Pressure(6.82929e-9), Density(3.18226e-15) },
    { 1040.0_km, 1000.00_K,  Pressure(6.73664e-9), Density(3.13235e-15) },
    { 1045.0_km, 1000.00_K,  Pressure(6.64524e-9), Density(3.08332e-15) },
    { 1050.0_km, 1000.00_K,  Pressure(6.55509e-9), Density(3.03513e-15) },
    { 1055.0_km, 1000.00_K,  Pressure(6.46618e-9), Density(2.98845e-15) },
    { 1060.0_km, 1000.00_K,  Pressure(6.37853e-9), Density(2.94230e-15) },
    { 1065.0_km, 1000.00_K,  Pressure(6.29213e-9), Density(2.89697e-15) },
    { 1070.0_km, 1000.00_K,  Pressure(6.20697e-9), Density(2.85243e-15) },
    { 1075.0_km, 1000.00_K,  Pressure(6.12307e-9), Density(2.80870e-15) },
    { 1080.0_km, 1000.00_K,  Pressure(6.04041e-9), Density(2.76575e-15) },
    { 1085.0_km, 1000.00_K,  Pressure(5.95901e-9), Density(2.72358e-15) },
    { 1090.0_km, 1000.00_K,  Pressure(5.87885e-9), Density(2.68219e-15) },
    { 1095.0_km, 1000.00_K,  Pressure(5.79994e-9), Density(2.64156e-15) },
    { 1100.0_km, 1000.00_K,  Pressure(5.72228e-9), Density(2.60170e-15) },
    { 1105.0_km, 1000.00_K,  Pressure(5.64588e-9), Density(2.56259e-15) },
    { 1110.0_km, 1000.00_K,  Pressure(5.57072e-9), Density(2.52423e-15) },
    { 1115.0_km, 1000.00_K,  Pressure(5.49681e-9), Density(2.48661e-15) },
    { 1120.0_km, 1000.00_K,  Pressure(5.42415e-9), Density(2.44972e-15) },
    { 1125.0_km, 1000.00_K,  Pressure(5.35273e-9), Density(2.41357e-15) },
    { 1130.0_km, 1000.00_K,  Pressure(5.28257e-9), Density(2.37814e-15) },
    { 1135.0_km, 1000.00_K,  Pressure(5.21366e-9), Density(2.34342e-15) },
    { 1140.0_km, 1000.00_K,  Pressure(5.14600e-9), Density(2.30942e-15) },
    { 1145.0_km, 1000.00_K,  Pressure(5.07958e-9), Density(2.27613e-15) },
    { 1150.0_km, 1000.00_K,  Pressure(5.01442e-9), Density(2.24353e-15) },
    { 1155.0_km, 1000.00_K,  Pressure(4.95050e-9), Density(2.21163e-15) },
    { 1160.0_km, 1000.00_K,  Pressure(4.88784e-9), Density(2.18042e-15) },
    { 1165.0_km, 1000.00_K,  Pressure(4.82642e-9), Density(2.14990e-15) },
    { 1170.0_km, 1000.00_K,  Pressure(4.76626e-9), Density(2.12006e-15) },
    { 1175.0_km, 1000.00_K,  Pressure(4.70734e-9), Density(2.09089e-15) },
    { 1180.0_km, 1000.00_K,  Pressure(4.64967e-9), Density(2.06240e-15) },
    { 1185.0_km, 1000.00_K,  Pressure(4.59325e-9), Density(2.03457e-15) },
    { 1190.0_km, 1000.00_K,  Pressure(4.53808e-9), Density(2.00740e-15) },
    { 1195.0_km, 1000.00_K,  Pressure(4.48416e-9), Density(1.98089e-15) },
    { 1200.0_km, 1000.00_K,  Pressure(4.43149e-9), Density(1.95503e-15) }
  };

  //=========================================================================//
  // "InterpolateUpper": A Helper Function:                                  //
  //=========================================================================//
  inline AtmConds InterpolateUpper
  (
    LenK                 a_z,
    int                  DEBUG_ONLY(a_n),
    UpperLayerInfo const a_table[],  // Of length "a_n"
    LenK                 a_step
  )
  {
    assert(a_n >= 2  && a_table != nullptr && IsPos(a_step));
    assert(a_table[0].m_baseZ < a_z && a_z <= a_table[a_n-1].m_baseZ);

    // Bracket the Table entry:
    int i = int(Ceil(double((a_z - a_table[0].m_baseZ) / a_step)));
    assert(1 <= i && i <= a_n-1);

    UpperLayerInfo const& left  = a_table[i-1];
    UpperLayerInfo const& right = a_table[i];
    assert(left.m_baseZ < a_z && a_z <= right.m_baseZ);

    // Now perform linear interpolation:
    double   x = double((a_z - left.m_baseZ) / a_step);
    assert(0.0 < x && x <= 1.0);
    Pressure p   = (1.0-x) * left.m_baseP   + x * right.m_baseP;
    Density  rho = (1.0-x) * left.m_baseRho + x * right.m_baseRho;
    AbsTemp  T   = (1.0-x) * left.m_baseT   + x * right.m_baseT;
    Vel      a   = SqRt(GammaAir * p / rho);
    return std::make_tuple(p, rho, T, a);
  }

  //=========================================================================//
  // "GetAtmConds":                                                          //
  //=========================================================================//
  // NB: here "a_z" is a Geometric Altitude, NOT The GeoPotential one.
  // Returns (P, Rho, T, SpeedOfSound):
  //
  AtmConds GetAtmConds(LenK a_z)
  {
    // XXX: For the moment, altitudes below the MSL are not allowed, so adjust
    // "a_z" if necessary:
    a_z = std::max(a_z, 0.0_km);

    if (a_z <= UpperLayers1[0].m_baseZ)
    {
      // Compute the GeoPotential altitude from the Geometric one. For this con-
      // version, some "conventional" Earth radius is used, which is close to,
      // but not exactly equal to, the Earth polar radius:
      constexpr LenK Rc = 6356.767_km;
      LenK   h  = Rc * a_z / (Rc + a_z);

      // The "Low" model must be applicable:
      assert(h <= a_z && h < LowLayers[NLowLayers-1].m_endH);

      for (LowLayerInfo const& l:  LowLayers)
        if (l.m_baseH <= h  && h <= l.m_endH)
        {
          // Found the "Low" Layer which "h" belongs to:
          Pressure p   =  l.P(h);
          AbsTemp  T   =  l.T(h);
          assert(IsPos(p) &&  IsPos(T));
          Density  rho =  p / (RAir * T);
          Vel      a   =  SqRt(GammaAir * p / rho);
          return std::make_tuple(p, rho, T, a);
        }
      // This point must not be reached:
      __builtin_unreachable();
    }
    else
    // "UpperLayers*" are just interpolation tables wrt "a_z":
    if (a_z <= UpperLayers1[NUpperLayers1-1].m_baseZ)
      return InterpolateUpper(a_z, NUpperLayers1, UpperLayers1, 1.0_km);
    else
    if (a_z <= UpperLayers2[NUpperLayers2-1].m_baseZ)
      return InterpolateUpper(a_z, NUpperLayers2, UpperLayers2, 2.0_km);
    else
    if (a_z <= UpperLayers5[NUpperLayers5-1].m_baseZ)
      return InterpolateUpper(a_z, NUpperLayers5, UpperLayers5, 5.0_km);

    // If we got here: We must be above all Layers, so presumably in the abso-
    // lute vacuum. So the Pressure and Density are 0; but the AbsTemperature
    // and the Speed of Sound are non-0:
    return
      std::make_tuple
      (
        Pressure(0.0),
        Density (0.0),
        UpperLayers5[NUpperLayers5-1].m_baseT,
        SqRt(GammaAir * UpperLayers5[NUpperLayers5-1].m_baseP  /
                        UpperLayers5[NUpperLayers5-1].m_baseRho)
      );
  }
}
// End namespace SpaceBallistics::EarthAtmosphereModel

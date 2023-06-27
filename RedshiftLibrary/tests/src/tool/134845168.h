// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/common/datatypes.h"

using namespace NSEpic;

static TFloat64List lambdaE = {
    12504.0, 12517.4, 12530.8, 12544.2, 12557.6, 12571.0, 12584.4, 12597.8,
    12611.2, 12624.6, 12638.0, 12651.4, 12664.8, 12678.2, 12691.6, 12705.0,
    12718.4, 12731.8, 12745.2, 12758.6, 12772.0, 12785.4, 12798.8, 12812.2,
    12825.6, 12839.0, 12852.4, 12865.8, 12879.2, 12892.6, 12906.0, 12919.4,
    12932.8, 12946.2, 12959.6, 12973.0, 12986.4, 12999.8, 13013.2, 13026.6,
    13040.0, 13053.4, 13066.8, 13080.2, 13093.6, 13107.0, 13120.4, 13133.8,
    13147.2, 13160.6, 13174.0, 13187.4, 13200.8, 13214.2, 13227.6, 13241.0,
    13254.4, 13267.8, 13281.2, 13294.6, 13308.0, 13321.4, 13334.8, 13348.2,
    13361.6, 13375.0, 13388.4, 13401.8, 13415.2, 13428.6, 13442.0, 13455.4,
    13468.8, 13482.2, 13495.6, 13509.0, 13522.4, 13535.8, 13549.2, 13562.6,
    13576.0, 13589.4, 13602.8, 13616.2, 13629.6, 13643.0, 13656.4, 13669.8,
    13683.2, 13696.6, 13710.0, 13723.4, 13736.8, 13750.2, 13763.6, 13777.0,
    13790.4, 13803.8, 13817.2, 13830.6, 13844.0, 13857.4, 13870.8, 13884.2,
    13897.6, 13911.0, 13924.4, 13937.8, 13951.2, 13964.6, 13978.0, 13991.4,
    14004.8, 14018.2, 14031.6, 14045.0, 14058.4, 14071.8, 14085.2, 14098.6,
    14112.0, 14125.4, 14138.8, 14152.2, 14165.6, 14179.0, 14192.4, 14205.8,
    14219.2, 14232.6, 14246.0, 14259.4, 14272.8, 14286.2, 14299.6, 14313.0,
    14326.4, 14339.8, 14353.2, 14366.6, 14380.0, 14393.4, 14406.8, 14420.2,
    14433.6, 14447.0, 14460.4, 14473.8, 14487.2, 14500.6, 14514.0, 14527.4,
    14540.8, 14554.2, 14567.6, 14581.0, 14594.4, 14607.8, 14621.2, 14634.6,
    14648.0, 14661.4, 14674.8, 14688.2, 14701.6, 14715.0, 14728.4, 14741.8,
    14755.2, 14768.6, 14782.0, 14795.4, 14808.8, 14822.2, 14835.6, 14849.0,
    14862.4, 14875.8, 14889.2, 14902.6, 14916.0, 14929.4, 14942.8, 14956.2,
    14969.6, 14983.0, 14996.4, 15009.8, 15023.2, 15036.6, 15050.0, 15063.4,
    15076.8, 15090.2, 15103.6, 15117.0, 15130.4, 15143.8, 15157.2, 15170.6,
    15184.0, 15197.4, 15210.8, 15224.2, 15237.6, 15251.0, 15264.4, 15277.8,
    15291.2, 15304.6, 15318.0, 15331.4, 15344.8, 15358.2, 15371.6, 15385.0,
    15398.4, 15411.8, 15425.2, 15438.6, 15452.0, 15465.4, 15478.8, 15492.2,
    15505.6, 15519.0, 15532.4, 15545.8, 15559.2, 15572.6, 15586.0, 15599.4,
    15612.8, 15626.2, 15639.6, 15653.0, 15666.4, 15679.8, 15693.2, 15706.6,
    15720.0, 15733.4, 15746.8, 15760.2, 15773.6, 15787.0, 15800.4, 15813.8,
    15827.2, 15840.6, 15854.0, 15867.4, 15880.8, 15894.2, 15907.6, 15921.0,
    15934.4, 15947.8, 15961.2, 15974.6, 15988.0, 16001.4, 16014.8, 16028.2,
    16041.6, 16055.0, 16068.4, 16081.8, 16095.2, 16108.6, 16122.0, 16135.4,
    16148.8, 16162.2, 16175.6, 16189.0, 16202.4, 16215.8, 16229.2, 16242.6,
    16256.0, 16269.4, 16282.8, 16296.2, 16309.6, 16323.0, 16336.4, 16349.8,
    16363.2, 16376.6, 16390.0, 16403.4, 16416.8, 16430.2, 16443.6, 16457.0,
    16470.4, 16483.8, 16497.2, 16510.6, 16524.0, 16537.4, 16550.8, 16564.2,
    16577.6, 16591.0, 16604.4, 16617.8, 16631.2, 16644.6, 16658.0, 16671.4,
    16684.8, 16698.2, 16711.6, 16725.0, 16738.4, 16751.8, 16765.2, 16778.6,
    16792.0, 16805.4, 16818.8, 16832.2, 16845.6, 16859.0, 16872.4, 16885.8,
    16899.2, 16912.6, 16926.0, 16939.4, 16952.8, 16966.2, 16979.6, 16993.0,
    17006.4, 17019.8, 17033.2, 17046.6, 17060.0, 17073.4, 17086.8, 17100.2,
    17113.6, 17127.0, 17140.4, 17153.8, 17167.2, 17180.6, 17194.0, 17207.4,
    17220.8, 17234.2, 17247.6, 17261.0, 17274.4, 17287.8, 17301.2, 17314.6,
    17328.0, 17341.4, 17354.8, 17368.2, 17381.6, 17395.0, 17408.4, 17421.8,
    17435.2, 17448.6, 17462.0, 17475.4, 17488.8, 17502.2, 17515.6, 17529.0,
    17542.4, 17555.8, 17569.2, 17582.6, 17596.0, 17609.4, 17622.8, 17636.2,
    17649.6, 17663.0, 17676.4, 17689.8, 17703.2, 17716.6, 17730.0, 17743.4,
    17756.8, 17770.2, 17783.6, 17797.0, 17810.4, 17823.8, 17837.2, 17850.6,
    17864.0, 17877.4, 17890.8, 17904.2, 17917.6, 17931.0, 17944.4, 17957.8,
    17971.2, 17984.6, 17998.0, 18011.4, 18024.8, 18038.2, 18051.6, 18065.0,
    18078.4, 18091.8, 18105.2, 18118.6, 18132.0, 18145.4, 18158.8, 18172.2,
    18185.6, 18199.0, 18212.4, 18225.8, 18239.2, 18252.6, 18266.0, 18279.4,
    18292.8, 18306.2, 18319.6, 18333.0, 18346.4, 18359.8, 18373.2, 18386.6,
    18400.0, 18413.4, 18426.8, 18440.2, 18453.6, 18467.0, 18480.4, 18493.8,
    18507.2};

static TFloat64List fluxE = {
    -3.3164543e-18, -9.50855e-18,   7.927449e-18,   -4.860559e-18,
    2.704523e-18,   2.4812736e-18,  -6.6802254e-18, 8.78693e-18,
    4.890033e-18,   -5.238228e-18,  3.368171e-18,   6.780399e-19,
    1.0978484e-18,  -1.1254897e-18, 2.3539774e-18,  5.526356e-18,
    5.775241e-18,   -3.5369131e-19, -2.86579e-18,   3.6498612e-18,
    8.770006e-18,   7.489351e-18,   8.00776e-18,    -4.7717293e-20,
    -4.1422227e-18, 3.3946342e-19,  -3.161433e-18,  6.128694e-18,
    -1.4106799e-18, 2.6908254e-18,  -2.5221564e-18, -4.6660522e-18,
    -2.497229e-19,  -1.4139764e-18, 3.6294952e-18,  -5.2861086e-19,
    3.5099755e-18,  -1.3764475e-18, 1.3228349e-18,  -5.3345266e-18,
    3.6892657e-18,  4.0006934e-18,  6.001309e-18,   1.0851079e-18,
    -1.741175e-18,  -1.6461417e-18, 1.4014372e-18,  2.2676003e-18,
    7.6008696e-19,  -7.71918e-18,   -6.757497e-20,  5.608371e-18,
    -2.3397503e-18, 3.0741551e-18,  6.387304e-18,   -4.913121e-18,
    1.3587652e-18,  5.2582296e-18,  4.6901927e-19,  -1.9964327e-18,
    9.309901e-19,   6.412932e-18,   2.2716082e-18,  -8.178568e-18,
    -6.920274e-19,  3.3853266e-18,  2.589927e-18,   -2.189216e-18,
    8.28628e-18,    5.351779e-18,   1.8667572e-18,  1.0848347e-18,
    -2.1241523e-18, -2.3208494e-18, 1.2676065e-18,  3.021456e-18,
    5.871011e-18,   8.321748e-18,   2.4760394e-18,  -4.9666463e-18,
    4.879734e-19,   -1.2305288e-20, 2.959369e-18,   1.475543e-18,
    1.4578141e-18,  9.541343e-19,   -5.069152e-19,  5.729362e-18,
    -4.2797175e-18, 3.7727807e-18,  -3.947708e-18,  -1.3719433e-18,
    5.058112e-19,   7.973074e-18,   2.3211217e-18,  -1.7450337e-18,
    1.678685e-18,   7.8070075e-18,  -6.227218e-19,  2.3835414e-18,
    3.47678e-18,    3.2337927e-18,  1.7620992e-18,  -9.156165e-19,
    -3.4938178e-18, 1.0694618e-18,  -1.727041e-18,  3.1914946e-18,
    2.6225559e-18,  -6.1063556e-18, -5.917207e-19,  -3.8382384e-18,
    2.7213589e-18,  4.1900424e-18,  6.5068395e-19,  -6.1310353e-19,
    3.94113e-18,    4.6265e-19,     4.508121e-19,   -1.5575869e-18,
    3.8102987e-18,  -1.7064555e-19, 1.4109981e-18,  7.1466265e-19,
    -3.8990788e-18, 2.6501355e-18,  2.884465e-18,   1.1964168e-17,
    -4.2375077e-18, 1.2862503e-18,  -3.4437398e-18, -2.0661454e-18,
    2.286797e-18,   -1.622151e-18,  -1.0588263e-18, 3.3800912e-18,
    2.7503214e-18,  -5.3780476e-19, -2.7366414e-18, 2.2218516e-19,
    3.640609e-18,   -3.446299e-18,  5.2128837e-21,  3.675561e-18,
    -2.7871921e-18, -3.347529e-18,  2.2079525e-18,  -9.405338e-19,
    -3.8788187e-19, 2.9747063e-18,  -3.89629e-19,   3.7849556e-18,
    1.3355489e-19,  -1.1904998e-18, 1.5829013e-18,  -2.443634e-18,
    2.6158063e-18,  -8.139586e-19,  -6.7449403e-19, 2.8588378e-19,
    -4.28631e-19,   5.067668e-18,   1.8787286e-18,  -3.4502826e-18,
    -5.843864e-18,  1.9200297e-18,  -1.7739562e-18, 9.995939e-19,
    2.0742828e-18,  5.2930427e-18,  4.387525e-18,   -2.6448056e-18,
    9.501946e-19,   5.062836e-18,   6.208692e-19,   6.2755996e-19,
    -3.9125556e-18, -2.7043926e-18, -3.9334523e-18, 2.2466073e-19,
    1.1219629e-18,  2.1695875e-19,  4.163647e-18,   6.0336597e-18,
    -1.2100177e-18, 2.9232778e-18,  8.294803e-19,   8.3770685e-19,
    1.8918266e-18,  3.0632229e-18,  3.886413e-18,   -3.4311717e-18,
    -1.0429714e-18, 8.818304e-19,   -2.462971e-18,  -2.2945651e-18,
    -6.653642e-19,  8.009738e-19,   -1.2350815e-18, 1.5106396e-18,
    9.784493e-19,   1.7528023e-18,  9.645465e-19,   -7.423104e-18,
    2.1608071e-18,  3.153454e-18,   1.0042749e-17,  2.2566384e-19,
    -4.2868087e-19, -5.1336286e-20, 3.1018121e-18,  5.4123544e-19,
    5.475417e-18,   5.78558e-19,    -1.1231198e-18, 2.5022207e-18,
    1.3290762e-18,  -2.804878e-18,  -1.6054267e-18, 3.75982e-18,
    -1.3547496e-18, 7.879798e-19,   1.7417103e-18,  2.5956615e-18,
    3.499134e-18,   4.1457035e-18,  3.771195e-18,   5.416691e-18,
    3.5995083e-18,  7.822691e-18,   5.100121e-18,   9.381267e-18,
    7.399297e-18,   8.29245e-18,    1.212533e-17,   3.5404873e-18,
    7.7647775e-18,  7.837518e-18,   2.1077672e-18,  4.419915e-18,
    3.2983386e-18,  5.6615674e-18,  -1.4285568e-18, 3.193497e-18,
    5.2401962e-18,  -5.1765256e-18, -1.1126763e-18, -2.6991303e-18,
    -2.24636e-18,   3.513318e-18,   4.166452e-18,   2.4291722e-18,
    8.92015e-19,    4.1616226e-18,  1.1409802e-18,  1.0213494e-18,
    3.0248192e-18,  4.2195724e-18,  2.0753445e-18,  5.9632265e-18,
    5.3178552e-18,  5.505043e-19,   4.898069e-18,   4.9055748e-18,
    3.6531183e-18,  -8.632689e-19,  8.52978e-19,    7.430364e-18,
    -6.1565564e-19, -1.2701776e-18, 2.483039e-18,   4.5497345e-18,
    -9.201721e-18,  5.423413e-19,   -9.775351e-19,  -7.1828386e-18,
    -1.940068e-19,  -3.536816e-18,  -2.1114178e-18, -7.351181e-19,
    3.1190507e-19,  5.9778788e-18,  6.275202e-19,   4.4203007e-20,
    -1.6659283e-18, 1.9130103e-18,  -2.6901522e-18, -2.8196114e-18,
    4.0583574e-18,  -4.3536837e-19, 1.6460635e-18,  1.7132297e-18,
    -2.217719e-18,  6.1968037e-18,  -2.2644998e-18, -2.2888743e-18,
    -3.1413398e-18, 5.2970697e-19,  -2.8770055e-18, 2.3033514e-18,
    4.5146773e-18,  -2.445541e-18,  -4.297836e-18,  4.097878e-18,
    -1.7722948e-19, 2.1015214e-18,  1.134999e-18,   5.915723e-19,
    -1.0125114e-18, 1.1517001e-18,  -1.8379068e-18, -6.7304386e-18,
    5.541154e-18,   -2.0314722e-18, 7.547345e-19,   3.4768837e-18,
    1.9356396e-18,  2.7975935e-18,  -1.7723443e-18, -3.4457444e-19,
    3.8949693e-18,  9.528665e-19,   3.5797167e-18,  -1.0536655e-18,
    6.564251e-19,   2.2845866e-18,  1.2141395e-19,  -3.2748566e-18,
    1.7824545e-18,  2.0049704e-18,  1.8551192e-18,  -3.043668e-18,
    5.1096523e-19,  -5.4512795e-19, 4.7175177e-18,  1.5131799e-18,
    4.751609e-18,   9.24084e-20,    -4.0888336e-19, 2.8557176e-18,
    3.2878357e-18,  1.1475713e-18,  -3.5430412e-18, 1.7183813e-18,
    -2.2020727e-18, -5.926522e-19,  -1.1544795e-18, -4.0477905e-19,
    2.6637292e-18,  -3.1396376e-18, -1.1858407e-18, 8.412668e-19,
    1.7625742e-18,  2.0790734e-18,  -5.982997e-19,  -4.677867e-18,
    7.812794e-19,   1.926261e-18,   2.272424e-18,   -2.6151414e-18,
    3.096784e-18,   -1.6863236e-18, -8.787915e-19,  1.9505096e-19,
    4.1804612e-19,  5.7014607e-18,  2.735607e-18,   2.4805565e-18,
    2.2673868e-18,  2.0434304e-18,  -9.645557e-19,  6.1415906e-18,
    7.882957e-19,   6.4269354e-19,  3.1018677e-18,  -4.391776e-19,
    7.0657205e-18,  -7.262734e-19,  5.4976037e-19,  -4.5177777e-19,
    1.5320625e-18,  -1.7140729e-18, -5.2844694e-18, 2.0397009e-19,
    2.0208647e-19,  4.4456466e-18,  -4.0758885e-19, 2.464714e-18,
    2.8219974e-18,  1.3157667e-19,  2.2650236e-18,  4.631635e-19,
    2.02962e-18,    1.8747472e-18,  3.899995e-18,   -4.2000318e-19,
    -3.0990084e-18, -1.9120497e-18, -5.0613843e-19, 1.6727952e-18,
    4.2872156e-19,  2.1976284e-18,  -7.876891e-20,  2.9592107e-18,
    -5.7002906e-18, 4.9714155e-18,  -1.058297e-18,  -3.8893974e-20,
    2.1768033e-18,  1.4380583e-18,  -2.2511824e-18, 5.388172e-19,
    1.980079e-18,   3.595305e-18,   6.05858e-19,    -4.2120567e-19,
    2.3841204e-18,  -2.7368043e-19, -2.9152466e-19, 3.270174e-18,
    1.2568196e-18,  3.1209504e-18,  8.0161556e-19,  -2.0960002e-18,
    -1.6493086e-18, 9.173737e-19,   5.1599e-21,     1.3769807e-19,
    -3.4951657e-18, 4.5733215e-18,  -1.1054583e-18, 7.4400355e-18,
    -2.2464606e-18, -8.182334e-19,  5.1766046e-18,  3.1489463e-18,
    3.114878e-19,   -8.054669e-20,  4.8045545e-18,  -6.1064365e-19,
    -2.4850865e-18, -6.544497e-19,  2.2808341e-18,  -2.1211475e-18,
    1.8717331e-18,  1.2674749e-18,  1.9711932e-19,  -1.4214297e-18,
    -1.5501066e-18};

static TFloat64List errorE = {
    4.4765956e-18, 4.446956e-18,  4.421172e-18,  4.4119812e-18, 4.3929734e-18,
    4.3596087e-18, 4.3327563e-18, 4.3294617e-18, 4.3026953e-18, 4.3120545e-18,
    4.2746858e-18, 4.2699824e-18, 4.2455595e-18, 4.2206303e-18, 4.2370486e-18,
    4.203658e-18,  4.1841488e-18, 4.185565e-18,  4.1674468e-18, 4.1647717e-18,
    4.158309e-18,  4.128793e-18,  4.108785e-18,  4.1039495e-18, 4.0829565e-18,
    4.0776547e-18, 4.0687286e-18, 4.0684386e-18, 4.045259e-18,  4.032103e-18,
    4.0258963e-18, 4.0224552e-18, 3.9864807e-18, 3.9697812e-18, 3.9743443e-18,
    3.9601867e-18, 3.949895e-18,  3.9621075e-18, 3.9452243e-18, 3.9259067e-18,
    3.9152096e-18, 3.9338497e-18, 3.918427e-18,  3.9054245e-18, 3.895635e-18,
    3.8902164e-18, 3.892634e-18,  3.898029e-18,  3.900532e-18,  3.8874755e-18,
    3.8957403e-18, 3.884371e-18,  3.8807526e-18, 3.8679147e-18, 3.863763e-18,
    3.872015e-18,  3.8537973e-18, 3.861002e-18,  3.850105e-18,  3.848995e-18,
    3.8344987e-18, 3.8232623e-18, 3.836815e-18,  3.8255155e-18, 3.8283242e-18,
    3.8178773e-18, 3.805055e-18,  3.813186e-18,  3.810481e-18,  3.799486e-18,
    3.790803e-18,  3.8018185e-18, 3.781881e-18,  3.768365e-18,  3.773004e-18,
    3.7648944e-18, 3.7672283e-18, 3.7612647e-18, 3.7575155e-18, 3.7513422e-18,
    3.743616e-18,  3.7524982e-18, 3.7464706e-18, 3.7442343e-18, 3.7294964e-18,
    3.743239e-18,  3.743036e-18,  3.7247294e-18, 3.7247873e-18, 3.7257187e-18,
    3.7251785e-18, 3.7124602e-18, 3.723054e-18,  3.7322373e-18, 3.706362e-18,
    3.70792e-18,   3.708176e-18,  3.707497e-18,  3.7039465e-18, 3.6894257e-18,
    3.7049647e-18, 3.685445e-18,  3.6747755e-18, 3.6782497e-18, 3.6734388e-18,
    3.6712906e-18, 3.6672176e-18, 3.666988e-18,  3.663834e-18,  3.663736e-18,
    3.6490887e-18, 3.6673863e-18, 3.6303725e-18, 3.633378e-18,  3.622372e-18,
    3.632355e-18,  3.6189747e-18, 3.624297e-18,  3.626855e-18,  3.6017533e-18,
    3.6145096e-18, 3.597016e-18,  3.057625e-18,  3.0487623e-18, 3.0443561e-18,
    3.049841e-18,  3.0432124e-18, 3.0328616e-18, 3.030999e-18,  3.032754e-18,
    3.0267145e-18, 3.013869e-18,  3.0177156e-18, 3.0093087e-18, 3.0081509e-18,
    3.0072556e-18, 3.0025713e-18, 2.9935735e-18, 2.9903551e-18, 2.9844559e-18,
    2.981177e-18,  2.975029e-18,  2.9662717e-18, 2.9635526e-18, 2.9638003e-18,
    2.959262e-18,  2.9623042e-18, 2.9552121e-18, 2.955983e-18,  2.965272e-18,
    2.9492788e-18, 2.9463266e-18, 2.9513192e-18, 2.957085e-18,  2.942376e-18,
    2.9394159e-18, 2.934404e-18,  2.9452024e-18, 2.9192165e-18, 2.9348023e-18,
    2.9230511e-18, 2.9115513e-18, 2.9104143e-18, 2.9173291e-18, 2.9186094e-18,
    2.911527e-18,  2.9104195e-18, 2.9027317e-18, 2.8994246e-18, 2.905929e-18,
    2.894851e-18,  2.900575e-18,  2.8857794e-18, 2.8788607e-18, 2.8902888e-18,
    2.877218e-18,  2.8825803e-18, 2.8788489e-18, 2.8740266e-18, 2.8670196e-18,
    2.8704267e-18, 2.86374e-18,   2.8654252e-18, 2.8596287e-18, 2.8505548e-18,
    2.8470806e-18, 2.8428947e-18, 2.8335537e-18, 2.8379103e-18, 2.8418015e-18,
    2.8266761e-18, 2.8283458e-18, 2.825382e-18,  2.817861e-18,  2.8181088e-18,
    2.815795e-18,  2.8059085e-18, 2.8058167e-18, 2.80164e-18,   2.8032793e-18,
    2.7970868e-18, 2.7952296e-18, 2.7851827e-18, 2.7891403e-18, 2.7852786e-18,
    2.784378e-18,  2.7805498e-18, 2.7789342e-18, 2.776761e-18,  2.7615433e-18,
    2.7606696e-18, 2.757504e-18,  2.7582216e-18, 2.760458e-18,  2.75805e-18,
    2.7496116e-18, 2.7484335e-18, 2.7499253e-18, 2.7543894e-18, 2.7562303e-18,
    2.7425794e-18, 2.7380505e-18, 2.7409322e-18, 2.7551906e-18, 2.7515568e-18,
    2.7440247e-18, 2.7466342e-18, 2.7457098e-18, 2.7465298e-18, 2.7486476e-18,
    2.7461745e-18, 2.7418481e-18, 2.749795e-18,  2.7445782e-18, 2.7399382e-18,
    2.7320465e-18, 2.730371e-18,  2.7197233e-18, 2.7197413e-18, 2.7154623e-18,
    2.7144385e-18, 2.7015067e-18, 2.6986193e-18, 2.6993282e-18, 2.6961644e-18,
    2.6904875e-18, 2.682976e-18,  2.6776184e-18, 2.6709281e-18, 2.6686654e-18,
    2.6675268e-18, 2.6621995e-18, 2.6611203e-18, 2.6614532e-18, 2.6644598e-18,
    2.6552618e-18, 2.6504101e-18, 2.6478732e-18, 2.6560629e-18, 2.6399105e-18,
    2.6487545e-18, 2.6443105e-18, 2.6520792e-18, 2.6418478e-18, 2.6402927e-18,
    2.635714e-18,  2.632607e-18,  2.6373405e-18, 2.6364711e-18, 2.6293055e-18,
    2.6272762e-18, 2.6240342e-18, 2.618436e-18,  2.6260844e-18, 2.6203241e-18,
    2.618733e-18,  2.6196434e-18, 2.6164988e-18, 2.6221843e-18, 2.614386e-18,
    2.625017e-18,  2.6077113e-18, 2.6091582e-18, 2.618075e-18,  2.6134738e-18,
    2.6107063e-18, 2.6074553e-18, 2.613787e-18,  2.6039935e-18, 2.5995859e-18,
    2.6061649e-18, 2.602984e-18,  2.6100675e-18, 2.6031715e-18, 2.5958179e-18,
    2.5980018e-18, 2.5922275e-18, 2.6005632e-18, 2.5884988e-18, 2.5915728e-18,
    2.5862792e-18, 2.5887794e-18, 2.5888613e-18, 2.5786698e-18, 2.5829364e-18,
    2.583202e-18,  2.5745664e-18, 2.5729749e-18, 2.5743774e-18, 2.5699941e-18,
    2.5778195e-18, 2.570773e-18,  2.5629263e-18, 2.5569799e-18, 2.572164e-18,
    2.5653567e-18, 2.568168e-18,  2.5641714e-18, 2.570463e-18,  2.559795e-18,
    2.5631273e-18, 2.5631986e-18, 2.5706552e-18, 2.5639441e-18, 2.564145e-18,
    2.5637365e-18, 2.5687217e-18, 2.557863e-18,  2.5613548e-18, 2.5641834e-18,
    2.5633242e-18, 2.5573856e-18, 2.5649709e-18, 2.562483e-18,  2.561887e-18,
    2.5589074e-18, 2.56801e-18,   2.5679142e-18, 2.5701482e-18, 2.569822e-18,
    2.5623274e-18, 2.5591188e-18, 2.5616092e-18, 2.5639741e-18, 2.56683e-18,
    2.5606648e-18, 2.5603537e-18, 2.5655238e-18, 2.5587062e-18, 2.5440399e-18,
    2.5517655e-18, 2.5407777e-18, 2.5313085e-18, 2.5240114e-18, 2.5144061e-18,
    2.5077556e-18, 2.498414e-18,  2.4880085e-18, 2.4803989e-18, 2.4776001e-18,
    2.466405e-18,  2.4489915e-18, 2.4484616e-18, 2.4457253e-18, 2.4380644e-18,
    2.436119e-18,  2.4315994e-18, 2.4297177e-18, 2.429212e-18,  2.4280913e-18,
    2.4275555e-18, 2.4277046e-18, 2.4292045e-18, 2.4302924e-18, 2.4203017e-18,
    2.4356623e-18, 2.4244906e-18, 2.4299187e-18, 2.42623e-18,   2.4267837e-18,
    2.4227611e-18, 2.4261735e-18, 2.4323585e-18, 2.4264419e-18, 2.4315373e-18,
    2.4242418e-18, 2.4281126e-18, 2.4203607e-18, 2.4268052e-18, 2.426715e-18,
    2.4302033e-18, 2.4288746e-18, 2.424281e-18,  2.427096e-18,  2.4301274e-18,
    2.4220436e-18, 2.4268977e-18, 2.4251217e-18, 2.432112e-18,  2.4293616e-18,
    2.4222313e-18, 2.4269829e-18, 2.427792e-18,  2.4303977e-18, 2.427647e-18,
    2.4240993e-18, 2.4298755e-18, 2.4388767e-18, 2.4328447e-18, 2.4319317e-18,
    2.422645e-18,  2.4314132e-18, 2.4358366e-18, 2.4340383e-18, 2.4320481e-18,
    2.4343648e-18, 2.4396284e-18, 2.4346186e-18, 2.4315187e-18, 2.4345081e-18,
    2.43178e-18,   2.4337885e-18, 2.4392076e-18, 2.4299115e-18, 2.4280865e-18,
    2.4342552e-18, 2.4389596e-18, 2.4355626e-18, 2.4335676e-18, 2.4369802e-18,
    2.4304345e-18, 2.4377002e-18, 2.4357762e-18, 2.434902e-18,  2.4483719e-18,
    2.4431203e-18, 2.4436222e-18, 2.4384892e-18, 2.4393403e-18, 2.4517944e-18,
    2.4474413e-18, 2.4614291e-18, 2.4674518e-18, 2.4793783e-18, 2.4833604e-18,
    2.495418e-18,  2.5062605e-18, 2.5250203e-18, 2.5474776e-18};

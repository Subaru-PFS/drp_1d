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
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/exception.h"
#include <cstdio>
using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Convolve_test)
TFloat64Range _lbdaRange = {200., 1299.};
Float64 lbdaStep = 1;
TFloat64List fluxcorr = { // 2.0_1
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.27387,  0.27387,  0.27387,
    0.27387,  0.280386, 0.280386, 0.280386, 0.280386, 0.268539, 0.268539,
    0.268539, 0.287627, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.262616, 0.282477, 0.282477, 0.282477, 0.269236, 0.269236, 0.269236,
    0.269236, 0.290753, 0.290753, 0.276684, 0.276684, 0.276684, 0.276684,
    0.276684, 0.262616, 0.262616, 0.285125, 0.285125, 0.285125, 0.285125,
    0.270119, 0.270119, 0.270119, 0.270119, 0.270119, 0.294772, 0.278694,
    0.278694, 0.278694, 0.278694, 0.278694, 0.278694, 0.278694, 0.262616,
    0.262616, 0.262616, 0.288588, 0.288588, 0.288588, 0.288588, 0.271273,
    0.271273, 0.271273, 0.271273, 0.271273, 0.271273, 0.271273, 0.271273,
    0.300132, 0.281374, 0.281374, 0.281374, 0.281374, 0.281374, 0.281374,
    0.281374, 0.281374, 0.281374, 0.262616, 0.262616, 0.262616, 0.262616,
    0.262616, 0.293311, 0.293311, 0.293311, 0.293311, 0.293311, 0.272847,
    0.272847, 0.272847, 0.272847, 0.272847, 0.272847, 0.272847, 0.272847,
    0.272847, 0.272847, 0.272847, 0.272847, 0.272847, 0.252384, 0.252384,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125, 0.285125,
    0.262616, 0.262616, 0.262616, 0.262616, 0.262616, 0.262616, 0.262616,
    0.262616, 0.262616, 0.262616, 0.262616, 0.262616, 0.262616, 0.300132,
    0.300132, 0.300132, 0.300132, 0.300132, 0.300132, 0.300132, 0.300132,
    0.300132, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121, 0.275121,
    0.275121, 0.275185, 0.275257, 0.275334, 0.27542,  0.275512, 0.27561,
    0.275715, 0.275829, 0.275948, 0.276074, 0.276208, 0.276348, 0.276496,
    0.27665,  0.276813, 0.276982, 0.277158, 0.277341, 0.277532, 0.27773,
    0.277935, 0.278148, 0.278368, 0.278595, 0.27883,  0.279072, 0.279321,
    0.279578, 0.279844, 0.280115, 0.280396, 0.280683, 0.306522, 0.306852,
    0.30719,  0.307538, 0.307894, 0.308258, 0.308632, 0.309014, 0.27073,
    0.271079, 0.271437, 0.271801, 0.272174, 0.272556, 0.272944, 0.273341,
    0.273745, 0.274159, 0.27458,  0.275009, 0.275446, 0.299539, 0.300031,
    0.300534, 0.301046, 0.301564, 0.302094, 0.302632, 0.303181, 0.303738,
    0.304306, 0.304882, 0.305469, 0.306063, 0.271454, 0.294053, 0.294651,
    0.29526,  0.295878, 0.296504, 0.297142, 0.297789, 0.298445, 0.299113,
    0.299789, 0.300476, 0.323763, 0.324523, 0.325293, 0.326077, 0.32687,
    0.293386, 0.294117, 0.294859, 0.295612, 0.296375, 0.318373, 0.319213,
    0.320066, 0.320929, 0.321804, 0.322693, 0.323591, 0.324503, 0.313745,
    0.314647, 0.315562, 0.316489, 0.317427, 0.318378, 0.319342, 0.340762,
    0.341813, 0.342879, 0.343958, 0.313995, 0.315002, 0.33537,  0.336466,
    0.337575, 0.338698, 0.339837, 0.340988, 0.361894, 0.33277,  0.333929,
    0.335101, 0.336289, 0.356241, 0.357523, 0.358823, 0.360139, 0.332934,
    0.352077, 0.353401, 0.354741, 0.356096, 0.357468, 0.377105, 0.350565,
    0.351948, 0.35335,  0.372215, 0.37372,  0.375243, 0.350291, 0.368491,
    0.370027, 0.371582, 0.390116, 0.365781, 0.367354, 0.368945, 0.386901,
    0.388602, 0.390323, 0.382951, 0.384672, 0.386411, 0.404124, 0.381559,
    0.383321, 0.40051,  0.402387, 0.419837, 0.397685, 0.399589, 0.416569,
    0.418592, 0.397499, 0.414068, 0.416117, 0.410196, 0.412253, 0.428623,
    0.430802, 0.424922, 0.42711,  0.429324, 0.423987, 0.426212, 0.442139,
    0.423614, 0.43919,  0.454935, 0.43674,  0.452164, 0.454643, 0.449879,
    0.452377, 0.467779, 0.463072, 0.465691, 0.461401, 0.464042, 0.479104,
    0.474871, 0.477638, 0.473776, 0.488483, 0.484684, 0.487575, 0.484109,
    0.498537, 0.495132, 0.498155, 0.495056, 0.509264, 0.506227, 0.520381,
    0.517409, 0.531524, 0.528622, 0.515538, 0.529365, 0.526987, 0.524844,
    0.538495, 0.536415, 0.544452, 0.558029, 0.556081, 0.554345, 0.567815,
    0.566151, 0.564683, 0.572742, 0.586088, 0.58476,  0.583615, 0.591702,
    0.590793, 0.598902, 0.598219, 0.597694, 0.605892, 0.605576, 0.613813,
    0.622031, 0.621983, 0.630254, 0.638522, 0.638737, 0.647069, 0.655407,
    0.663757, 0.659948, 0.668362, 0.676797, 0.685254, 0.682016, 0.690554,
    0.699125, 0.696371, 0.70503,  0.709726, 0.71845,  0.723273, 0.732073,
    0.737029, 0.735401, 0.751014, 0.756212, 0.755027, 0.760513, 0.766099,
    0.771782, 0.777564, 0.783442, 0.795477, 0.801494, 0.807612, 0.819715,
    0.825986, 0.829128, 0.841328, 0.839144, 0.851451, 0.855211, 0.873032,
    0.876978, 0.881127, 0.885475, 0.895168, 0.904927, 0.909721, 0.919679,
    0.929712, 0.944704, 0.947323, 0.957625, 0.96538,  0.975883, 0.983901,
    0.992077, 1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1,        1,        1,        1,        1,        1,        1,
    1};

TFloat64List fluxcorr_25_4 = {
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074272,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074272,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074272,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074272,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074272,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074272,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074272,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074272,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074272,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074272,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074272,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,  0.074438,
    0.074438,  0.074438,  0.074272,  0.0741085, 0.073948,  0.0737895, 0.073633,
    0.0734795, 0.073328,  0.0731795, 0.073033,  0.0728895, 0.072748,  0.0726085,
    0.072472,  0.072338,  0.072206,  0.0720765, 0.071949,  0.071824,  0.071702,
    0.071581,  0.071464,  0.0713485, 0.071235,  0.071124,  0.071016,  0.07091,
    0.070806,  0.070704,  0.070605,  0.070508,  0.070413,  0.070321,  0.070231,
    0.0701425, 0.070057,  0.069974,  0.069893,  0.069814,  0.069738,  0.0696635,
    0.069592,  0.069522,  0.069455,  0.0693895, 0.069327,  0.0692655, 0.069207,
    0.069151,  0.069097,  0.069046,  0.068996,  0.0689495, 0.068904,  0.068862,
    0.068821,  0.0687835, 0.068747,  0.068714,  0.068683,  0.0686535, 0.068627,
    0.0686025, 0.06858,   0.06856,   0.068542,  0.068527,  0.068513,  0.0685025,
    0.068494,  0.0684875, 0.068484,  0.068482,  0.068483,  0.0684855, 0.068491,
    0.068498,  0.068508,  0.06852,   0.068535,  0.0685515, 0.068571,  0.0685925,
    0.068616,  0.0686425, 0.068671,  0.0687025, 0.068736,  0.0687715, 0.06881,
    0.068851,  0.068894,  0.0689395, 0.068987,  0.0690375, 0.069091,  0.069146,
    0.069204,  0.069264,  0.069327,  0.0693925, 0.06946,   0.069531,  0.069604,
    0.06968,   0.069758,  0.0698385, 0.069922,  0.070008,  0.070097,  0.070188,
    0.070282,  0.070379,  0.070478,  0.07058,   0.070685,  0.0707925, 0.070903,
    0.0710155, 0.071132,  0.0712505, 0.071372,  0.071496,  0.071624,  0.071754,
    0.071887,  0.072023,  0.072162,  0.0723035, 0.072449,  0.0725965, 0.072748,
    0.0729015, 0.073059,  0.073219,  0.073382,  0.073549,  0.073719,  0.073891,
    0.074068,  0.074247,  0.07443,   0.0746155, 0.074805,  0.074998,  0.075194,
    0.0753935, 0.075596,  0.075803,  0.076013,  0.0762265, 0.076443,  0.076664,
    0.076888,  0.0771165, 0.077348,  0.077584,  0.077823,  0.078066,  0.078313,
    0.078564,  0.078819,  0.0790775, 0.07934,   0.0796065, 0.079878,  0.0801525,
    0.080432,  0.080715,  0.081003,  0.081295,  0.081591,  0.0818915, 0.082196,
    0.0825055, 0.08282,   0.0831385, 0.083462,  0.08379,   0.084122,  0.0844595,
    0.084802,  0.085149,  0.085501,  0.0858585, 0.08622,   0.0865875, 0.08696,
    0.0873375, 0.08772,   0.0881085, 0.088502,  0.0889015, 0.089306,  0.089716,
    0.090132,  0.0905535, 0.090981,  0.091414,  0.091853,  0.092298,  0.092749,
    0.0932065, 0.09367,   0.0941405, 0.094617,  0.0950995, 0.095589,  0.0960845,
    0.096587,  0.0970965, 0.097613,  0.0981365, 0.098667,  0.099204,  0.099749,
    0.100301,  0.10086,   0.1014275, 0.102002,  0.1025845, 0.103174,  0.103773,
    0.104379,  0.104993,  0.105615,  0.1062465, 0.106886,  0.107534,  0.10819,
    0.108856,  0.10953,   0.1102135, 0.110906,  0.111608,  0.11232,   0.1130405,
    0.113771,  0.1145115, 0.115262,  0.1160235, 0.116794,  0.1175755, 0.118367,
    0.1191695, 0.119983,  0.1208075, 0.121643,  0.1224905, 0.123349,  0.124219,
    0.125101,  0.1259955, 0.126901,  0.12782,   0.128751,  0.129695,  0.130652,
    0.1316215, 0.132605,  0.1336015, 0.134612,  0.1356365, 0.136675,  0.1377275,
    0.138795,  0.1398775, 0.140975,  0.142087,  0.143215,  0.1443585, 0.145518,
    0.146694,  0.147886,  0.1490955, 0.150321,  0.1515645, 0.152825,  0.154104,
    0.1554,    0.156715,  0.158049,  0.1594015, 0.160773,  0.162164,  0.163575,
    0.165007,  0.166458,  0.167931,  0.169425,  0.1709405, 0.172477,  0.1740375,
    0.175619,  0.177224,  0.178852,  0.1805045, 0.18218,   0.1838815, 0.185607,
    0.187358,  0.189134,  0.1909375, 0.192767,  0.194624,  0.196508,  0.19842,
    0.20036,   0.20233,   0.204329,  0.206358,  0.208417,  0.210508,  0.21263,
    0.214784,  0.216971,  0.2191905, 0.221444,  0.2237325, 0.226055,  0.2284135,
    0.230809,  0.23324,   0.235709,  0.2382165, 0.240762,  0.243348,  0.245973,
    0.2486405, 0.251348,  0.2540995, 0.256893,  0.2597305, 0.262613,  0.2655415,
    0.268515,  0.271538,  0.274607,  0.277726,  0.280894,  0.2841135, 0.287384,
    0.2907085, 0.294085,  0.297518,  0.301004,  0.3045495, 0.308151,  0.311812,
    0.315532,  0.3193145, 0.323158,  0.327066,  0.331037,  0.3350755, 0.339179,
    0.343353,  0.347594,  0.3519085, 0.356293,  0.3607535, 0.365288,  0.3698995,
    0.374588,  0.379358,  0.384207,  0.3891405, 0.394156,  0.3992605, 0.404449,
    0.40973,   0.4151,    0.4205655, 0.426123,  0.4317795, 0.437532,  0.443388,
    0.449344,  0.4554075, 0.461574,  0.4678525, 0.47424,   0.480743,  0.487359,
    0.4940955, 0.50095,   0.507929,  0.515031,  0.5222645, 0.529624,  0.5371215,
    0.544751,  0.552523,  0.560433,  0.5684905, 0.576692,  0.585049,  0.593555,
    0.6022225, 0.611046,  0.620037,  0.629192,  0.6385205, 0.648021,  0.657701,
    0.667561,  0.6776105, 0.687845,  0.698278,  0.708904,  0.7197365, 0.730772,
    0.742022,  0.753483,  0.765169,  0.777075,  0.7892165, 0.798279,  0.798141,
    0.799053,  0.8001325, 0.801556,  0.802042,  0.802845,  0.8040795, 0.805918,
    0.806903,  0.807897,  0.807164,  0.811202,  0.810483,  0.809762,  0.815367,
    0.814665,  0.813961,  0.813256,  0.821197,  0.820519,  0.8198385, 0.819157,
    0.818474,  0.817789,  0.817103,  0.828991,  0.828346,  0.8277,    0.827052,
    0.826403,  0.825752,  0.8251,    0.8244465, 0.823791,  0.823135,  0.822477,
    0.839955,  0.839955,  0.8393585, 0.83876,   0.8381605, 0.83756,   0.836957,
    0.836353,  0.8357485, 0.835142,  0.834534,  0.833924,  0.8333135, 0.832701,
    0.832088,  0.831473,  0.830857,  0.830239,  0.8296205, 0.829,     0.828378,
    0.827755,  0.8271305, 0.8573205, 0.8573205, 0.856801,  0.8562805, 0.855758,
    0.855235,  0.85471,   0.8541845, 0.853658,  0.8531295, 0.852601,  0.85207,
    0.851538,  0.8510055, 0.850471,  0.8499365, 0.8494,    0.848862,  0.848323,
    0.847783,  0.847242,  0.8466995, 0.846156,  0.845611,  0.845065,  0.8445175,
    0.843969,  0.8434195, 0.842869,  0.8423165, 0.841763,  0.8412085, 0.840653,
    0.840096,  0.839538,  0.838979,  0.838418,  0.8378565, 0.837294,  0.8367295,
    0.836164,  0.835598,  0.83503,   0.8344615, 0.833891,  0.8333195, 0.832747,
    0.8321735, 0.831599,  0.8310225, 0.830445,  0.8298665, 0.829287,  0.890402,
    0.890402,  0.8900285, 0.889655,  0.88928,   0.888904,  0.888528,  0.888151,
    0.8877725, 0.887394,  0.887014,  0.886633,  0.8862515, 0.88587,   0.8854865,
    0.885103,  0.884718,  0.884332,  0.883946,  0.883559,  0.8831705, 0.882782,
    0.882392,  0.882001,  0.88161,   0.881218,  0.8808245, 0.880431,  0.8800355,
    0.87964,   0.879244,  0.878847,  0.878449,  0.87805,   0.87765,   0.877249,
    0.8768475, 0.876446,  0.8760425, 0.875638,  0.875234,  0.874828,  0.8744215,
    0.874014,  0.873606,  0.873197,  0.8727875, 0.872377,  0.871965,  0.871553,
    0.8711395, 0.870725,  0.8703105, 0.869895,  0.8694785, 0.869061,  0.8686425,
    0.868223,  0.867803,  0.867383,  0.866961,  0.866538,  0.8661155, 0.865691,
    0.865266,  0.86484,   0.864413,  0.863986,  0.8635575, 0.863128,  0.862698,
    0.862267,  0.861835,  0.861403,  0.860969,  0.860535,  0.8600995, 0.859664,
    0.8592265, 0.858789,  0.8583505, 0.857911,  0.857471,  0.85703,   0.8565875,
    0.856145,  0.8557015, 0.855257,  0.8548115, 0.854365,  0.853918,  0.85347,
    0.853021,  0.852572,  0.8521215, 0.85167,   0.8512175, 0.850765,  0.850311,
    0.849856,  0.8494005, 0.848944,  0.848487,  0.848029,  0.84757,   0.84711,
    0.8466495, 0.846188,  0.8457255, 0.845262,  0.8447985, 0.844333,  0.8438675,
    0.843401,  0.8429335, 0.842465,  0.841996,  0.841526,  0.8410555, 0.840584,
    0.840111,  0.839638,  0.8391635, 0.838688,  0.838212,  0.837736,  0.837258,
    0.836779,  0.8363005, 0.83582,   0.835339,  0.834857,  0.834374,  0.833891,
    0.8334065, 0.832921,  0.8324345, 0.831948,  0.83146,   0.830971,  0.830482,
    0.829992,  0.8295005, 0.829008,  0.8285155, 0.828022,  0.827527,  0.827032,
    0.8265355, 0.826038,  0.8255405, 0.825042,  0.824542,  0.824041,  0.82354,
    0.823038,  0.822535,  0.822031,  0.8215265, 0.821021,  0.8205145, 0.820007,
    0.819499,  0.81899,   0.8184805, 0.81797,   0.8174585, 0.816946,  0.816433,
    0.815919,  0.8154045, 0.814889,  0.814372,  0.813855,  0.813337,  0.812818,
    0.812298,  0.811777,  0.8112555, 0.810733,  0.81021,   0.809686,  0.809161,
    0.808636,  0.8081095, 0.807582,  0.8070535, 0.806525,  0.805995,  0.805464,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1,         1,         1,         1,         1,         1,         1,
    1};

//-----------------------------------------------------------------------------
bool correctMessage(const GlobalException &ex) {
  BOOST_CHECK_EQUAL(
      ex.what(),
      std::string("Cannot convolve: either kernel or array is empty. "));
  return true;
}

CSpectrumFluxCorrectionMeiksin
fakeFluxCorrectionMeiksin(TFloat64List igmLambdas, TFloat64List fluxcorr) {
  std::vector<MeiksinCorrection> corrections;
  corrections.push_back(MeiksinCorrection(igmLambdas, {fluxcorr}));
  CSpectrumFluxCorrectionMeiksin fluxMeiksinObj(corrections);
  return fluxMeiksinObj;
}

BOOST_AUTO_TEST_CASE(Convolve_test) {
  TFloat64List arr = {1, 4, 2, 5}, kernel = {3, 4, 1};
  CSpectrumFluxCorrectionMeiksin spc = fakeFluxCorrectionMeiksin(arr, arr);
  TFloat64List conv = spc.convolve(arr, kernel);
  TFloat64List expectedRes = {8, 21, 25, 27};
  BOOST_CHECK_EQUAL_COLLECTIONS(conv.begin(), conv.end(), expectedRes.begin(),
                                expectedRes.end());
}

BOOST_AUTO_TEST_CASE(Convolve_test_empty) {
  TFloat64List arr = {}, kernel = {3, 4, 1};
  CSpectrumFluxCorrectionMeiksin spc = fakeFluxCorrectionMeiksin(arr, arr);
  BOOST_CHECK_EXCEPTION(spc.convolve(arr, kernel), GlobalException,
                        correctMessage);
}

BOOST_AUTO_TEST_CASE(Convolve_test_identity) {
  TFloat64List arr = {1, 4, 2, 5}, kernel = {1};
  CSpectrumFluxCorrectionMeiksin spc = fakeFluxCorrectionMeiksin(arr, arr);
  TFloat64List conv = spc.convolve(arr, kernel);
  BOOST_CHECK_EQUAL_COLLECTIONS(conv.begin(), conv.end(), conv.begin(),
                                conv.end());
}

BOOST_AUTO_TEST_CASE(correction_Convolve_test) {
  // create a gaussian kernel
  Float64 sigma = 1.;
  TFloat64List kernel = {0, 1, 0};
  TFloat64List arr = {1, 4, 2, 5};
  CSpectrumFluxCorrectionMeiksin spc = fakeFluxCorrectionMeiksin(arr, arr);
  TFloat64List conv = spc.convolve(fluxcorr, kernel);

  BOOST_CHECK_EQUAL(conv.size(), fluxcorr.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(fluxcorr.begin(), fluxcorr.end(), conv.begin(),
                                conv.end());
}

BOOST_AUTO_TEST_CASE(correction_multiply_test) {
  // create a gaussian kernel
  Float64 sigma = 1.;
  Float64 z_center = 2 - (2 - 1.5) / 2;

  std::string lsfType = "GaussianConstantWidth";
  Float64 width = 1.09;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.width", width);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantWidthArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  TFloat64List igmLambdas = _lbdaRange.SpreadOver(lbdaStep);
  CSpectrumFluxCorrectionMeiksin fluxMeiksinObj =
      fakeFluxCorrectionMeiksin(igmLambdas, fluxcorr);
  fluxMeiksinObj.m_convolRange = _lbdaRange;
  CRange<Float64> lbdaRange(fluxMeiksinObj.getLambdaMin(),
                            fluxMeiksinObj.getLambdaMax()); // 200..1299
  BOOST_CHECK_EQUAL(igmLambdas.size(), fluxcorr.size());

  TFloat64List conv =
      fluxMeiksinObj.applyAdaptativeKernel(fluxcorr, z_center, lsf, igmLambdas);
  /*
    FILE* fspc = fopen( "convolvedIGM.csv", "w+" );
    fprintf(fspc, "lambdas,IGMCurve,IGMConvolvedCurve\n");
    for (Int32 i = 0; i < igmLambdas.size(); i++)
    {
        fprintf( fspc, "%f,%f,%f\n", igmLambdas[i], fluxcorr[i], conv[i]);
    }
    fclose( fspc );
  */
  BOOST_CHECK_EQUAL(conv.size(), fluxcorr.size());
  // BOOST_CHECK_EQUAL_COLLECTIONS(fluxcorr.begin(), fluxcorr.end(),
  // conv.begin(), conv.end());
}
BOOST_AUTO_TEST_CASE(correction_multiply_test_CteResolution) {
  // create a gaussian kernel
  Float64 sigma = 1.;
  Float64 z_center = 2 - (2 - 1.5) / 2;

  std::string lsfType = "GaussianConstantResolution";
  Float64 resolution = 2350;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.resolution", resolution);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  TFloat64List igmLambdas = _lbdaRange.SpreadOver(lbdaStep);
  CSpectrumFluxCorrectionMeiksin fluxMeiksinObj =
      fakeFluxCorrectionMeiksin(igmLambdas, fluxcorr);
  fluxMeiksinObj.m_convolRange = _lbdaRange;
  CRange<Float64> lbdaRange(fluxMeiksinObj.getLambdaMin(),
                            fluxMeiksinObj.getLambdaMax()); // 200..1299
  BOOST_CHECK_EQUAL(igmLambdas.size(), fluxcorr.size());

  TFloat64List conv =
      fluxMeiksinObj.applyAdaptativeKernel(fluxcorr, z_center, lsf, igmLambdas);
  /*
    FILE* fspc = fopen( "convolvedIGM_cteResolution.csv", "w+" );
    fprintf(fspc, "lambdas,IGMCurve,IGMConvolvedCurve\n");
    for (Int32 i = 0; i < igmLambdas.size(); i++)
    {
        fprintf( fspc, "%f,%f,%f\n", igmLambdas[i], fluxcorr[i], conv[i]);
    }
    fclose( fspc );
  */
  BOOST_CHECK_EQUAL(conv.size(), fluxcorr.size());
  // BOOST_CHECK_EQUAL_COLLECTIONS(fluxcorr.begin(), fluxcorr.end(),
  // conv.begin(), conv.end());
}

BOOST_AUTO_TEST_CASE(correction_multiply_test_CteResolution25_4) {
  // create a gaussian kernel
  Float64 sigma = 1.;
  Float64 z_center = 2.5 - (0.5) / 2;

  std::string lsfType = "GaussianConstantResolution";
  Float64 resolution = 4300;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.resolution", resolution);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  TFloat64List igmLambdas = _lbdaRange.SpreadOver(lbdaStep);
  CSpectrumFluxCorrectionMeiksin fluxMeiksinObj =
      fakeFluxCorrectionMeiksin(igmLambdas, fluxcorr_25_4);
  fluxMeiksinObj.m_convolRange = _lbdaRange;
  CRange<Float64> lbdaRange(fluxMeiksinObj.getLambdaMin(),
                            fluxMeiksinObj.getLambdaMax()); // 200..1299

  BOOST_CHECK_EQUAL(igmLambdas.size(), fluxcorr_25_4.size());

  TFloat64List conv = fluxMeiksinObj.applyAdaptativeKernel(
      fluxcorr_25_4, z_center, lsf, igmLambdas);

  FILE *fspc = fopen("convolvedIGM_cteResolution_Z25_Curv4.csv", "w+");
  fprintf(fspc, "lambdas,IGMCurve,IGMConvolvedCurve\n");
  for (Int32 i = 0; i < igmLambdas.size(); i++) {
    fprintf(fspc, "%f,%f,%f\n", igmLambdas[i], fluxcorr_25_4[i], conv[i]);
  }
  fclose(fspc);

  BOOST_CHECK_EQUAL(conv.size(), fluxcorr_25_4.size());
  // BOOST_CHECK_EQUAL_COLLECTIONS(fluxcorr_25_4.begin(), fluxcorr_25_4.end(),
  // conv.begin(), conv.end());
}
/*  CRange<Float64> range(-3.*sigma, 3.*sigma); Float64 step = 1.;
  TFloat64List values = range.SpreadOver(step), conv, kernel;
  Float64 sigmaKernel = 1;

  for(Float64 x : values)
      kernel.push_back(
  1/(sqrt(2*M_PI)*sigmaKernel)*exp(-0.5*x*x/sigmaKernel/sigmaKernel) );
*/

BOOST_AUTO_TEST_CASE(correction_multiply_test_CteResolution25_4_incontext) {
  // create a gaussian kernel
  Float64 sigma = 1.;
  Float64 z_center = 2.5 - (0.5) / 2;

  std::string lsfType = "GaussianConstantResolution";
  Float64 resolution = 4300;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.resolution", resolution);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  TFloat64List igmLambdas = _lbdaRange.SpreadOver(lbdaStep);
  CSpectrumFluxCorrectionMeiksin fluxMeiksinObj =
      fakeFluxCorrectionMeiksin(igmLambdas, fluxcorr_25_4);
  CRange<Float64> lbdaRange(fluxMeiksinObj.getLambdaMin(),
                            fluxMeiksinObj.getLambdaMax()); // 200..1299
  BOOST_CHECK_EQUAL(igmLambdas.size(), fluxcorr_25_4.size());
  fluxMeiksinObj.m_convolRange = _lbdaRange;
  fluxMeiksinObj.m_corrections.resize(1);
  fluxMeiksinObj.m_corrections[0].lbda = igmLambdas;
  fluxMeiksinObj.m_corrections[0].fluxcorr.push_back(fluxcorr_25_4);

  TFloat64List conv = fluxMeiksinObj.applyAdaptativeKernel(
      fluxcorr_25_4, z_center, lsf, igmLambdas);
  /*
    FILE* fspc = fopen( "convolvedIGM_cteResolution_Z25_Curv4_context.csv", "w+"
    ); fprintf(fspc, "lambdas,IGMCurve,IGMConvolvedCurve\n"); for (Int32 i = 0;
    i < igmLambdas.size(); i++)
    {
        fprintf( fspc, "%f,%f,%f\n", igmLambdas[i], fluxcorr_25_4[i], conv[i]);
    }
    fclose( fspc );
  */
  BOOST_CHECK_EQUAL(conv.size(), fluxcorr_25_4.size());
  // BOOST_CHECK_EQUAL_COLLECTIONS(fluxcorr_25_4.begin(), fluxcorr_25_4.end(),
  // conv.begin(), conv.end());
}

BOOST_AUTO_TEST_CASE(correction_test) {
  // create a gaussian kernel
  Float64 sigma = 1.;
  Float64 z_center = 2.5 - (0.5) / 2;

  std::string lsfType = "GaussianConstantResolution";
  Float64 resolution = 4300;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.resolution", resolution);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  CRange<Float64> lbdaRange(200, 210);
  TFloat64List igmLambdas = lbdaRange.SpreadOver(lbdaStep);
  TFloat64List fluxsim(igmLambdas.size());
  Int32 count = fluxsim.size();
  for (Int32 i = 0; i < count; i++) {
    fluxsim[i] = exp(-0.5 * i * i / 2 / 2);
  }
  BOOST_CHECK_EQUAL(igmLambdas.size(), fluxsim.size());

  CSpectrumFluxCorrectionMeiksin fluxMeiksinObj =
      fakeFluxCorrectionMeiksin(igmLambdas, {fluxsim});
  fluxMeiksinObj.m_convolRange = _lbdaRange;
  fluxMeiksinObj.m_corrections.resize(1);
  fluxMeiksinObj.m_corrections[0].lbda = igmLambdas;
  fluxMeiksinObj.m_corrections[0].fluxcorr.push_back(fluxsim);

  TFloat64List conv =
      fluxMeiksinObj.applyAdaptativeKernel(fluxsim, z_center, lsf, igmLambdas);
  /*
    FILE* fspc = fopen( "convolvedRandom.csv", "w+" );
    fprintf(fspc, "lambdas,IGMCurve,IGMConvolvedCurve\n");
    for (Int32 i = 0; i < igmLambdas.size(); i++)
    {
        fprintf( fspc, "%f,%f,%f\n", igmLambdas[i], fluxsim[i], conv[i]);
    }
    fclose( fspc );
  */
  BOOST_CHECK_EQUAL(conv.size(), fluxsim.size());
  // BOOST_CHECK_EQUAL_COLLECTIONS(fluxcorr_25_4.begin(), fluxcorr_25_4.end(),
  // conv.begin(), conv.end());
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

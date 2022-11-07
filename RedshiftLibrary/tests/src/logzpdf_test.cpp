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
#include "RedshiftLibrary/operator/logZPdfResult.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(logzpdfTest)

BOOST_AUTO_TEST_CASE(mixedGrid_test) {

  TZGridListParams zparams(2);
  zparams[0] = ZGridParameters(0, 10, 3, NAN);
  zparams[1] = ZGridParameters(2, 5, 1, 3);
  const TFloat64List mixedGrid =
      CLogZPdfResult::buildLogMixedZPdfGrid(false, zparams);

  TFloat64List resultingVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_test_zcenterNotincluded) {

  TZGridListParams zparams(2);
  zparams[0] = ZGridParameters(0, 10, 3, NAN); // 0, 3, 6, 9
  zparams[1] = ZGridParameters(2, 5, 1, 4);
  const TFloat64List mixedGrid =
      CLogZPdfResult::buildLogMixedZPdfGrid(false, zparams);

  TFloat64List resultingVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_test_zcenterNotincluded2) {

  TZGridListParams zparams(2);
  zparams[0] = ZGridParameters(1, 10, 3, NAN); // 1, 4, 7, 10
  zparams[1] = ZGridParameters(2, 5, 0.8, 4);
  const TFloat64List mixedGrid =
      CLogZPdfResult::buildLogMixedZPdfGrid(false, zparams);

  TFloat64List resultingVect{1, 2.4, 3.2, 4, 4.8, 7, 10};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(finegrid) {

  TZGridListParams zparams(2);
  zparams[0] = ZGridParameters(1, 10, 3, NAN);  // 1, 4, 7, 10
  zparams[1] = ZGridParameters(2.5, 5, 1.5, 4); // 1, 2.5, 4, 7, 10
  TFloat64List valproba = {10, 15, 20, 30, 40};
  TFloat64List resulting_zgridFine = {1, 2.5, 4, 5.5, 7, 8.5, 10};
  TFloat64List resulting_proba = {10, 15, 20, 25, 30, 35, 40};
  const TPdf _pdf = CLogZPdfResult::getLogZPdf_fine(false, zparams, valproba);
  BOOST_CHECK(_pdf.zgrid.size() == _pdf.probaLog.size());

  BOOST_CHECK(_pdf.zgrid.size() == resulting_zgridFine.size());
  BOOST_CHECK(_pdf.zgrid == resulting_zgridFine);

  BOOST_CHECK(_pdf.probaLog.size() == resulting_proba.size());
  BOOST_CHECK(_pdf.probaLog == resulting_proba);
}
BOOST_AUTO_TEST_CASE(iscoherent) {
  const TFloat64List zstep = {0.01, 0.001};
  BOOST_CHECK(CLogZPdfResult::isZGridCoherent(zstep) == false);
  const TFloat64List zstep_coherent = {0.001, 0.001};
  BOOST_CHECK(CLogZPdfResult::isZGridCoherent(zstep_coherent) == true);
  const TFloat64List zstep_coherent_single = {0.001};
  BOOST_CHECK(CLogZPdfResult::isZGridCoherent(zstep_coherent_single) == true);
}

BOOST_AUTO_TEST_CASE(genericBuildZPDFGrid) {
  TZGridListParams zparams(2);
  zparams[0] = ZGridParameters(0, 10, 3, NAN); // 0, 3, 6, 9
  zparams[1] = ZGridParameters(2, 5, 1, 4);
  const TFloat64List coarseGrid =
      CLogZPdfResult::buildLogZPdfGrid(false, zparams, COARSE);

  TFloat64List resultingCoarseVect{0, 3, 6, 9};
  BOOST_CHECK(coarseGrid.size() == resultingCoarseVect.size());
  BOOST_CHECK(coarseGrid == resultingCoarseVect);

  const TFloat64List fineGrid =
      CLogZPdfResult::buildLogZPdfGrid(false, zparams, FINE);
  TFloat64List resultingFineVect{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  BOOST_CHECK(fineGrid.size() == resultingFineVect.size());
  BOOST_CHECK(fineGrid == resultingFineVect);

  const TFloat64List mixedGrid =
      CLogZPdfResult::buildLogZPdfGrid(false, zparams, MIXED);
  TFloat64List resultingMixedVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingMixedVect.size());
  BOOST_CHECK(mixedGrid == resultingMixedVect);
}

BOOST_AUTO_TEST_SUITE_END()
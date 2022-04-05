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
#include "RedshiftLibrary/operator/tplcombinationresult.h"
#include <cfloat>
using namespace NSEpic;

void CTplCombinationResult::Init(Int32 n, Int32 EbmvListSize,
                                 Int32 MeiksinListSize, Int32 componentSize) {
  ChiSquare.resize(n);
  ChiSquarePhot.resize(n);
  FitEbmvCoeff.resize(n);
  FitMeiksinIdx.resize(n);
  FitCOV.resize(n); // covariance
  // LogPrior.resize(n);
  Redshifts.resize(n);
  Overlap.resize(n);
  Status.resize(n);
  SNR.resize(n);

  ChiSquareIntermediate.clear();
  IsmEbmvCoeffIntermediate.clear();
  IgmMeiksinIdxIntermediate.clear();

  std::vector<TFloat64List> _chi2ListList(
      EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX));
  std::vector<TFloat64List> _ismListList(EbmvListSize,
                                         TFloat64List(MeiksinListSize, NAN));
  std::vector<TInt32List> _igmListList(EbmvListSize,
                                       TInt32List(MeiksinListSize, -1));

  ChiSquareIntermediate =
      std::vector<std::vector<TFloat64List>>(n, _chi2ListList);
  IsmEbmvCoeffIntermediate =
      std::vector<std::vector<TFloat64List>>(n, _ismListList);
  IgmMeiksinIdxIntermediate =
      std::vector<std::vector<TInt32List>>(n, _igmListList);

  FitAmplitude = std::vector<TFloat64List>(n, TFloat64List(componentSize, NAN));
  FitAmplitudeError =
      std::vector<TFloat64List>(n, TFloat64List(componentSize, NAN));
  FitAmplitudeSigma =
      std::vector<TFloat64List>(n, TFloat64List(componentSize, NAN));
}

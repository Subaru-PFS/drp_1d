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
#ifndef _REDSHIFT_OPERATOR_LINEMODELRESULT_
#define _REDSHIFT_OPERATOR_LINEMODELRESULT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/processflow/result.h"

#include "RedshiftLibrary/continuum/indexes.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/continuummodelsolution.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/linemodel/linemodelsolution.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/priorhelper.h"
#include <memory>

namespace NSEpic {

class CTemplatesFitStore;

class CLineModelResult : public COperatorResult {
public:
  CLineModelResult() = default;

  void Init(TFloat64List redshifts, CLineCatalog::TLineVector restLines,
            Int32 nTplContinuum, Int32 nTplratios,
            TFloat64List tplratiosPriors);

  Int32 getNLinesOverCutThreshold(Int32 solutionIdx, Float64 snrThres,
                                  Float64 fitThres) const;
  TBoolList getStrongLinesPresence(
      Int32 filterType,
      const std::vector<CLineModelSolution> &linemodelsols) const;
  TBoolList getStrongestLineIsHa(
      const std::vector<CLineModelSolution> &linemodelsols) const;
  TInt32List getNLinesAboveSnrcut(
      const std::vector<CLineModelSolution> &linemodelsols) const;

  Float64 getMinChiSquare() const;
  Float64 getMaxChiSquare() const;
  Int32 ResizeChisquareTplRatios(Int32 nTplratios, Int32 nRedshifts);
  void SetChisquareTplContinuumResult(
      Int32 index,
      const std::shared_ptr<const CTemplatesFitStore> &tplFitStore);
  void SetChisquareTplContinuumResultFromPrevious(Int32 index);
  void
  SetChisquareTplratioResult(Int32 index,
                             std::shared_ptr<CTplratioManager> tplratioManager);
  void
  SetChisquareTplratioResult(Int32 index, const TFloat64List &chisquareTplratio,
                             const TFloat64List &scaleMargCorrTplratio,
                             const TBoolList &strongEmissionLinePresentTplratio,
                             const TBoolList &strongHalphaELPresentTplratios,
                             const TInt32List &nLinesAboveSNRTplratio,
                             const TFloat64List &priorLinesTplratio);

  TFloat64List getChisquareTplContinuumResult(Int32 index_z);
  TFloat64List getChisquareTplratioResult(Int32 index_z);
  TFloat64List getScaleMargCorrTplratioResult(Int32 index_z);
  TBoolList getStrongELPresentTplratioResult(Int32 index_z);
  TBoolList getHaELPresentTplratioResult(Int32 index_z);
  TInt32List getNLinesAboveSNRTplratioResult(Int32 index_z);
  TFloat64List getPriorLinesTplratioResult(Int32 index_z);

  // Merit results
  TFloat64List Redshifts;           // z axis
  TFloat64List ChiSquare;           // min chi2
  TFloat64List ScaleMargCorrection; // margCorrection for min chi2

  std::vector<TFloat64List>
      ChiSquareTplratios;      // full chi2 results (for each tplratio)
  TFloat64List PriorTplratios; // model prior (for each tplratio)
  std::vector<TFloat64List>
      PriorLinesTplratios; // lines priors (for each tplratio)
  std::vector<TFloat64List>
      ScaleMargCorrectionTplratios; // full scale marginalization correction
                                    // results (for each tplratio)
  std::vector<TBoolList>
      StrongELPresentTplratios; // full strongELPresent results (for each
                                // tplratio)
  std::vector<TBoolList>
      StrongHalphaELPresentTplratios; // full strongHalphaPresent results (for
                                      // each tplratio)
  std::vector<TInt32List>
      NLinesAboveSNRTplratios;     // full n_lines_above_snr results (for each
                                   // tplratio)
  TFloat64List ChiSquareContinuum; // chi2 result for the continuum
  TFloat64List
      ScaleMargCorrectionContinuum; //  scale marginalization correction result
                                    //  for the continuum
  std::vector<TFloat64List>
      ChiSquareTplContinuum; // chi2 for all continuum templates fited

  std::vector<CLineModelSolution> LineModelSolutions;
  std::vector<CContinuumModelSolution> ContinuumModelSolutions;

  COperator::TStatusList Status;
  CLineCatalog::TLineVector restLineList;
  Int32 nSpcSamples = 0;
  Float64 dTransposeD = 0.0;
  Float64 cstLog = 0.0;
};

} // namespace NSEpic

#endif

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
#ifndef _REDSHIFT_HYBRID_FITTER_
#define _REDSHIFT_HYBRID_FITTER_

#include "RedshiftLibrary/linemodel/svdfitter.h"

namespace NSEpic

{

// class CRegulament;
class CHybridFitter : public CSvdFitter {
public:
  CHybridFitter(const std::shared_ptr<CLMEltListVector> &elementsVector,
                const CCSpectrumVectorPtr &inputSpcs,
                const CTLambdaRangePtrVector &lambdaRanges,
                const CSpcModelVectorPtr &spectrumModels,
                const CLineMap &restLineList,
                const CSpectraGlobalIndex &spcIndex,
                bool enableAmplitudeOffsets = false,
                bool enableLambdaOffsetsFit = false);

protected:
  virtual void doFit(Float64 redshift) override;
  bool m_opt_enable_improveBalmerFit = false;
  void fitAmplitudesHybrid(Float64 redshift);
  void improveBalmerFit(Float64 redshift);
  virtual bool isIndividualFitEnabled() const {
    return !m_enableAmplitudeOffsets;
  };
  // Initialization of tags
  TStringList initEmissionBalmerTags();
  TStringList initAbsorptionBalmerTags();
  std::vector<TStringList> initAdditionalTags();

  // Validation
  bool validateTagSizes(const TStringList &E, const TStringList &A,
                        const std::vector<TStringList> &More);

  // Main per-pair processing
  void processBalmerPair(Int32 itag, Float64 redshift,
                         const TStringList &linetagsE,
                         const TStringList &linetagsA,
                         const std::vector<TStringList> &linetagsMore);

  // Validity checks
  bool isValidBalmerPair(Int32 iEltE, Int32 iEltA, Int32 lineE_id,
                         Int32 lineA_id);

  // Collect additional lines
  std::pair<TInt32List, TInt32List>
  collectAdditionalLines(Int32 itag,
                         const std::vector<TStringList> &linetagsMore);

  // Width condition
  bool widthConditionSatisfied(Int32 iEltE, Int32 lineE_id, Int32 iEltA,
                               Int32 lineA_id, Float64 redshift);

  // Refit logic
  void attemptBalmerRefit(Int32 iEltA, Int32 lineA_id, Int32 iEltE,
                          Int32 lineE_id, const TInt32List &ilinesMore,
                          const TInt32List &idsMore, Float64 redshift);
  std::pair<Float64, Float64> getAmplitudeAndError(Int32 iElt, Int32 lineId);
  std::pair<TFloat64List, TFloat64List>
  collectAmplitudes(const TInt32List &ilines);
  void
  restoreAmplitudes(Int32 iEltA, Int32 lineA_id, Float64 ampA, Float64 errA,
                    Int32 iEltE, Int32 lineE_id, Float64 ampE, Float64 errE,
                    const TInt32List &ilinesMore, const TInt32List &idsMore,
                    const TFloat64List &ampsMore, const TFloat64List &errsMore);
};
} // namespace NSEpic
#endif

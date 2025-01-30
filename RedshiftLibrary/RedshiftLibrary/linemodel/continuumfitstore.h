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
#ifndef _REDSHIFT_LINEMODEL_CONTINUUMFITSTORE_
#define _REDSHIFT_LINEMODEL_CONTINUUMFITSTORE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/continuumfitstore.h"
#include "RedshiftLibrary/linemodel/continuummodelsolution.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
namespace NSEpic {

class CContinuumFitStore {
public:
  CContinuumFitStore(const TFloat64List &redshifts);

  Int32 GetRedshiftIndex(Float64 z) const;
  Int32 getClosestLowerRedshiftIndex(Float64 z) const;
  const TFloat64List &GetRedshiftList() const;

  virtual const CContinuumModelSolution &
  GetFitValues(Int32 idxz, Int32 continuumCandidateRank) const;
  virtual const CContinuumModelSolution &
  GetFitValues(Float64 redshiftVal, Int32 continuumCandidateRank) const;
  std::vector<std::vector<CContinuumModelSolution>> GetFitValuesVector();
  virtual Int32 getContinuumCount() const = 0;
  void initFitValues();
  std::pair<Float64, CContinuumModelSolution const>
  FindMaxAmplitudeSigma() const;
  CContinuumModelSolution const &FindMinReducedChi2() const;

protected:
  TFloat64List m_redshiftgrid;
  std::vector<std::vector<CContinuumModelSolution>>
      m_fitValues; //[nz][m_nContinuumCandidates]
  virtual Float64
  getFracAmplitudeSigma(CContinuumModelSolution const &continuum) const = 0;

private:
};

} // namespace NSEpic

#endif

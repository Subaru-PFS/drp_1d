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
#ifndef _REDSHIFT_STATISTICS_PRIORHELPER_
#define _REDSHIFT_STATISTICS_PRIORHELPER_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <boost/format.hpp>
#include <cfloat>
#include <string>
#include <vector>

namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CPriorHelper {

public:
  struct SPriorTZE {
    Float64 priorTZE; // p(ebv, tpl | z)
    Float64 A_sigma;
    Float64 A_mean;
    Float64 betaA;
    Float64 betaTE;
    Float64 betaZ;

    Float64 logprior_precompA;  // log(1/sqrt(2 pi sigma_a^2))
    Float64 logprior_precompTE; // log(p(ebv, tpl|z))
    Float64 logprior_precompZ;  // log(p(z)/dz)
  };
  typedef std::vector<SPriorTZE> TPriorEList;
  typedef std::vector<TPriorEList> TPriorZEList;
  typedef std::vector<TPriorZEList> TPriorTZEList;

  CPriorHelper();
  ~CPriorHelper();

  bool Init(std::string priorDirPath, Int32 type);

  bool GetTplPriorData(const std::string &tplname,
                       const TRedshiftList &redshifts,
                       TPriorZEList &zePriorData,
                       Int32 outsideZRangeExtensionMode = 0) const;

  bool GetTZEPriorData(const std::string &tplname, Int32 EBVIndexfilter,
                       Float64 redshift, SPriorTZE &tzePrioData,
                       Int32 outsideZRangeExtensionMode = 0) const;

  bool SetBetaA(Float64 beta);
  bool SetBetaTE(Float64 beta);
  bool SetBetaZ(Float64 beta);

  bool mInitFailed = false;

private:
  bool LoadFileEZ(const char *filePath, std::vector<TFloat64List> &data);
  bool LoadFileZ(const char *filePath, TFloat64List &data);

  bool SetSize(Int32 size);
  bool SetTNameData(Int32 k, std::string tname);
  bool SetEZTData(Int32 k, const std::vector<TFloat64List> &ezt_data);
  bool SetAGaussmeanData(Int32 k,
                         const std::vector<TFloat64List> &agaussmean_data);
  bool SetAGausssigmaData(Int32 k,
                          const std::vector<TFloat64List> &agausssigma_data);
  bool SetPzData(const TFloat64List &z_data);

  Int32 m_type; // 0=Continuum, 1=Lines

  TPriorTZEList m_data;
  TFloat64List m_data_pz;
  TStringList m_tplnames;

  Int32 m_nZ = 24;
  Float64 m_dz = 0.25;
  Float64 m_z0 = 0.0;

  Int32 m_nEbv = 10;

  Float64 m_deltaA = (1e-13 - 1e-20);

  Float64 m_betaTE = -1;
  Float64 m_betaA = -1;
  Float64 m_betaZ = -1;

  Float64 m_priorminval = 0.0; // DBL_MIN;
};

} // namespace NSEpic

#endif

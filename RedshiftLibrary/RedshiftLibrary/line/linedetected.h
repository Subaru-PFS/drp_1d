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
#ifndef _REDSHIFT_LINE_LINEDETECTED_
#define _REDSHIFT_LINE_LINEDETECTED_

#include "RedshiftLibrary/line/line.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 * Represent a Single Line.
 */
class CLineDetected : public CLine {

public:
  CLineDetected() = default;
  CLineDetected(const std::string &name, Float64 pos, EType type,
                CLineProfile_ptr &&profile, EForce force, Float64 amp = -1.0,
                Float64 width = -1.0, Float64 cut = -1.0, Float64 posErr = -1.0,
                Float64 sigmaErr = -1.0, Float64 ampErr = -1.0)
      : CLine(name, pos, type, std::move(profile), force, 0., false, undefStr,
              1.0, undefStr, undefIdx, undefStr),
        m_Amp(amp), m_Width(width), m_Cut(cut), m_PosFitErr(posErr),
        m_SigmaFitErr(sigmaErr), m_AmpFitErr(ampErr){};

  ~CLineDetected() = default;
  CLineDetected(const CLineDetected &other) = default;
  CLineDetected(CLineDetected &&other) = default;
  CLineDetected &operator=(const CLineDetected &other) = default;
  CLineDetected &operator=(CLineDetected &&other) = default;

  bool operator<(const CLineDetected &rhs) const {
    if (m_Pos == rhs.m_Pos)
      return (m_Amp < rhs.m_Amp);
    else
      return (m_Pos < rhs.m_Pos);
  };

  bool operator!=(const CLineDetected &rhs) const {
    if (m_Pos == rhs.m_Pos)
      return (m_Amp != rhs.m_Amp);
    else
      return (m_Pos != rhs.m_Pos);
  };

  Float64 GetAmplitude() const { return m_Amp; };
  Float64 GetWidth() const { return m_Width; };
  Float64 GetCut() const { return m_Cut; };
  Float64 GetPosFitError() const { return m_PosFitErr; };
  Float64 GetSigmaFitError() const { return m_SigmaFitErr; };
  Float64 GetAmpFitError() const { return m_AmpFitErr; };

private:
  Float64 m_Amp = 0.;
  Float64 m_Width = 0.;
  Float64 m_Cut = 0.;

  // fit err
  Float64 m_PosFitErr = 0.;
  Float64 m_SigmaFitErr = 0.;
  Float64 m_AmpFitErr = 0.;
};

using CLineDetectedVector = std::vector<CLineDetected>;
using CLineDetectedMap = std::map<Int32, CLineDetected>;

} // namespace NSEpic

#endif

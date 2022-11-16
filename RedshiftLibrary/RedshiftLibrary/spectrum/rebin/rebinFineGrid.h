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
#ifndef _REDSHIFT_SPECTRUM_REBIN_REBINFINEGRID_
#define _REDSHIFT_SPECTRUM_REBIN_REBINFINEGRID_

#include "RedshiftLibrary/spectrum/rebin/rebin.h"

namespace Rebin {
class rebinFineGrid_test;
} // namespace Rebin
namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CRebinFineGrid : public CRebin {

public:
  using CRebin::CRebin;
  CRebinFineGrid(CRebin &&other) : CRebin(std::move(other)){};

  void rebin(CSpectrumFluxAxis &rebinedFluxAxis, const TFloat64Range &range,
             const CSpectrumSpectralAxis &targetSpectralAxis,
             CSpectrum &rebinedSpectrum, CMask &rebinedMask,
             const std::string opt_error_interp, const TAxisSampleList &Xsrc,
             const TAxisSampleList &Ysrc, const TAxisSampleList &Xtgt,
             const TFloat64List &Error, Int32 &cursor) override;
  void reset() override { clearFineGrid(); };

  const std::string &getType() override { return m_type; };

protected:
  friend class Rebin::rebinFineGrid_test;
  void rebinFineGrid();
  void clearFineGrid();

  const std::string m_type = "precomputedfinegrid";
};

} // namespace NSEpic

#endif
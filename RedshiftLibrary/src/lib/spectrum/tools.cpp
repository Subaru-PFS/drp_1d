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
#include "RedshiftLibrary/spectrum/tools.h"

#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <math.h>

using namespace NSEpic;

CSpectrumTools::CSpectrumTools() {}

CSpectrumTools::~CSpectrumTools() {}

void CSpectrumTools::Interpolate(const CSpectrumAxis &axisXorg,
                                 const CSpectrumAxis &axisYorg, Int32 offsetOrg,
                                 Int32 nOrg, const CSpectrumAxis &axisXint,
                                 CSpectrumAxis &axisYint, CMask &mask) {
  Int32 nInt = axisXint.GetSamplesCount();
  Int32 i; // move over original axis
  Int32 j; // move over interpolated axis

  const Float64 *Xorg = axisXorg.GetSamples() + offsetOrg;
  const Float64 *Yorg = axisYorg.GetSamples() + offsetOrg;
  const Float64 *Xint = axisXint.GetSamples();
  Float64 *Yint = axisYint.GetSamples();

  // While original interpolated axis doesn't overlap original axis,
  // mask this area
  for (i = 0, j = 0; j < nInt && Xint[j] < Xorg[i]; j++) {
    mask[j] = 0;
  }

  // For each sample in original axis
  for (i = 0; i < nOrg - 1; i++) {
    // While spectrum points are between 2 template points
    while (j < nInt && Xint[j] < Xorg[i + 1]) {
      // Perform linear interpolation of template flux value
      Float64 t = (Xint[j] - Xorg[i]) / (Xorg[i + 1] - Xorg[i]);

      TFloat64Range r = TFloat64Range(Yorg[i], Yorg[i + 1]);
      Yint[j] = r.Interpolate(t);
      mask[j] = 1;
      j++;
    }
  }

  for (; j < nInt; j++) {
    // DebugAssert( Xint[j]>Xorg[Norg-1] );

    if (Xint[j] > Xorg[nOrg - 1]) {
      mask[j] = 0;
    } else {
      mask[j] = 1;
      Yint[j] = Yorg[j]; // No interpolation performed
    }
  }
}

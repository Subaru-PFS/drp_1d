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
#include "RedshiftLibrary/operator/peakdetection.h"

#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/operator/peakdetectionresult.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;

CPeakDetection::CPeakDetection(Float64 windowSize, Float64 cut,
                               Int32 medianSmoothHalfWidth, Int32 enlargeRate,
                               Float64 detectionnoiseoffset) {
  m_winsize = windowSize;
  m_cut = cut;
  m_medianSmoothHalfWidth = medianSmoothHalfWidth;
  m_enlargeRate = enlargeRate;
  m_detectionnoiseoffset = detectionnoiseoffset;
}

CPeakDetection::~CPeakDetection() {}

std::shared_ptr<const CPeakDetectionResult>
CPeakDetection::Compute(const CSpectrum &spectrum) {
  auto result = std::make_shared<CPeakDetectionResult>();

  const CSpectrumFluxAxis &fluxAxis = spectrum.GetFluxAxis();
  const CSpectrumSpectralAxis &spectralAxis = spectrum.GetSpectralAxis();

  CSpectrumFluxAxis smoothedFluxAxis = fluxAxis;

  if (m_medianSmoothHalfWidth) {
    smoothedFluxAxis.ApplyMedianSmooth(m_medianSmoothHalfWidth);
  }

  TInt32RangeList peaksBorders;
  FindPossiblePeaks(smoothedFluxAxis, spectralAxis, peaksBorders);

  // No Peak detected, exit
  if (peaksBorders.size() == 0) {
    return nullptr;
  }

  RedefineBorders(peaksBorders, spectralAxis, smoothedFluxAxis, fluxAxis);

  result->PeakList = peaksBorders;
  TInt32RangeList peaksBordersEnlarged = peaksBorders;
  if (m_enlargeRate) {
    for (Int32 i = 0; i < ssize(peaksBorders); i++) {
      TInt32Range fitRange = FindGaussianFitStartAndStop(
          i, peaksBorders, m_enlargeRate, spectralAxis.GetSamplesCount());
      peaksBordersEnlarged[i] = fitRange;
    }
  }
  result->EnlargedPeakList = peaksBordersEnlarged;

  return result;
}

TInt32Range
CPeakDetection::FindGaussianFitStartAndStop(Int32 i,
                                            const TInt32RangeList &peaksBorders,
                                            Int32 enlargeRate, Int32 len) {
  Int32 fitStart = peaksBorders[i].GetBegin();
  Int32 fitStop = peaksBorders[i].GetEnd() + 1;

  Float64 width = fitStop - fitStart;
  fitStart = max(0, fitStart - (int)(enlargeRate * width));
  fitStop = min(len, fitStop + (int)(enlargeRate * width));

  if (i > 0) {
    if (peaksBorders[i - 1].GetEnd() > -1)
      fitStart = max(peaksBorders[i - 1].GetEnd(), fitStart);
  }

  if (i < ssize(peaksBorders) - 1) {
    if (peaksBorders[i + 1].GetBegin() > -1)
      fitStop = min(peaksBorders[i + 1].GetBegin(), fitStop);
  }

  return TInt32Range(fitStart, fitStop);
}

void CPeakDetection::GetNewBorder(const TFloat64List &smoothFluxData,
                                  Int32 &new_border, Int32 &old_border,
                                  bool isRightSide) {
  Int32 border = isRightSide ? smoothFluxData.size() - 1 : 0;
  while (smoothFluxData[new_border] <= smoothFluxData[old_border]) {
    if (new_border == border) {
      old_border = isRightSide ? new_border : 0;
      break;
    }
    old_border = new_border;
    new_border = isRightSide ? max(0, new_border + 1) : max(0, new_border - 1);
  }
}

/**
 * Not fully implemented.
 * Does not seam to do any vital things (CF ez.function.lines.EZELfind: 697)
 *
 * _find_borders8, (CF ez.function.lines.EZELfind: 365)
 * 1. find center = max_position in the current range
 * 2. if center (=max position) is  on the left or on the right -> disable this
 * peak
 * 3. todo: use waves and fluxAxis for proper logging
 */
void CPeakDetection::RedefineBorders(TInt32RangeList &peakList,
                                     const CSpectrumAxis &waves,
                                     const CSpectrumAxis &smoothFluxAxis,
                                     const CSpectrumAxis &fluxAxis) {
  Int32 nPeaksInitial = peakList.size();
  for (Int32 iPeak = nPeaksInitial - 1; iPeak >= 0; iPeak--) {
    int centerPos = -1;
    Float64 centerVal = -1e12;
    // find position of the maximum on the smoothed flux
    Int32 start = peakList[iPeak].GetBegin();
    Int32 stop = peakList[iPeak].GetEnd() + 1;

    for (Int32 i = start; i < stop; i++)
      if (centerVal < smoothFluxAxis[i]) {
        centerVal = smoothFluxAxis[i];
        centerPos = i;
      }

    int old_left = centerPos;
    int old_right = centerPos;
    int new_left = std::max(0, centerPos - 1);
    int new_right =
        std::min((int)smoothFluxAxis.GetSamplesCount(), centerPos + 1);

    if (new_left == 0 || new_right == smoothFluxAxis.GetSamplesCount())
      // if the max is found on the border, then erase this range
      peakList.erase(peakList.begin() + iPeak);
    else {
      GetNewBorder(smoothFluxAxis.GetSamplesVector(), new_left, old_left,
                   false);
      GetNewBorder(smoothFluxAxis.GetSamplesVector(), new_right, old_right,
                   true);
      if (new_right == new_left)
        peakList.erase(peakList.begin() + iPeak);
    }
  }
}

void CPeakDetection::FindPossiblePeaks(
    const CSpectrumAxis &fluxAxis, const CSpectrumSpectralAxis &spectralAxis,
    TInt32RangeList &peakList) {
  peakList.clear();
  Int32 s = fluxAxis.GetSamplesCount();
  TFloat64List med(s);
  TFloat64List xmad(s);

  // Compute median value for each sample over a window of size
  // windowSampleCount
  const TFloat64List &fluxVector = fluxAxis.GetSamplesVector();
  CMedian<Float64> medianFilter;
  for (Int32 i = 0; i < s; i++) {
    // old: regular sampling hypthesis
    // Int32 halfWindowSampleCount = windowSampleCount / 2;
    // Int32 start = std::max( 0, i - halfWindowSampleCount );
    // Int32 stop = std::min( (Int32) fluxAxis.GetSamplesCount(), i +
    // halfWindowSampleCount ); irregular sampling compatible
    Int32 start = std::max(0, spectralAxis.GetIndexAtWaveLength(
                                  spectralAxis[i] - m_winsize / 2.0));
    Int32 stop = std::min(s, spectralAxis.GetIndexAtWaveLength(
                                 spectralAxis[i] + m_winsize / 2.0));

    med[i] = medianFilter.Find(fluxVector.begin() + start,
                               fluxVector.begin() + stop);
    xmad[i] =
        XMad(fluxVector.begin() + start, fluxVector.begin() + stop, med[i]);
    xmad[i] +=
        m_detectionnoiseoffset; // add a noise level, useful for simulation data
  }

  /*//debug:
  // save median and xmad,  flux data
  FILE* f = fopen( "peakdetection_dbg_median.txt", "w+" );
  for( Int32 i=0; i<fluxAxis.GetSamplesCount(); i++ )
  {
      if( med[i] < 0.0001 ){
          fprintf( f, "%e %e %e %e %e\n", spectralAxis[i], med[i], xmad[i],
  med[i]+0.5*m_cut*xmad[i], fluxData[i]); }else{ fprintf( f, "%f %f %f %f %f\n",
  spectralAxis[i], med[i], xmad[i], med[i]+0.5*m_cut*xmad[i], fluxData[i]);
      }
  }
  fclose( f );
  //*/

  // Detect each point whose value is over the median precomputed median
  TInt32List points(s + 1);
  Int32 j = 0;

  for (Int32 i = 0; i < s; i++)
    if (fluxVector[i] > med[i] + 0.5 * m_cut * xmad[i])
      points[j++] = i;

  // No potential peak detected , exit
  if (j == 0)
    return;

  // Find contiguous sample
  Int32 start = points[0];
  Int32 current = points[0];
  Int32 stop = points[0];

  for (Int32 i = 1; i < j + 1; i++) {
    // Index stored in points[] are contiguous (i.e: 4,5,6)
    if (points[i] - 1 == current)
      current = points[i];

    // if we hit a discontinuity, store the previous range of contiguous points
    // representing a potential peak
    else {
      stop = current;
      peakList.push_back(TInt32Range(start, stop));
      start = points[i];
      current = start;
    }
  }
}

Float64 CPeakDetection::XMad(const TFloat64List::const_iterator &begin,
                             const TFloat64List::const_iterator &end,
                             Float64 median) {
  Int32 n = distance(begin, end);
  TFloat64List xdata(n);
  Float64 xmadm = 0.0;

  for (Int32 i = 0; i < n; i++) {
    xdata[i] = fabs(begin[i] - median);
  }

  std::sort(xdata.begin(), xdata.end());

  Int32 n2 = n / 2;
  if (2 * n2 == n) {
    Int32 i1 = n2 - 1;
    Int32 i2 = n2;
    xmadm = 0.5 * (xdata[i1] + xdata[i2]);
  } else {
    xmadm = xdata[n2];
  }

  return xmadm;
}

Float64 CPeakDetection::XMad(const TFloat64List &x, Float64 median) {
  return XMad(x.begin(), x.end(), median);
}

/*
def Xmad(data, xmed):
    """
    The XMAD subroutine calculates the Median Absolute Deviation from
    the sample median. The median, M , is subtracted from each
    ORDERED statistic and then the absolute value is taken. This new
    set of of statistics is then resorted so that they are ORDERED
    statistics. The MAD is then defined to be the median of this
    new set of statistics and is returned as XMADM. The MAD can
    be defined:

                   XMADM = median{ abs(x(i) - M) }

    where the x(i) are the values passed in the array XDATA, and
    the median, M, is passed in the array XLETTER. The set of stats
    in the brackets is assumed to be resorted. For more information
    see page 408 in UREDA.
    """
    ndat = len(data)
    xdata = []
    for item in data:
        xdata.append(abs(item-xmed))
    xdata.sort()
    xdata.insert(0, 0.)
    if (float(ndat)/2. - int(ndat/2.)) == 0:
        i1 = ndat/2
        i2 = ndat/2 + 1
        xmadm = 0.5*(xdata[i1]+xdata[i2])
    else:
        i1 = int(ndat/2) + 1
        xmadm = xdata[i1]
    return xmadm
*/

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
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

#include <cmath>
#include <float.h>
#include <iostream>
#include <numeric>
#include <stdio.h>

#define PEAKS_MIN_THRESHOLD (3)
#define PEAKS_SMOOTH_LIMIT (20)

using namespace NSEpic;
using namespace std;


/**
 * Member attribution constructor.
 */
CExtremum::CExtremum( UInt32 maxPeakCount,
                      Float64 peakSeparation, 
                      Float64 meritcut,
                      bool invertForMinSearch,
                      bool allow_extrema_at_border, 
                      const TFloat64Range& xRange) :
    m_MaxPeakCount( maxPeakCount ),
    m_extrema_separation(peakSeparation),
    m_meritCut(meritcut),
    m_SignSearch(invertForMinSearch ? -1.0 : 1.0 ),
    m_allow_extrema_at_border(allow_extrema_at_border),
    m_PeakSeparationActive(peakSeparation > 0),
    m_XRange( xRange )
{
}

/**
 *  Sets m_XRange to r.
 */
void CExtremum::SetXRange( const TFloat64Range& r )
{
    m_XRange = r;
}

/**
 * Sets m_MaxPeakCount to n.
 */
void CExtremum::SetMaxPeakCount( UInt32 n )
{
    m_MaxPeakCount = n;
}

/**
 * Sets m_meritCut to n.
 */
void CExtremum::SetMeritCut( Float64 n )
{
    m_meritCut = n;
}

/**
 * create an index vector sorting maxY elements
 */
void CExtremum::SortIndexes(TFloat64List&  maxY) const
{
    m_sortedIndexes.resize(maxY.size());
    iota(m_sortedIndexes.begin(), m_sortedIndexes.end(), 0);
    sort(m_sortedIndexes.begin(), m_sortedIndexes.end(),
       [&maxY](Int32 i1, Int32 i2) {return maxY[i1] > maxY[i2];});
}

/**
 * Set a list as default extremum as extreumum, i.e., independent from their values
*/
bool CExtremum::DefaultExtremum( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) 
{
  if(!xAxis.size())
    return false; 
  for (Int32 ke = 0; ke < xAxis.size(); ke++){
      maxPoint.push_back(SPoint(xAxis[ke], yAxis[ke]));
  }
  return true;
}

/**
 * Wrapper around InternalFind, this function validates and sets the search interval limits.
 */
bool CExtremum::Find( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint) const
{
    Int32 n = xAxis.size();
    
    if( n == 0 ){
      throw GlobalException(INTERNAL_ERROR,"CExtremum::Find, input X vector is empty");
    }

    if ( n != yAxis.size()){
      Log.LogError("CExtremum::Find, input X and Y vector have not the same size");
      throw ("CExtremum::Find, input X and Y vector have not the same size");
    }

    Int32 BeginIndex = 0;
    Int32 EndIndex = n-1;

    // Find index in xAxis that correspond to the boundary specified by m_XRange
    if( !m_XRange.GetIsEmpty() )
    {
        bool rangeok;
        rangeok = m_XRange.getClosedIntervalIndices(xAxis, BeginIndex , EndIndex);
        if (!rangeok){
          Log.LogError("CExtremum::Find, bad range [%f, %f] for Xaxis: [%f,%f]", 
              m_XRange.GetBegin(), m_XRange.GetEnd(), xAxis.front(), xAxis.back());
          throw GlobalException(INTERNAL_ERROR,"CExtremum::Find, bad range");
        }
    }
    
    TFloat64List  maxX, minX;
    TFloat64List  maxY, minY;
    bool method = FindAllPeaks( xAxis, yAxis, BeginIndex, EndIndex, maxX, maxY );    

    if(maxX.size() == 0){ 
      // we should not raise an exception here 
      // (we can accept a missing candidate in second pass window)
      // The boolean return has to be tested by the caller
      Log.LogWarning("          CExtremum::Find: FindAllPeaks returned empty MaxX");
      return false;
    }

    //Look for all local minima
    bool activateprominence = false;
    if(activateprominence){
      bool method_min = FindAllPeaks( xAxis, yAxis, BeginIndex, EndIndex, minX, minY, true); // invert signSearch
    
      for(Int32 i = 0; i <minX.size(); i++){
        minY[i] = -1*minY[i];
      }
    
      //Calculate prominence and remove "low" prominence peaks
      //this is ALSO useful for eliminating neighboring peaks in some cases
      // Deactivate Cut_Prominence_Merit for the moment until stat tests are ready
      TFloat64List ret_prominences = Cut_Prominence_Merit(maxX, maxY, minX, minY);    

      if(maxX.size() == 0){
        Log.LogWarning("        CExtremum::Find: Cut_Prominence_Merit returns empty MaxX");
        return false;
      }
    }

    //TODO: add a boolean referring to the metric to use for the sort
    //by default, using the pdf value
    SortIndexes(maxY);
    
    Int32 keepMinN = 1;
    if(m_meritCut>0.0 && maxX.size()>keepMinN){ 
        bool v = Cut_Threshold(maxX, maxY, keepMinN);
    }
    
    //refine: eliminate very close candidates when possible. 
    bool b;
    if(m_PeakSeparationActive)
    {
      b = FilterOutNeighboringPeaksAndTruncate(maxX, maxY, keepMinN, maxPoint);
      
      //verify that peaks are well separated by at least secondpassradius
      bool verified = verifyPeakSeparation(maxPoint);

    }else{
      Truncate(maxX, maxY, maxPoint);
    }
    return true;
}

bool CExtremum::verifyPeakSeparation( TFloat64List& maxX) const{
  bool verified = true;
 
  std::sort(maxX.begin(), maxX.end());
  for(Int32 i = 0; i< maxX.size() -1; i++)  {
    TFloat64Range windowh(maxX[i+1] - m_extrema_separation/2*(1+maxX[i+1]), maxX[i+1] + m_extrema_separation/2*(1+maxX[i+1]));
    windowh.IntersectWith(TFloat64Range(maxX));
    TFloat64Range windowl(maxX[i] - m_extrema_separation/2*(1+maxX[i]), maxX[i] + m_extrema_separation/2*(1+maxX[i]));
    windowl.IntersectWith(TFloat64Range(maxX));
    Float64 overlap;
    overlap = windowh.GetBegin() - windowl.GetEnd();
    if(overlap < 0){
      verified = false;
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"CExtremum::verifyPeakSeparation: Peaks " <<maxX[i]<<" and "<<maxX[i+1] <<" are not enough separated.");
    }
  }
  return verified;
}

bool CExtremum::verifyPeakSeparation( TPointList& maxPoint) const{
  TFloat64List maxX(maxPoint.size());
  for(Int32 i = 0; i<maxPoint.size(); i++)
    maxX[i] = maxPoint[i].X;
  return verifyPeakSeparation(maxX);
}

/**
 * prominence: vertical distance between a summit and the key col, i.e., the closest (horizontal) minima
 * joint prominence_merit cut
 * Method: 
 * 1. identify the key col for consecutive peaks.
 * 2. calculate prominence for the current peak
 * 3. if prominence == peak.y, then peak is a very high peak
 * 4. Remove low prominence peaks : TO decide on "low" value
 * 
 * To identify key col for each peak:
 * 1. Extend to the right and to the left of the peak until reaching higher peaks (ya_current < ya_right_after) //skip peaks of same value and continue extending
 * 2. Find the highest minima within the range identified in (1), considered as the key col.
 * 
 * //returned @prominences could be used as a metric for truncating candidates instead of pdf values
*/
TFloat64List CExtremum::Cut_Prominence_Merit( TFloat64List& maxX, TFloat64List& maxY, TFloat64List& minX, TFloat64List& minY) const
{
  //find highest peak:
  Float64 maxV = *std::max_element(maxY.begin(), maxY.end()); //max of maxY
  Float64 minV = *std::min_element(minY.begin(), minY.end()); //min of minY
  Float64 ref_prominence = maxV - minV; //used to normalize obtained prominence TODO: requires discussion based on results from stats
  Float64 prominence_thresh = 0.1; //heuristic value, TODO: requires discussion (maybe passed in param.json)

  Log.LogDebug("CExtremum::Find: setting prominenceCut to %d and normalization value to %d \n", prominence_thresh, ref_prominence);

  if(maxX.size() == 0){
    throw GlobalException(INTERNAL_ERROR,"    CExtremum::Cut_Prominence_Merit:empty MaxX");
  }

  TFloat64List prominence(maxX.size()), tmpX, tmpY; 

  if(maxX.size() == 1){
    prominence[0] = DBL_MAX; //there is no point of calculating the prominence; Force keeping this only candidate
    return prominence;
  }

  for(Int32 i = 0; i<maxX.size(); i++){
    //current peak elevation
    Float64 ymax = maxY[i];
    TFloat64List rangex, rangey;//To remove: kept for debugging
    Float64 rangex_low = -1, rangex_high = -1;
    //extend to the left until reaching a higher peak
    bool right = true, left = true;
    //extend to the right until reaching a higher peak
    Int32 j = i + 1;
    while(right){
      //a higher peak or last element is reached, leave
      if(maxY[j] > ymax || j == maxX.size()){
        right = false;
        if(j == maxX.size())
          rangex_high = std::max(minX[minX.size() - 1], maxX[maxX.size() - 1]);
        else
          rangex_high = maxX[j];
        break;
      } 
      rangex.push_back(maxX[j]);
      rangey.push_back(maxY[j]);
      j++;
    }
    j = i - 1;
    while(left){
      //a higher peak  or first element is reached, leave
      if(maxY[j] > ymax || j == -1){
        left = false;
        if(j<0)
          rangex_low = std::min(minX[0], maxX[0]);
        else
          rangex_low = maxX[j];
        break;
      }
      rangex.push_back(maxX[j]);
      rangey.push_back(maxY[j]); 
      j--;
    }
    if(rangex_low == -1 || rangex_high == -1){
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Problem in range determination "<< i);
    }
    //look into minX for key_col within rangex
    //key_col is the highest minima between the lowest minimum to the right and left of the peak
    Float64 key_coly_r = DBL_MAX, key_coly_l = DBL_MAX, key_colx;
    Int32 r, l; 
    if(maxX[0]<minX[0]){ //signal starts with a peak
      l = i - 1;
      r = i;
      if(l == -1) //consider only right range
        key_coly_l = -DBL_MAX;
    } else { //signal starts with a minima
      l = i;
      r = i + 1; 
      if(r == minX.size()) //consider only left range
        key_coly_r = -DBL_MAX;  
    }
    //usually there should be as much minima as maxima, for this reason I start j = i and check inclusion
    //case where a peak is the first of last element (needs more work)
    while(r < minX.size()){
      if(minX[r] < rangex_low || minX[r] > rangex_high){ //minima outside range
        break;
      }
      if(minY[r] < key_coly_r){
        key_coly_r = minY[r];
        key_colx = minX[r];//not useful, but just to check!
      }
      r++;
    }

    while(l > -1){
      if(minX[l] < rangex_low || minX[l] > rangex_high){ //minima outside range
        break;
      }
      if(minY[l] < key_coly_l){
        key_coly_l = minY[l];
        key_colx = minX[l];
      }
      l--;
    }
    Float64 tmpProminence = (maxY[i] - std::max(key_coly_l, key_coly_r))/ ref_prominence;
    //keep peaks whose height is almost equal to their prominence
    if(tmpProminence > prominence_thresh || (m_meritCut &&(maxV - maxY[i] < m_meritCut)) ){ //heuristic value
      tmpX.push_back(maxX[i]);
      tmpY.push_back(maxY[i]);
      prominence[i] = tmpProminence;
    }
  }

  maxX.clear(); maxY.clear();
  for(Int32 i = 0; i<tmpX.size(); i++){
    maxX.push_back(tmpX[i]);
    maxY.push_back(tmpY[i]);
  }
  return prominence;
}
/**
 * \Brief: removes extrema based on meritCut.
 * First: Sort based on Y values
 * Second: Remove non-relevant candidates 
 * Third: Sort based on X values to return candidates following their initial order
 * It returns the new list of extrema
*/
bool CExtremum::Cut_Threshold( TFloat64List& maxX, TFloat64List& maxY, Int32 keepMinN) const
{
  if(maxX.size() == 0){
        throw GlobalException(INTERNAL_ERROR,"      CExtremum::Cut_threshold: empty MaxX arg");
  }
  if(maxX.size()<=keepMinN){
    return true; 
  }

  TFloat64List tmpX(maxX), tmpY(maxY);
  maxX.clear(); maxY.clear();
  for(Int32 i = 0; i<tmpX.size(); i++){
    Float64 meritDiff =  tmpY[m_sortedIndexes[0]] - tmpY[i];
    if( meritDiff <= m_meritCut ){
      maxX.push_back(tmpX[i]);
      maxY.push_back(tmpY[i]); 
    }
  }
  if(maxX.size()<keepMinN)
    for (Int32 isort=maxX.size(); isort<keepMinN; isort++){
      Int32 i = m_sortedIndexes[isort];
      auto it=std::lower_bound(maxX.begin(), maxX.end(), tmpX[i]);
      Int32 index = it - maxX.begin();
      maxX.insert(it, tmpX[i]);
      maxY.insert(maxY.begin()+index , tmpY[i]);  
    }

  return true;
}

bool CExtremum::FilterOutNeighboringPeaksAndTruncate(TFloat64List& maxX, TFloat64List& maxY, UInt32 keepmin, TPointList& maxPoint)const
{
  if(maxX.size()<= keepmin){
    for(Int32 i = 0; i<maxX.size(); i++){
      maxPoint.push_back(SPoint(maxX[m_sortedIndexes[i]], m_SignSearch * maxY[m_sortedIndexes[i]]) );
    }
    return true;
  }

  TBoolList peakTracking(maxX.size(), true);
  Float64 wind_high, wind_low;
  Int32 nkeep = 0; 

  for(Int32 i:m_sortedIndexes){
    if (nkeep == m_MaxPeakCount){
        break;
    }

    if(!peakTracking[i]){
      continue;
    }
    nkeep++;
    maxPoint.push_back(SPoint(maxX[i], m_SignSearch * maxY[i]) );
    wind_high = maxX[i] + (maxX[i] + 1) * m_extrema_separation/ (1 - m_extrema_separation/2); 
    wind_low = maxX[i] - (maxX[i] + 1) * m_extrema_separation/ (1 + m_extrema_separation/2);
    TFloat64Range window(wind_low, wind_high);
    Int32 i_min = i, i_max = i;
    bool ret = window.getClosedIntervalIndices(maxX, i_min, i_max, false);    
    for (Int32 j = i_min; j<=i_max; j++){
        if(j == i) continue;
        peakTracking[j] = false;
    }
  }
  return true;
}

bool CExtremum::FindAllPeaks(const TFloat64List & xAxis, const TFloat64List & yAxis,
                             Int32 BeginIndex, Int32 EndIndex, 
                             TFloat64List& maxX, TFloat64List& maxY,
                             bool invertSearch) const
{
  const Float64 SignSearch = invertSearch ? -1.0*m_SignSearch : m_SignSearch;
  TFloat64List tmpY(yAxis);
  for (Float64 & val: tmpY) val *= SignSearch;

  auto goup = [&tmpY](Int32 i){return tmpY[i] < tmpY[i+1];};
  auto godown = [&tmpY](Int32 i){return tmpY[i] > tmpY[i+1];};
  auto goflat = [&tmpY](Int32 i){return tmpY[i] == tmpY[i+1];};

  Int32 maxCount = 0;

  // find first and last non Nan element
  TFloat64List::const_iterator firstNonNan;
  for (; BeginIndex!=EndIndex; ++BeginIndex)
    if (!std::isnan(yAxis[BeginIndex])) break; 
  
  for (; BeginIndex!=EndIndex; --EndIndex)
    if (!std::isnan(yAxis[EndIndex])) break;

  // check at least 3 points left (to get an extrema)
  if ( (EndIndex - BeginIndex) <2 ){
    throw GlobalException(INTERNAL_ERROR,"CExtremum::FindAllPeaks, less than 3 contiguous non nan values");
  }

  bool plank = false; Int32 cnt_plk = 0; 

  // First element
  if (m_allow_extrema_at_border){
    if (godown(BeginIndex)){
      maxX.push_back(xAxis[BeginIndex]);
      maxY.push_back(tmpY[BeginIndex]);
      maxCount++;
    }else if (goflat(BeginIndex)) {
      plank = true;
      cnt_plk++;
    }
  }

  for (Int32 i = BeginIndex + 1; i < EndIndex; i++) {
    if ( goup(i-1)  &&  godown(i)) {
      maxX.push_back(xAxis[i]);
      maxY.push_back(tmpY[i]);
      maxCount++;
      continue;
    }
    if( goup(i-1) && goflat(i)){
      //plank start point
      plank = true;
      cnt_plk++;
      continue;
    }
    if( goflat(i-1)) { 
      //high plank: signal is decreasing after the plank. Peak identified
      if ( godown(i) && plank){ //check if we already identified a high plank
        cnt_plk++;
        Int32 idx_plk = i - round(cnt_plk/2);
        //plank end point
        maxX.push_back(xAxis[idx_plk]);
        maxY.push_back(tmpY[idx_plk]);
        maxCount++; 
        plank = false;
        cnt_plk = 0;
        continue;
      }
      //low plank: pdf is increasing after the plank. No peak here!
      if( goup(i)){
        plank = false; //end
        cnt_plk = 0;
      } else { //Plank is getting larger
        cnt_plk++;
      }
    }
  }

  // last element: check if plank is extended
  if (m_allow_extrema_at_border){
    if ( goup(EndIndex-1) || goflat(EndIndex-1) ) {
      Int32 idx;
      if(plank){
        cnt_plk++;
        idx = EndIndex - round(cnt_plk/2);
      }else{ 
        idx = EndIndex;
      }
      maxX.push_back(xAxis[idx]);
      maxY.push_back(tmpY[idx]);
      maxCount++;
    }
  }

  return  true;
}

/**
 * Brief: Reduce number of peaks based on maxCount passed from param.json if present
*/
bool CExtremum::Truncate( TFloat64List& maxX, TFloat64List& maxY, TPointList& maxPoint) const {
  UInt32 n = maxX.size();
  n = std::min(m_MaxPeakCount, n);
  for (Int32 j = 0; j < n; j++) { 
    maxPoint.push_back(SPoint(maxX[m_sortedIndexes[j]], m_SignSearch * maxY[m_sortedIndexes[j]]) );
  }   
  return true;
}

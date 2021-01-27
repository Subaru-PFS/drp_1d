#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/common/range.h>

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
 * Constructs CExtremum with default values, and SetSignSearch ( -1.0 ) if argument is true, SetSignSearch ( 1.0 ) otherwise.
 */
CExtremum::CExtremum( Bool invertForMinSearch ) :
    m_MaxPeakCount( 10 ),
    m_XRange( 0.0, 0.0 ),
    m_SignSearch( 1.0 ),
    m_extrema_separation(0.005),
    m_meritCut(-2)
{
    if( invertForMinSearch )
      {
        SetSignSearch( -1.0 );
      }
    else
      {
        SetSignSearch( 1.0 );
      }
}

/**
 * Member attribution constructor.
 */
CExtremum::CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount, Float64 peakSeparation, Bool invertForMinSearch, Bool usePeakSeparation ) :
    m_MaxPeakCount( maxPeakCount ),
    m_XRange( xRange ),
    m_SignSearch( 1.0 ),
    m_extrema_separation(peakSeparation),
    m_meritCut(-2),
    m_PeakSeparationActive(usePeakSeparation)
{
    if( invertForMinSearch )
      {
        SetSignSearch( -1.0 );
      }
    else
      {
        SetSignSearch( 1.0 );
      }
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
 * Sets m_SignSearch to val.
 */
void CExtremum::SetSignSearch( Float64 val )
{
    m_SignSearch = val;
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
Bool CExtremum::DefaultExtremum( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) 
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
Bool CExtremum::Find( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint) const
{
    Int32 n = xAxis.size();
    const Float64* selectedXAxis = xAxis.data();
    const Float64* selectedYAxis = yAxis.data();

    Int32 rangeXBeginIndex = -1;
    Int32 rangeXEndIndex = -1;

    if( n == 0 )
        return false;

    // Find index in xAxis that correspond to the boundary specified by m_XRange
    if( m_XRange.GetIsEmpty() )
    {
        rangeXBeginIndex = 0;
        rangeXEndIndex = n-1;
    }
    else
    {
    	// Find index range for the given lambda range
        for( UInt32 i=0; i<n; i++ )
        {
            if( rangeXBeginIndex==-1 && xAxis[i]>=m_XRange.GetBegin() )
            {
                rangeXBeginIndex = i;
            }

            if( xAxis[i]<=m_XRange.GetEnd() )
            {
                rangeXEndIndex = i;
            }
        }
    }
    
    vector < Float64 >  maxX, minX;
    vector < Float64 >  maxY, minY;
    Bool method = FindAllPeaks( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, maxX, maxY );    

    //Look for all local minima
    Bool activateprominence = false;
    if(activateprominence){
      Bool method_min = FindAllPeaks( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, minX, minY, -1*m_SignSearch); 
    
      for(Int32 i = 0; i <minX.size(); i++){
        minY[i] = -1*minY[i];
      }
    }

    if(maxX.size() == 0){
      Log.LogError("          CExtremum::Find: FindAllPeaks returned empty MaxX");
      throw runtime_error("CExtremum::Find: FindAllPeaks returned empty MaxX");
    }

    //Calculate prominence and remove "low" prominence peaks
    //this is ALSO useful for eliminating neighboring peaks in some cases
    // Deactivate Cut_Prominence_Merit for the moment until stat tests are ready
    if(activateprominence){
      TFloat64List ret_prominences = Cut_Prominence_Merit(maxX, maxY, minX, minY);    

      if(maxX.size() == 0){
        Log.LogError("        CExtremum::Find: Cut_Prominence_Merit returns empty MaxX");
        throw runtime_error("CExtremum::Find: Cut_Prominence_Merit returns empty MaxX");
      }
    }

    //TODO: add a boolean referring to the metric to use for the sort
    //by default, using the pdf value
    SortIndexes(maxY);
    
    Int32 keepMinN = 1;
    if(m_meritCut>0.0 && maxX.size()>keepMinN){ 
        Bool v = Cut_Threshold(maxX, maxY, keepMinN);
    }
    //refine: eliminate very close candidates when possible. 
    Bool b;
    if(m_PeakSeparationActive)
      b = FilterOutNeighboringPeaksAndTruncate(maxX, maxY, keepMinN, maxPoint);
    else{
      Truncate(maxX, maxY, maxPoint);
    }
    //verify that peaks are well separated by at least secondpassradius
    Bool verified = verifyPeakSeparation(maxPoint);
    return true;
}

Bool CExtremum::verifyPeakSeparation( TFloat64List& maxX) const{
  Bool verified = true;
 
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
      Log.LogError("  CExtremum::verifyPeakSeparation: Peaks %f and %f are not enough separated.", maxX[i], maxX[i+1]);
      throw runtime_error("  CExtremum::verifyPeakSeparation: two consecutive peaks are not enough separated. Abort");
    }
  }
  return verified;
}

Bool CExtremum::verifyPeakSeparation( TPointList& maxPoint) const{
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
    Log.LogError("    CExtremum::Cut_Prominence_Merit:empty MaxX");
    throw runtime_error("CExtremum::Cut_Prominence_Merit:empty MaxX");
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
      Log.LogError("Problem in range determination %d", i);
      throw runtime_error("Problem in range determination");
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
Bool CExtremum::Cut_Threshold( TFloat64List& maxX, TFloat64List& maxY, Int32 keepMinN) const
{
  if(maxX.size() == 0){
        Log.LogError("      CExtremum::Cut_threshold: empty MaxX arg");
        throw runtime_error(" CExtremum::Cut_threshold: empty MaxX arg");
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

Bool CExtremum::FilterOutNeighboringPeaksAndTruncate(TFloat64List& maxX, TFloat64List& maxY, UInt32 keepmin, TPointList& maxPoint)const
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
    window.IntersectWith(TFloat64Range(maxX));
    Int32 i_min, i_max;
    bool ret = window.getClosedIntervalIndices(maxX, i_min, i_max);    
    for (Int32 j = i_min; j< i; j++){
        peakTracking[j] = false;
    }
    for (Int32 j = i+1; j<=i_max; j++){
        peakTracking[j] = false;
    }
  }
  return true;
}

/**
 * Brief: Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 * Extended to include distance between peaks and a threshold to eliminate peaks
 */
Bool CExtremum::FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, TFloat64List& maxX, TFloat64List& maxY) const {
  if (n == 0)
    return false;

  vector < Float64 > tmpX(n);
  vector < Float64 > tmpY(n);

  for (UInt32 t = 0; t < n; t++) {
    tmpX[t] = xAxis[t];
    tmpY[t] = m_SignSearch * yAxis[t];
  }

  Int32 maxCount = 0;

    // find first and last non Nan element
    Int32 firstNonNanInd = 0;
    for (Int32 iFirst = 0; iFirst < n - 1; iFirst++) {
      if (!std::isnan((double) tmpY[iFirst])) {
        firstNonNanInd = iFirst;
        break;
      }
    }
    Int32 lastNonNanInd = n - 1;
    for (Int32 iLast = n - 1; iLast > 0; iLast--) {
      if (!std::isnan(tmpY[iLast])) {
        lastNonNanInd = iLast; 
        break;
      }
    }
    
    bool plank = false; Int32 cnt_plk = 0; 
    // First element
    if (tmpY[firstNonNanInd] > tmpY[firstNonNanInd + 1]) {
      maxX.push_back(tmpX[firstNonNanInd]);
      maxY.push_back(tmpY[firstNonNanInd]);
      maxCount++;
    }else{
      if(tmpY[firstNonNanInd] == tmpY[firstNonNanInd + 1]) {
        plank = true;
        cnt_plk++;
      }
    }

    for (Int32 i = firstNonNanInd + 1; i < lastNonNanInd; i++) {
     if ((tmpY[i] > tmpY[i - 1]) && (tmpY[i] > tmpY[i + 1])) {
        maxX.push_back(tmpX[i]);
        maxY.push_back(tmpY[i]);
        maxCount++;
        continue;
      }
      if((tmpY[i] > tmpY[i - 1]) && (tmpY[i] == tmpY[i + 1])){
        //plank start point
        plank = true;
        cnt_plk++;
        continue;
      }
      if((tmpY[i] == tmpY[i - 1])){ 
        //high plank: signal is decreasing after the plank. Peak identified
        if (tmpY[i] > tmpY[i + 1] && plank){ //check if we already identified a high plank
          cnt_plk++;
          Int32 idx_plk = i - round(cnt_plk/2);
          //plank end point
          maxX.push_back(tmpX[idx_plk]);
          maxY.push_back(tmpY[idx_plk]);
          maxCount++; 
          plank = false;
          cnt_plk = 0;
          continue;
        }
        //low plank: pdf is increasing after the plank. No peak here!
        if(tmpY[i] < tmpY[i + 1]){
          plank = false; //end
          cnt_plk = 0;
        } else { //Plank is getting larger
          cnt_plk++;
        }
      }
    }
    // last element: check if plank is extended
    if (tmpY[lastNonNanInd - 1] <= tmpY[lastNonNanInd]) {
      Int32 idx;
      if(plank){
        cnt_plk++;
        idx = lastNonNanInd - round(cnt_plk/2);
      }else{ 
        idx = lastNonNanInd;
      }

      maxX.push_back(tmpX[idx]);
      maxY.push_back(tmpY[idx]);
      maxCount++;
    }
  return  true;
}

Bool CExtremum::FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, TFloat64List& maxX, TFloat64List& maxY, Float64 SignSearch) const 
{
  if (n == 0)
    return false;

  vector < Float64 > tmpX(n);
  vector < Float64 > tmpY(n);

  for (UInt32 t = 0; t < n; t++) {
    tmpX[t] = xAxis[t];
    tmpY[t] = SignSearch * yAxis[t];
  }

  Int32 maxCount = 0;

    // find first and last non Nan element
    Int32 firstNonNanInd = 0;
    for (Int32 iFirst = 0; iFirst < n - 1; iFirst++) {
      if (!std::isnan((double) tmpY[iFirst])) {
        firstNonNanInd = iFirst;
        break;
      }
    }
    Int32 lastNonNanInd = n - 1;
    for (Int32 iLast = n - 1; iLast > 0; iLast--) {
      if (!std::isnan(tmpY[iLast])) {
        lastNonNanInd = iLast; 
        break;
      }
    }
    
    bool plank = false; Int32 cnt_plk = 0; 
    // First element
    if (tmpY[firstNonNanInd] > tmpY[firstNonNanInd + 1]) {
      maxX.push_back(tmpX[firstNonNanInd]);
      maxY.push_back(tmpY[firstNonNanInd]);
      maxCount++;
    }else{
      if(tmpY[firstNonNanInd] == tmpY[firstNonNanInd + 1]) {
        plank = true;
        cnt_plk++;
      }
    }

    for (Int32 i = firstNonNanInd + 1; i < lastNonNanInd; i++) {
     if ((tmpY[i] > tmpY[i - 1]) && (tmpY[i] > tmpY[i + 1])) {
        maxX.push_back(tmpX[i]);
        maxY.push_back(tmpY[i]);
        maxCount++;
        continue;
      }
      if((tmpY[i] > tmpY[i - 1]) && (tmpY[i] == tmpY[i + 1])){
        //plank start point
        plank = true;
        cnt_plk++;
        continue;
      }
      if((tmpY[i] == tmpY[i - 1])){ 
        //high plank: signal is decreasing after the plank. Peak identified
        if (tmpY[i] > tmpY[i + 1] && plank){ //check if we already identified a high plank
          cnt_plk++;
          Int32 idx_plk = i - round(cnt_plk/2);
          //plank end point
          maxX.push_back(tmpX[idx_plk]);
          maxY.push_back(tmpY[idx_plk]);
          maxCount++; 
          plank = false;
          cnt_plk = 0;
          continue;
        }
        //low plank: pdf is increasing after the plank. No peak here!
        if(tmpY[i] < tmpY[i + 1]){
          plank = false; //end
          cnt_plk = 0;
        } else { //Plank is getting larger
          cnt_plk++;
        }
      }
    }
    // last element: check if plank is extended
    if (tmpY[lastNonNanInd - 1] <= tmpY[lastNonNanInd]) {
      Int32 idx;
      if(plank){
        cnt_plk++;
        idx = lastNonNanInd - round(cnt_plk/2);
      }else{ 
        idx = lastNonNanInd;
      }

      maxX.push_back(tmpX[idx]);
      maxY.push_back(tmpY[idx]);
      maxCount++;
    }
  return  true;
}
/**
 * Brief: Reduce number of peaks based on maxCount passed from param.json if present
*/
Bool CExtremum::Truncate( TFloat64List& maxX, TFloat64List& maxY, TPointList& maxPoint) const {
  UInt32 n = maxX.size();
  n = std::min(m_MaxPeakCount, n);
  for (Int32 j = 0; j < n; j++) { 
    maxPoint.push_back(SPoint(maxX[m_sortedIndexes[j]], m_SignSearch * maxY[m_sortedIndexes[j]]) );
  }   
  return true;
}

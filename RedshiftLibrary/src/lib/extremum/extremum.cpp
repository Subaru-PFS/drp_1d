#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>

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
    m_Radius(0.005),
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
CExtremum::CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount, Float64 radius, Bool invertForMinSearch) :
    m_MaxPeakCount( maxPeakCount ),
    m_XRange( xRange ),
    m_SignSearch( 1.0 ),
    m_Radius(radius),
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
 * Empty destructor.
 */
CExtremum::~CExtremum()
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

void CExtremum::DeactivateSlidingWindow()
{
    m_slidingWindowactive = false;
}
/**
 * Sets m_SignSearch to val.
 */
void CExtremum::SetSignSearch( Float64 val )
{
    m_SignSearch = val;
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
    
    //Find does 4 things: findAllPeaks; Cut_Threshold; PeakRefinement using sliding windows; Truncation based on allowed peak number
    //Didier proposes to set extremaCount=100 (hardcoded value) when called from linemodelsolve.
    //And that extremaCount, as defined in param.json, is taken into consideration elsewhere, at the end of the processing
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
    m_sortedIndexes.resize(maxX.size());
    iota(m_sortedIndexes.begin(), m_sortedIndexes.end(), 0);
    sort(m_sortedIndexes.begin(), m_sortedIndexes.end(),
       [&maxY](Float64 i1, Float64 i2) {return maxY[i1] > maxY[i2];});

    //m_boolTracking correspond to peaks of sortedIndexes
    Int32 keepMinN = 2;
    if(m_meritCut>0.0 && maxX.size()>keepMinN){ 
        Bool v = Cut_Threshold(maxX, maxY, keepMinN);
    }else{
      TFloat64List tmpX(maxX), tmpY(maxY);
      maxX.clear(); maxY.clear();
      for (Int32 i = 0; i<m_sortedIndexes.size(); i++){
        maxX.push_back(tmpX[m_sortedIndexes[i]]);
        maxY.push_back(tmpY[m_sortedIndexes[i]]);
      }
    }
    //refine using sliding windows: avoiding duplicate candidates when possible. 
    Bool b;
    if(m_slidingWindowactive)
      b = FilterOutNeighboringPeaks_2(maxX, maxY, keepMinN, maxPoint);
    //verify that peaks are well separated by at least secondpassradius
    Bool verified = verifyPeakSeparation(maxPoint);
    return true;
}
//passing by copy to not risk changing maxX order without changing maxY
//considers that maxX are sorted 
Bool CExtremum::verifyPeakSeparation( TFloat64List& maxX) const{
  Bool verified = true;
 
  std::sort(maxX.begin(), maxX.end());
  for(Int32 i = 0; i< maxX.size() -1; i++)  {
    Float64 overlap;
    overlap = std::max( maxX[i+1] - m_Radius*(1+maxX[i+1]), 0.0 ) - (maxX[i] + m_Radius*(1+maxX[i]));
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
    Float64 meritDiff =  tmpY[m_sortedIndexes[0]]- tmpY[m_sortedIndexes[i]];
    if( meritDiff > m_meritCut && i>=keepMinN){
      //truncate m_sortedIndexes
      m_sortedIndexes.resize(i-1);
      break; 
    }
    else{
        maxX.push_back(tmpX[m_sortedIndexes[i]]);
        maxY.push_back(tmpY[m_sortedIndexes[i]]);
    }
  }
  return true;
}

/**
 * Brief: extending find peak function to inlude distance, threshold in peak selection
 * @xAxis and @yAxis represents local maxima from 
 * Eliminate candidates with a very low Y value, comparing to others?
 * Eliminate close candidates: close, i.e., in a window of 0.005(1+z) around a z candidate of strong value???
 * //TODO: Need to implement the keepmin: probably no need for this!!
*/
Bool CExtremum::FilterOutNeighboringPeaks(TFloat64List& maxX, TFloat64List& maxY, UInt32 keepmin)const
{
  if(maxX.size() <= keepmin)
    return true;

  vector < Float64 >  tmpX;
  vector < Float64 >  tmpY;
  //copy maxX/Y into temperory arrays
  for (Int32 i = 0; i<maxX.size(); i++){
    tmpX.push_back(maxX[i]);
    tmpY.push_back(maxY[i]);
  }

  //clear vectors to fill later
  maxX.clear(); maxY.clear();

  //starting from the first element
  //minimal distance between two close extrema should be >= secondpass_range used to find extended redshifts
  Float64 wind_high = tmpX[0] + 2*m_Radius*(1+tmpX[0]);
  vector<Int32> idxList; 
  Int32 i = 0, imax = -1; 
  Float64 maxPDF = -INFINITY;
  Bool fullwdw = false;

  while(i<tmpX.size()){
    //if element belongs to the selected half-window
    if(tmpX[i]<=wind_high){
      idxList.push_back(i);
      if(tmpY[i]>maxPDF){//keep track of maxPDF element
        maxPDF = tmpY[i];
        imax = i;
      }
      //case where the last element belongs to the current window
      //we need to push 
      if(i==tmpX.size()-1 && imax>-1){         
        if(!maxX.size()){
              maxX.push_back(tmpX[imax]);
              maxY.push_back(tmpY[imax]);
        }else{
              Float64 overlap;
              overlap = (tmpX[imax] - m_Radius*(1+tmpX[imax]) ) - (maxX.back() + m_Radius*(1+maxX.back()));
              if(overlap <0) {
                if(tmpY[imax]>maxY.back()){
                  //push the best PDF as a best local candidate before calculating a new window
                  maxX[maxX.size()-1] = tmpX[imax];
                  maxY[maxY.size()-1] = tmpY[imax];
                }//else //nothing to do cause best cand is already saved in maxX
              }else{
                maxX.push_back(tmpX[imax]);
                maxY.push_back(tmpY[imax]);
              }
        }
        break;
      }else{
        i++;
        continue;
      }
    }

    if(tmpX[i]>wind_high){
     
      if(fullwdw){
        //before pushing check if the newly identified best is very close to the previos bestZ
        //if the two following best candid do not respect the strict distance, keep only one of them
        Float64 overlap;
        if(!maxX.size()){
          maxX.push_back(tmpX[imax]);
          maxY.push_back(tmpY[imax]);
        }else{
          overlap = (tmpX[imax] - m_Radius*(1+tmpX[imax]) ) - (maxX.back() + m_Radius*(1+maxX.back()));
          if(overlap <0) {
            if(tmpY[imax]>maxY.back()){
              //push the best PDF as a best local candidate before calculating a new window
              maxX[maxX.size()-1] = tmpX[imax];
              maxY[maxY.size()-1] = tmpY[imax];
            }//else //nothing to do cause best cand is already saved in maxX
          }else{
            maxX.push_back(tmpX[imax]);
            maxY.push_back(tmpY[imax]);
          }
        }
        //reinitialize all
        fullwdw = false; 
        imax = -1; maxPDF = -INFINITY;
        //calculate a new window based on the current i
        wind_high = tmpX[i] + 2*m_Radius*(1+tmpX[i]);
      }else {
        if(imax == -1)
          imax = i; 
        wind_high = tmpX[imax] + 2*m_Radius*(1+tmpX[imax]);//add half wdw to reach a full wdw
        fullwdw = true; 
        continue;
      }
      
    }
  }

  return true;
}

//todo use the tracking mecanism
Bool CExtremum::FilterOutNeighboringPeaks_2(TFloat64List& maxX, TFloat64List& maxY, UInt32 keepmin, TPointList& maxPoint)const
{
  if(maxX.size()<= keepmin || maxX.size()<= m_MaxPeakCount){
    for(Int32 i = 0; i<maxX.size(); i++){
      maxPoint.push_back(SPoint(maxX[i], m_SignSearch * maxY[i]) );
    }
    return true;
  }

  //Int32 n = std::count(m_boolTracking.begin(), m_boolTracking.end(), true);
  TBoolList peakTracking(maxX.size(), true);
  Float64 wind_high, wind_low;
  Int32 i = 0; 
  /**
   * start from the highest proba candidate 
   * construct a wdw to its right and left
   * check presence of candidates in this wdw
   * eliminate neighboring candidates  
   * update tmpX and tmpY
  */

  while(i<maxX.size()){
    //skipped peak
    if(!peakTracking[i]){
      i++;
      continue;
    }
    //compute wdw_high with respect to a slightly higher redshift(Z+), i.e., z+ = Z0 + delta*(1 + z+)
    wind_high = (maxX[i] + m_Radius)/ (1 - m_Radius); //tmpX + m_Radius*(1+tmpX);
    wind_low = std::max( maxX[i] - m_Radius*(1+maxX[i]), 0.0);
    for(Int32 j = i+1; j< maxX.size(); j++){
      //if j is eliminated, skip it and move to the next index
      if(!peakTracking[j])
        continue;
      Float64 overlap = 0.0;
      if(maxX[j] > maxX[i]){
        overlap = -wind_high + std::max( maxX[j] - m_Radius*(1+maxX[j]), 0.0 );
        if(overlap < 0.0){
          peakTracking[j] = false;
        }
        continue; 
      }
      if(maxX[j] < maxX[i]){
        overlap = wind_low - (maxX[j] + m_Radius*(1+maxX[j]));
        if(overlap < 0.0){
          peakTracking[j] = false;
        }
        continue;  
      }
    }
    maxPoint.push_back(SPoint(maxX[i], m_SignSearch * maxY[i]) );
    //if at index i, we have already identified our MaxPeakCount as valid Peaks, break the loop
    //no need to continue if we have reached the required MaxPeakCount
    if(maxPoint.size() == m_MaxPeakCount)
      break;
    i++;
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
 * considering that maxY values are sorted in decreasing order
*/
Bool CExtremum::Truncate( TFloat64List& maxX, TFloat64List& maxY, TPointList& maxPoint) const {
  UInt32 n = maxX.size();
  n = std::min(m_MaxPeakCount, n);
  for (Int32 j = 0; j < n; j++) { 
    maxPoint.push_back(SPoint(maxX[j], m_SignSearch * maxY[j]) );
  }   
  return true;
}

Bool CExtremum::Truncate( TFloat64List& maxX, TFloat64List& maxY, Int32 maxCount, TPointList& maxPoint) const {
  
  Int32 n = maxX.size();
  maxCount = std::min(maxCount, n);

  TInt32List sortedIndexes(n);
  iota(sortedIndexes.begin(), sortedIndexes.end(), 0);
  sort(sortedIndexes.begin(), sortedIndexes.end(),
       [&maxY](Float64 i1, Float64 i2) {return maxY[i1] > maxY[i2];});
  // save candidates as SPoints
  for (Int32 j = 0; j < maxCount; j++) {
     Int32 jidx = sortedIndexes[j];
    if( !std::isnan( maxY[jidx] ) ){
      maxPoint.push_back(SPoint(maxX[jidx], m_SignSearch * maxY[jidx]) );
    }
  }
  return true;
}
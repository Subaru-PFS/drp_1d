#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>

#include <cmath>
#include <float.h>
#include <iostream>

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
void CExtremum::SetMeritCut( UInt32 n )
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
    Bool method_min = FindAllPeaks( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, minX, minY, -1*m_SignSearch); 
    for(Int32 i = 0; i <minX.size(); i++){
      minY[i] = -1*minY[i];
    }
    //Calculate prominence and remove "low" prominence peaks
    //this is ALSO useful for eliminating neighboring peaks in some cases
    Cut_Prominence_Merit(maxX, maxY, minX, minY);

    //refine using sliding windows: aiming at avoiding duplicate candidates when possible. 
    FilterOutNeighboringPeaks(maxX, maxY, keepMinN);//keep at least keepMinN candidates

    //Cut_Threshold is optional
  /*  if(m_meritCut>0.0){ 
      Bool v = Cut_Threshold(maxX, maxY);
    }*/
    //truncate: reduces size of candidate list and also prepares the maxPoint List
    Truncate(maxX, maxY, m_MaxPeakCount, maxPoint);

    return true;
}
/**
 * prominence: vertical distance between a summit and the key col, i.e., the closest (horizontal) minima
 * joint prominence_merit cut
 * Method: 
 * 1. identify for consecutive peaks the key col.
 * 2. calculate prominence for the current peak
 * 3. if prominence == peak.y, then peak is a very high peak
 * 4. Remove low prominence peaks : TO decide on "low" value
 * 
 * To identify key col for each peak:
 * 1. Extend to the right and to the left of the peak until reaching higher peaks (ya_current < ya_right_after) //skip peaks of same value and continue extending
 * 2. Find the highest minima within the range identified in (1), considered as the key col.
*/
Bool CExtremum::Cut_Prominence_Merit( vector <Float64>& maxX, vector <Float64>& maxY, vector <Float64>& minX, vector <Float64>& minY) const{
  //find highest peak:
  Float64 maxV = *std::max_element(maxY.begin(), maxY.end()); //max of maxY
  Float64 minV = *std::min_element(minY.begin(), minY.end()); //min of minY
  Float64 ref_prominence = maxV - minV; //used to normalize obtained prominence TODO: requires discussion
  Float64 prominence_thresh = 0.1; //heuristic value, TODO: requires discussion (maybe passed in param.json)
  vector <Float64> prominence(maxX.size()), tmpX, tmpY; 

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
    //key_col is the highest minima between the lowest minimum to the right of the peak abd to the left of the peak
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
    prominence[i] = (maxY[i] - std::max(key_coly_l, key_coly_r))/ ref_prominence;
    //keep peaks whose height is almost equal to their prominence
    if(prominence[i] > prominence_thresh || (m_meritCut &&(maxV - maxY[i] < m_meritCut)) ){ //heuristic value
      tmpX.push_back(maxX[i]);
      tmpY.push_back(maxY[i]);
    }
  }
  maxX.clear(); maxY.clear();
  for(Int32 i = 0; i<tmpX.size(); i++){
    maxX.push_back(tmpX[i]);
    maxY.push_back(tmpY[i]);
  }
  return true;
}
/**
 * \Brief: removes extrema based on meritCut.
 * First: Sort based on Y values
 * Second: Remove non-relevant candidates 
 * Third: Sort based on X values to return candidates following their initial order
 * It returns the new list of extrema
*/
Bool CExtremum::Cut_Threshold( vector <Float64>& maxX, vector <Float64>& maxY, Int32 keepMinN) const{
  Int32 n = maxX.size();
  //create pairs of X and Y
  vector<pair<Float64,Float64> > vp, vp_;
  vp.reserve(n);
  for (Int32 i = 0 ; i < n ; i++) {
    vp.push_back(make_pair(maxY[i], maxX[i]));
  }
  std::sort(vp.rbegin(), vp.rend()); //sort descending order
  vp_.push_back(make_pair(vp[0].second, vp[0].first)); //save best one
  
  maxX.clear(); maxY.clear();
  for(Int32 i = 1; i<n-1; i++){
    Float64 meritDiff = vp[0].first - vp[i].first;
    if( meritDiff > m_meritCut && i>=keepMinN){
      break; //no need to continue iterating since vp is sorted!
    }else{
      vp_.push_back(make_pair(vp[i].second, vp[i].first));
    }
  }
  std::sort(vp_.rbegin(), vp_.rend()); //sort based on maxX values
  Int32 s = vp_.size();
  for(Int32 i = 1; i<s+1; i++){
    maxX.push_back(vp_[s-i].first);
    maxY.push_back(vp_[s-i].second);
  }
  return true;
}

/**
 * Brief: extending find peak function to inlude distance, threshold in peak selection
 * @xAxis and @yAxis represents local maxima from 
 * Eliminate candidates with a very low Y value, comparing to others?
 * Eliminate close candidates: close, i.e., in a window of 0.005(1+z) around a z candidate of strong value???
 * //TODO: Need to implement the keepmin
*/
Bool CExtremum::FilterOutNeighboringPeaks(vector <Float64>& maxX, vector <Float64>& maxY, UInt32 keepmin)const
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

  //starting from first element
  Float64 wind_high = tmpX[0] + m_Radius*(1+tmpX[0]);
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
      //i++; 
      //case where the last element belongs to the current window
      //we need to push 
      if(i==tmpX.size()-1 && imax>-1){
          maxX.push_back(tmpX[imax]);
          maxY.push_back(tmpY[imax]);
          break;
      }else{
        i++;
        continue;
      }
    }

    if(tmpX[i]>wind_high){
     
      if(fullwdw){
        //push the best PDF as a best local candidate before calculating a new window
        maxX.push_back(tmpX[imax]);
        maxY.push_back(tmpY[imax]);
        //reinitialize all
        fullwdw = false; 
        imax = -1; maxPDF = -INFINITY;
        //calculate a new window based on the current i
        wind_high = tmpX[i] + m_Radius*(1+tmpX[i]);
      }else {
        //Check if last elet. if yes, push it to maxX/Y
        if(i==tmpX.size()-1 && imax>-1){
          if(imax>-1){
            maxX.push_back(tmpX[imax]);
            maxY.push_back(tmpY[imax]);
          }
          //Maybe there is here a need to compare with the last identified
          break; //end of loop
        }else{
          //calculate a new window based on the imax we found till now
          //keep maxPDF set to compare others with it
          wind_high = tmpX[imax] + m_Radius*(1+tmpX[imax]);
          fullwdw = true; 
          continue;
        }
      }
      
    }
  }

  return true;
}

/**
 * Brief: Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 * Extended to include distance between peaks and a threshold to eliminate peaks
 */
Bool CExtremum::FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, vector <Float64>& maxX, vector <Float64>& maxY) const {
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

Bool CExtremum::FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, vector <Float64>& maxX, vector <Float64>& maxY, Float64 SignSearch) const {
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
 * Reduction happens based on PDF values (we eliminate peaks with the smallest values)
*/
Bool CExtremum::Truncate( vector <Float64>& maxX, vector <Float64>& maxY, Int32 maxCount, TPointList& maxPoint) const {
  
  Int32 n = maxX.size();
  maxCount = std::min(maxCount, n);

  CQuickSort<Float64> sort;
  vector<Int32> sortedIndexes(n);
  sort.SortIndexes( maxY.data(), sortedIndexes.data(), sortedIndexes.size() );
 
  // save candidates as SPoints
  for (Int32 j = 0; j < maxCount; j++) {
     Int32 jidx = sortedIndexes[(sortedIndexes.size()- j -1)];
    if( !std::isnan( maxY[jidx] ) ){
      maxPoint.push_back(SPoint(maxX[jidx], m_SignSearch * maxY[jidx]) );
    }
  }
  return true;
}
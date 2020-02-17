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
    vector < Float64 >  maxX;
    vector < Float64 >  maxY;
    Bool method = FindAllPeaks( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, maxX, maxY );
 
    Int32 keepMinN = 2;

    //Cut_Threshold is optional
    if(m_meritCut>0.0){
      Bool v = Cut_Threshold(maxX, maxY, keepMinN);
    }

    //refine using sliding windows: aiming at avoiding duplicate candidates when possible. 
    FilterOutNeighboringPeaks(maxX, maxY, keepMinN);//keep at least keepMinN candidates

    //truncate: reduces size of candidate list and also prepares the maxPoint List
    Truncate(maxX, maxY, m_MaxPeakCount, maxPoint);

    return true;
}
/**
 * \Brief: removes extrema based on meritCut.
 * First: sort based on Y values
 * Second: remove non-relevant candidates 
 * Third: sort based on X values to return candidates following their initial order
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
  vector < Float64 >  tmpX;
  vector < Float64 >  tmpY;
  //copy maxX/Y into temperory arrays
  for (Int32 i = 0; i<maxX.size(); i++){
    tmpX.push_back(maxX[i]);
    tmpY.push_back(maxY[i]);
  }

  //clear vectors to fill later
  maxX.clear();
  maxY.clear();

  Float64 wind_low = tmpX[1] - m_Radius*(1+tmpX[1]);
  Float64 wind_high = tmpX[1] + m_Radius*(1+tmpX[1]);
  vector<Int32> idxList, missedList; 
  Int32 i = 0, imax = -1; 
  Float64 maxPDF = -INFINITY;

  while(i<tmpX.size()){
    //if element belongs to the selected window
    if(tmpX[i]>wind_low && tmpX[i]<=wind_high){
      idxList.push_back(i);
      if(tmpY[i]>maxPDF){//keep track of maxPDF element
        maxPDF = tmpY[i];
        imax = i;
      }
      i++;
      continue;
    }
    //we should enter below only after updating window size
    if(tmpX[i]<wind_low){
      //TODO:test if its PDF is relevant maybe, ??
      //otherwise add it by default (is dangerous!)
      maxX.push_back(tmpX[i]);
      maxY.push_back(tmpY[i]);

      missedList.push_back(i);
      i++;
    }

    if(tmpX[i]>wind_high){
      //push the best PDF as a best local candidate before calculating a new window
      maxX.push_back(tmpX[imax]);
      maxY.push_back(tmpY[imax]);
      
      if(i == tmpX.size()-1 ){ //PUSH ANYWAY!! 
        maxX.push_back(tmpX[i]);
        maxY.push_back(tmpY[i]);
        idxList.push_back(i);//update listed indexes for later checks
        break; //leave the while lopp
      }
      //reinitialize imax and maxPDF
      imax = -1; maxPDF = -INFINITY;
      //calculate a new window based on two candidates ahead of the current one
      Int32 startIdx;
      if(i+1 < tmpX.size()){
        startIdx = i+1;
      }else{
        if(i < tmpX.size()){
          startIdx = i;
        }
      }
      wind_low = tmpX[startIdx] - m_Radius*(1+tmpX[startIdx]);
      wind_high = tmpX[startIdx] + m_Radius*(1+tmpX[startIdx]);
      //dont increment so that we can check the value of tmpX with the newly calculated window
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
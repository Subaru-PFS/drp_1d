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
    m_MaxPeakCount( 5 ),
    m_RefreshCount( 1 ),
    m_XRange( 0.0, 0.0 ),
    m_SignSearch( 1.0 ),
    m_Radius(0.005)
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
CExtremum::CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount, Bool invertForMinSearch, UInt32 refreshCount, Float64 radius ) :
    m_MaxPeakCount( maxPeakCount ),
    m_RefreshCount( refreshCount ),
    m_XRange( xRange ),
    m_SignSearch( 1.0 ),
    m_Radius(radius)
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
 * Sets m_RefreshCount to n.
 */
void CExtremum::SetRefreshCount( UInt32 n )
{
    m_RefreshCount = n;
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
    TPointList maxPoint_ ;
    Bool method1 = InternalFind( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, maxPoint_ );
         
    
    Bool method2 = InternalFind_refact_ext( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, maxPoint );
    return method2;
}
/**
 * Brief: removes extrema based on meritCut.
 * It returns the new list of extrema
*/
/*Bool CExtremum::Cut_Threshold(const TFloat64List& xAxis, const TFloat64List& yAxis, UInt32 meritCut, Int32 keepMinN, TPointList& maxPoint ) const{

}*/

Bool CExtremum::Cut_Threshold( TPointList& maxPoint, UInt32 meritCut, Float64 bestPDF, Int32 keepMinN) const{
  Int32 n = maxPoint.size();
  Int32 iExtremumFinalList = 0;
  for (Int32 i = 0; i < n; i++){
    Float64 meritDiff = bestPDF - maxPoint[iExtremumFinalList].Y;
    if( meritDiff > meritCut && i>=keepMinN){
        Log.LogInfo("  Extremum: Candidates selection by proba cut: removing i=%d, final_i=%d, e.X=%f, e.Y=%e",
                            i,
                            iExtremumFinalList,
                            maxPoint[iExtremumFinalList].X,
                            maxPoint[iExtremumFinalList].Y);
        maxPoint.erase(maxPoint.begin() + iExtremumFinalList);
        }else{
          iExtremumFinalList++;
        }
  }
  return true;
}

/**
 * Brief: extending find peak function to inlude distance, threshold in peak selection
 * @xAxis and @yAxis represents local maxima from 
 * Eliminate candidates with a very low Y value, comparing to others?
 * Eliminate close candidates: close, i.e., in a window of 0.005(1+z) around a z candidate of strong value???
 * //TODO: should rethink about the order for applying  these conditions
*/
Bool CExtremum::FindPeaks_extended(const TFloat64List& xAxis, const TFloat64List& yAxis, UInt32 n, TPointList& maxPoint ) const
{
  if (n == 0)
    return false;
  vector < Float64 >  tmpX;
  vector < Float64 >  tmpY;
  vector < Float64 >  maxX;
  vector < Float64 >  maxY;


  CQuickSort<Float64> sort_;
  vector<Int32> sortedIndexes_(n);
  sort_.SortIndexes( yAxis.data(), sortedIndexes_.data(), sortedIndexes_.size() );
  //best value corresponds to the last in sortedIndexes
  Float64 best_pdf = yAxis[sortedIndexes_[n-1]];
  //do a first selection based on maxY value relative to the best tmpY
  Int32 cnt = 0;  

//below doesnt do anything given that meritCut is set to max value
  Float64 /*pdfDiff = DBL_MAX,*/ meritCut = DBL_MAX;//30;
  for (Int32 i = 0; i< n; i++){
    if(i == sortedIndexes_[n-1]){
      tmpX.push_back(xAxis[i]);
      tmpY.push_back(yAxis[i]);
      continue;
    }
    /*
    //looking for the smallest diff, i.e., 10^7..too big
    if( pdfDiff>(best_pdf - yAxis[i])){
      pdfDiff = best_pdf - yAxis[i];
    }*/
    if((best_pdf - yAxis[i]) < meritCut){
      tmpX.push_back(xAxis[i]);
      tmpY.push_back(yAxis[i]);
    }
  }

  //Mira: TODO: should pass radius as an optional argument probably
  //Float64 radius = 0.005;//similar to when we extend around zcand
  if(tmpX.size() == 0){
     return false;
  }
 
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
 
//we need to check if we didnt skip any candidate because of window borders
  Int32 s = missedList.size();

  Int32 nbPeaks = m_MaxPeakCount;
  if( maxX.size()<m_MaxPeakCount ){
    nbPeaks = maxX.size();
  }
  maxPoint.resize( nbPeaks );
 // Extract best results
  CQuickSort<Float64> sort;
  vector<Int32> sortedIndexes( maxX.size() );
  sort.SortIndexes( maxY.data(), sortedIndexes.data(), sortedIndexes.size() );

  // save candidates as SPoints
  for (Int32 j = 0; j < nbPeaks; j++) {
     Int32 jidx = sortedIndexes[(sortedIndexes.size()- j -1)];
    if( !std::isnan( maxY[jidx] ) ){
      maxPoint[j] = SPoint(maxX[jidx], m_SignSearch * maxY[jidx]);
    }
  }
  return true;
}

/**
 * Brief: Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 * Extended to include distance between peaks and a threshold to eliminate peaks
 */
Bool CExtremum::InternalFind_refact_ext(const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint) const {
  if (n == 0)
    return false;

  vector < Float64 > maxX;
  vector < Float64 > maxY;
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
        if (tmpY[i] > tmpY[i + 1]){
          cnt_plk++;
          //TODO: talk to Didier: which element of the plank should be kept as Zcand, last elt or plank center?
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

    //compute average and variance of yaxis: not used
    Float64 avg = 0.0, var = 0.0;
    /*for (Int32 i = firstNonNanInd; i < lastNonNanInd + 1; i++) {
      avg = avg + yAxis[i]; 
    }
    avg = avg/(lastNonNanInd - firstNonNanInd + 1);

    for (Int32 i = firstNonNanInd; i < lastNonNanInd + 1; i++) {
      var += (yAxis[i] - avg)*(yAxis[i] - avg); 
    }
    var /= (lastNonNanInd - firstNonNanInd + 1);*/
    
  //select candidates
  return  FindPeaks_extended( maxX, maxY, maxCount, maxPoint );
}
/**
 * Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 */
Bool CExtremum::InternalFind( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    if( n == 0 )
        return false;

    //Method 1, use only 1 extremum
    /*
    maxPoint.resize( 1 );

    Float64 max = DBL_MIN ;
    Int32 maxIndex = 0;
    for( Int32 i=0; i<n; i++ )
    {
    	if( yAxis[i] > max ) {
    		max = yAxis[i];
    		maxIndex = i;
    	}
    }

    maxPoint[0].X = xAxis[maxIndex];
    maxPoint[0].Y = yAxis[maxIndex];

    return true;
    //*/

    vector<Float64> maxX( n );
    vector<Float64> maxY( n );

    // Tmp array can be considered as the "input" of each iteration.
    vector<Float64> tmpX( n );
    vector<Float64> tmpY( n );
    UInt32 tmpSize = n;
    for( UInt32 t=0; t<n; t++ )
    {
        tmpX[t] = xAxis[t];
        tmpY[t] = m_SignSearch*yAxis[t];
    }


    for( UInt32 count=0; count<m_RefreshCount; count++ )
    {
        Int32 maxCount = 0;

        // find first and last non Nan element
        Int32 firstNonNanInd = 0;
        for( Int32 iFirst=0; iFirst<tmpSize-1; iFirst++ )
	  {
            if( !std::isnan( (double) tmpY[iFirst] ) )
	      {
                firstNonNanInd = iFirst;
                break;
	      }
	  }
        Int32 lastNonNanInd = tmpSize-1;
        for( Int32 iLast=tmpSize-1; iLast>0; iLast-- )
	  {
            if( !std::isnan( tmpY[iLast] ) )
	      {
                lastNonNanInd = iLast;
                break;
	      }
	  }

        // First element
        if( tmpY[firstNonNanInd] > tmpY[firstNonNanInd+1] )
	  {
            maxX[maxCount] = tmpX[firstNonNanInd];
            maxY[maxCount] = tmpY[firstNonNanInd];
            maxCount++;
	  }

        for( Int32 i=firstNonNanInd+1; i<lastNonNanInd; i++ )
	  {
            if( ( tmpY[i]>tmpY[i-1] ) && ( tmpY[i]>tmpY[i+1] ) )
	      {
                maxX[maxCount] = tmpX[i];
                maxY[maxCount] = tmpY[i];
                maxCount++;
	      }
	  }

        // last elements
        if( tmpY[lastNonNanInd-1]<tmpY[lastNonNanInd] )
	  {
            maxX[maxCount] = tmpX[lastNonNanInd];
            maxY[maxCount] = tmpY[lastNonNanInd];
            maxCount++;
	  }

        //tmpX = vector<Float64>( maxCount );
        //tmpY = vector<Float64>( maxCount );

        tmpSize = maxCount;
        // Prepare for next iteration by storing every maximum found in this iteration in the tmp array used as input for the next iteration
        for( Int32 t=0; t<maxCount; t++ )
	  {
            tmpX[t]=maxX[t];
            tmpY[t]=maxY[t];
	  }

        /*//debug:
        // save median and xmad,  flux data
        FILE* f = fopen( "extremum_dbg.txt", "w+" );
        for( Int32 t=0;t<maxCount;t++)
        {
            fprintf( f, "%d %f %f\n", t, tmpX[t], tmpY[t]);
        }
        fclose( f );
        //*/

        if( maxCount == 0 )
	  {
            break;
	  }

        if( tmpSize <= PEAKS_SMOOTH_LIMIT )
	  {
            break;
	  }
    }

    Int32 nbPeaks = m_MaxPeakCount;
    if( tmpSize<m_MaxPeakCount )
      {
        nbPeaks = tmpSize;
      }

    // Extract best results
    CQuickSort<Float64> sort;

    vector<Int32> sortedIndexes( tmpSize );
    sort.SortIndexes( tmpY.data(), sortedIndexes.data(), sortedIndexes.size() );

    maxPoint.resize( nbPeaks );
    Int32 k = 0;

    for( Int32 i=0; i<nbPeaks; i++ )
      {
        Int32 j = sortedIndexes[(sortedIndexes.size()-1)-i];
        if( !std::isnan( tmpY[j] ) )
	  {
            maxPoint[k++] = SPoint( tmpX[j], m_SignSearch*tmpY[j] );
	  }
      }

    return true;
}


/**
 * Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 */
Bool CExtremum::InternalFind2( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    if( n == 0 )
        return false;

    vector<Float64> tmpX( n );
    vector<Float64> tmpY( n );
    Int32 tmpSize = n;
    for( UInt32 t=0; t<n; t++ )
      {
        if( !std::isnan( yAxis[t] ) )
	  {
            tmpX[t] = xAxis[t];
            tmpY[t] = yAxis[t];
	  }
      }

    // Extract best results
    CQuickSort<Float64> sort;

    vector<Int32> sortedIndexes( tmpSize );
    sort.SortIndexes( tmpY.data(), sortedIndexes.data(), tmpSize );

    maxPoint.resize( m_MaxPeakCount );
    Int32 k = 0;

    for( UInt32 i=0; i<m_MaxPeakCount; i++ )
      {
        UInt32 j = sortedIndexes[(sortedIndexes.size()-1)-i];

        maxPoint[k++] = SPoint( tmpX[j], tmpY[j] );
      }

    return true;
}

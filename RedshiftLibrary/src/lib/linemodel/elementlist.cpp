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
#include "RedshiftLibrary/linemodel/elementlist.h"

using namespace NSEpic;


CLineModelElementList::CLineModelElementList(const CRayCatalog::TRayVector &restRayList):
  m_RestRayList(restRayList)
{}

/**
 * \brief Returns the number of m_Elements that fail IsOutsideLambdaRange().
 **/
Int32 CLineModelElementList::GetModelValidElementsNDdl()
{
    Int32 nddl = 0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }

        nddl++;
    }
    return nddl;
}

/**
 * \brief Returns the number of elements that have only subelements with non-positive amplitude.
 **/
Int32 CLineModelElementList::GetModelNonZeroElementsNDdl()
{
    Int32 nddl = 0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        bool isAllZero=true;
        for(Int32 ie=0; ie<m_Elements[iElts]->GetSize(); ie++){
            if(m_Elements[iElts]->GetFittedAmplitude(ie) > 0.0){
                isAllZero=false;
            }
        }

        if(isAllZero==false){
            nddl++;
        }
    }
    return nddl;
}

/**
 * \brief Returns the list of indexes of elements that fail IsOutsideLambdaRange.
 **/
std::vector<UInt32> CLineModelElementList::GetModelValidElementsIndexes()
{
    std::vector<UInt32> nonZeroIndexes;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        if(IsElementIndexInDisabledList(iElts))
        {
            continue;
        }

        nonZeroIndexes.push_back(iElts);
    }
    return nonZeroIndexes;
}

bool CLineModelElementList::IsElementIndexInDisabledList(Int32 index)
{
    for( UInt32 i=0; i<m_elementsDisabledIndexes.size(); i++ )
    {
        if( m_elementsDisabledIndexes[i]== index){
            return true;
        }
    }
    return false;
}

/**
 * @brief CLineModelElementList::SetElementIndexesDisabledAuto
 * Disables all the elements that have all sub-elements (lines) amplitudes equal to zero
 */
void CLineModelElementList::SetElementIndexesDisabledAuto()
{
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        bool isAllZero=true;
        for(Int32 ie=0; ie<m_Elements[iElts]->GetSize(); ie++){
            if(m_Elements[iElts]->GetFittedAmplitude(ie) > 0.0){
                isAllZero=false;
            }
        }

        if(isAllZero==true){
            m_elementsDisabledIndexes.push_back(iElts);
        }
    }
}

void CLineModelElementList::ResetElementIndexesDisabled()
{
    m_elementsDisabledIndexes.clear();
}

/**
 * \brief Returns the list of groups, with each group being a set of line indexes with the velcocity to be jointly
 * TEMPORARY-DEV: return all the indexes individually as  agroup
**/
std::vector<std::vector<Int32>> CLineModelElementList::GetModelVelfitGroups( Int32 lineType )
{
    bool verbose = false;
    if(verbose)
    {
        Log.LogDebug("    model: group tags for lineType=%d", lineType);
    }
    std::vector<std::string> tags;
    std::vector<UInt32> nonGroupedLines;

    std::vector<UInt32> nonZeroIndexes = GetModelValidElementsIndexes();
    for(Int32 i=0; i<nonZeroIndexes.size(); i++)
    {
        Int32 iElts = nonZeroIndexes[i];
        Int32 nRays = m_Elements[iElts]->GetSize();
        for(Int32 iSubElts=0; iSubElts<nRays; iSubElts++)
        {
            if(verbose)
            {
                Log.LogDebug("    model: group tags - lineType=%d", lineType);
                Log.LogDebug("    model: group tags - m_Elements[iElts]->m_Rays[iSubElts].GetType()=%d", m_Elements[iElts]->m_Rays[iSubElts].GetType());
            }
            if(lineType == m_Elements[iElts]->m_Rays[iSubElts].GetType())
            {
                std::string _tag = m_Elements[iElts]->m_Rays[iSubElts].GetVelGroupName();
                if(_tag != "-1"){
                    tags.push_back(_tag);
                }else{
                    nonGroupedLines.push_back(iElts);
                }
            }
        }
    }

    if(verbose)
    {
        Log.LogDebug("    model: group tags non unique found=%d", tags.size());
        for( Int32 itag = 0; itag<tags.size(); itag++){
            Log.LogDebug("    model: non unique tag %d/%d = %s", itag+1, tags.size(), tags[itag].c_str());
        }
    }

    // create the group tag set by removing duplicates
    std::sort( tags.begin(), tags.end() );
    tags.erase( std::unique( tags.begin(), tags.end() ), tags.end() );
    //*
    if(verbose)
    {
        //print the tags
        for( Int32 itag = 0; itag<tags.size(); itag++){
            Log.LogDebug("    model: velfit group Tag %d/%d = %s", itag+1, tags.size(), tags[itag].c_str());
        }
    }
    //*/

    //add the grouped lines
    std::vector<std::vector<Int32>> groups;
    std::vector<std::string> groupsTags;
    for( Int32 itag = 0; itag<tags.size(); itag++)
    {
        std::vector<Int32> _group;
        for(Int32 i=0; i<nonZeroIndexes.size(); i++)
        {
            Int32 iElts = nonZeroIndexes[i];
            Int32 nRays = m_Elements[iElts]->GetSize();
            for(Int32 iSubElts=0; iSubElts<nRays; iSubElts++)
            {
                if(lineType == m_Elements[iElts]->m_Rays[iSubElts].GetType())
                {
                    std::string _tag = m_Elements[iElts]->m_Rays[iSubElts].GetVelGroupName();
                    if(_tag == tags[itag]){
                        _group.push_back(iElts);
                    }
                }
            }
        }
        //add the grouped lines, no duplicates
        std::sort( _group.begin(), _group.end() );
        _group.erase( std::unique( _group.begin(), _group.end() ), _group.end() );

        groups.push_back(_group);
        groupsTags.push_back(tags[itag]);
    }
    //add the non grouped lines, no duplicates
    std::sort( nonGroupedLines.begin(), nonGroupedLines.end() );
    nonGroupedLines.erase( std::unique( nonGroupedLines.begin(), nonGroupedLines.end() ), nonGroupedLines.end() );
    for( Int32 i = 0; i<nonGroupedLines.size(); i++)
    {
        std::vector<Int32> _group;
        _group.push_back(nonGroupedLines[i]);
        groups.push_back(_group);
        groupsTags.push_back("-1");
    }

    if(true)
    {
        //print the groups
        for( Int32 igr = 0; igr<groups.size(); igr++){
            Log.LogDebug("    model: Group %d/%d: nlines=%d, tag=%s", igr+1, groups.size(), groups[igr].size(), groupsTags[igr].c_str());
            for(Int32 i=0; i<groups[igr].size(); i++)
            {
                Log.LogDebug("    model: \t%d: iElt=%d", i+1, groups[igr][i]);
            }
        }
    }

    //Override velGroups from Catalog: => Individual lines as groups
    /*
    std::vector<std::vector<Int32>> groups;
    std::vector<Int32> nonZeroIndexes = GetModelValidElementsIndexes();
    for(Int32 i=0; i<nonZeroIndexes.size(); i++)
    {
        if(lineType == m_Elements[nonZeroIndexes[i]]->m_Rays[0].GetType())
        {
            std::vector<Int32> gr;
            gr.push_back(nonZeroIndexes[i]);
            groups.push_back(gr);
            //Log.LogInfo("Group %d, idx=%d", groups.size(), groups[groups.size()-1][0]);
        }
    }
    //*/

    return groups;
}

/**
 * \brief Returns a sorted, de-duplicated list of indices of lines whose support overlap ind's support and are not listed in the argument excludedInd.
 **/

std::vector<UInt32> CLineModelElementList::getOverlappingElements(UInt32 ind, const std::vector<UInt32> & excludedInd,Float64 redshift,Float64 overlapThres)
{
    std::vector<UInt32> indexes;

    if(m_Elements[ind]->IsOutsideLambdaRange()){
        indexes.push_back(ind);
        return indexes;
    }

    std::vector<CRay> raysRef = m_Elements[ind]->GetRays();
    Int32 linetypeRef = raysRef[0].GetType();

    Int32 xinf=0;
    Int32 yinf=0;
    Int32 xsup=0;
    Int32 ysup=0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        //check linetype
        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetypeRef){
            continue;
        }

        //check if outside lambdarange
        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            continue;
        }

        //check if in exclusion list
        bool excluded=false;
        for( UInt32 iexcl=0; iexcl<excludedInd.size(); iexcl++ )
        {
            if(iElts == excludedInd[iexcl]){
               excluded = true;
               break;
            }
        }
        if(excluded){
            continue;
        }

        std::vector<CRay> raysElt = m_Elements[iElts]->GetRays();

        for( UInt32 iRayElt=0; iRayElt<raysElt.size(); iRayElt++ )
        {
            for( UInt32 iRayRef=0; iRayRef<raysRef.size(); iRayRef++ )
            {
                Float64 muRef = raysRef[iRayRef].GetPosition()*(1+redshift);
                std::shared_ptr<const CLineProfile> profileRef = raysRef[iRayRef].GetProfile();
                Float64 cRef = m_Elements[ind]->GetLineWidth(muRef, redshift, raysRef[iRayRef].GetIsEmission());
                Float64 winsizeRef = profileRef->GetNSigmaSupport()*cRef;
                Float64 overlapSizeMin = winsizeRef*overlapThres;
                xinf = muRef-winsizeRef/2.0;
                xsup = muRef+winsizeRef/2.0;

                Float64 muElt = raysElt[iRayElt].GetPosition()*(1+redshift);
                std::shared_ptr<const CLineProfile> profileElt = raysElt[iRayElt].GetProfile();
                Float64 cElt = m_Elements[iElts]->GetLineWidth(muElt, redshift, raysElt[iRayElt].GetIsEmission());
                Float64 winsizeElt = profileElt->GetNSigmaSupport()*cElt;
                yinf = muElt-winsizeElt/2.0;
                ysup = muElt+winsizeElt/2.0;

                Float64 max = std::max(xinf,yinf);
                Float64 min = std::min(xsup,ysup);
                if( max-min < -overlapSizeMin ){
                    indexes.push_back(iElts);
                    break;
                }
            }
        }
    }

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    return indexes;
}


/**
 * \brief If argument j is a valid index of m_Elements, updates the element in that index calling its SetFittedAmplitude with arguments a and snr.
 **/
void CLineModelElementList::SetElementAmplitude(Int32 j, Float64 a, Float64 snr)
{
    if(j>=0 && j<m_Elements.size())
    {
        m_Elements[j]->SetFittedAmplitude(a, snr);
    }
    return;
}

/**
 * \brief If j is a valid index of m_Elements, returns a call to that element's GetElementAmplitude. If not, returns -1.
 **/
Float64 CLineModelElementList::GetElementAmplitude(Int32 j)
{
    Float64 a=-1.0;
    if(j>=0 && j<m_Elements.size())
    {
        a = m_Elements[j]->GetElementAmplitude();
    }
    return a;
}


void CLineModelElementList::SetSourcesizeDispersion(Float64 sizeArcsec)
{
    for(Int32 j=0; j<m_Elements.size(); j++)
    {
        m_Elements[j]->SetSourcesizeDispersion(sizeArcsec);
    }
}

/**
 * \brief Returns the first index of m_Elements where calling the element's FindElementIndex method with LineCatalogIndex argument does not return -1.
 **/
Int32 CLineModelElementList::FindElementIndex(Int32 LineCatalogIndex)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->FindElementIndex(LineCatalogIndex) !=-1){
            idx = iElts;
            break;
        }
    }
    return idx;
}

/**
 * \brief Returns the first index of m_Elements where calling the element's FindElementIndex method with LineTagStr argument does not return -1.
 **/
Int32 CLineModelElementList::FindElementIndex(std::string LineTagStr, Int32 linetype, Int32& lineIdx )
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        lineIdx = m_Elements[iElts]->FindElementIndex(LineTagStr) ;
        if( lineIdx!=-1 ){
            if( linetype!=-1 )
            {
                if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetType() != linetype){
                    continue;
                }
            }
            idx = iElts;
            break;
        }
    }
    return idx;
}

/**
 * \brief Returns the error of the support for subelements under the element with the argument eltId as index.
 * Accumulate "fit", the squared difference between model and spectrum, divided by the square of the m_ErrorNoContinuum value.
 * Accumulate "sumErr" 1 / square of the m_ErrorNoContinuum value.
 * return the square root of fit / sumErr.
 **/
Float64 CLineModelElementList::getModelErrorUnderElement( UInt32 eltId ,
					      const CSpectrumFluxAxis& spcFluxAxis,
					      const CSpectrumFluxAxis& modelFluxAxis
					      )
{
  //before elementlistcutting this variable was CElementList::m_ErrorNoContinuum, a reference initialized twice in CElementList constructor, first init to m_spcFluxAxisNoContinuum.GetError() and after to spectrumFluxAxis.GetError
  const CSpectrumNoiseAxis& errorNoContinuum=spcFluxAxis.GetError();

  
    if(eltId<0)
    {
        return -1.0;
    }

    Int32 numDevs = 0;
    Float64 fit = 0.0;
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 sumErr=0.0;

    TInt32RangeList support;
    UInt32 iElts=eltId;
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            return 0.0;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            support.push_back(s[iS]);
        }
    }


    Float64 w=0.0;
    for( UInt32 iS=0; iS<support.size(); iS++ )
    {
        for( UInt32 j=support[iS].GetBegin(); j<support[iS].GetEnd(); j++ )
        {
            numDevs++;
            diff = (Yspc[j] - Ymodel[j]);
            w = 1.0 / (errorNoContinuum[j]*errorNoContinuum[j]);
            fit += (diff*diff) * w;
            sumErr += w;
        }
    }
    return sqrt(fit/sumErr);
}

/**
 * \brief Returns a sorted set of samples indices present in the supports of the argument.
 * For each EltsIdx entry, if the entry is not outside lambda range, get the support of each subelement.
 * For each selected support, get the sample index. Sort this list and remove multiple entries. Return this clean list.
 **/
std::vector<UInt32> CLineModelElementList::getSupportIndexes( const std::vector<UInt32> & EltsIdx )
{
    std::vector<UInt32> indexes;

    TInt32RangeList support;
    for( UInt32 i=0; i<EltsIdx.size(); i++ )
    {
        Int32 iElts = EltsIdx[i];

        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            continue;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            support.push_back(s[iS]);
        }
    }

    for( UInt32 iS=0; iS<support.size(); iS++ )
    {
        for( UInt32 j=support[iS].GetBegin(); j<support[iS].GetEnd(); j++ )
        {
            indexes.push_back(j);
        }
    }

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    return indexes;
}


Int32 CLineModelElementList::getIndexAmpOffset(UInt32 xIndex)
{
  Int32 idxAmpOffset = -1;

  for( Int32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
      if(xIndex >= m_ampOffsetsIdxStart[i] && xIndex<=m_ampOffsetsIdxStop[i])
	{
	  idxAmpOffset=i;
	  break;
	}
    }
  Log.LogDetail("AmplitudeOffset enabled: idx offset = %d", idxAmpOffset);
  return idxAmpOffset;
}
    
void CLineModelElementList::setAmplitudeOffsetsCoeffsAt(UInt32 index,Float64 x0,Float64 x1,Float64 x2)
{
  m_ampOffsetsX0[index] = x0;
  m_ampOffsetsX1[index] = x1;
  m_ampOffsetsX2[index] = x2;
}

bool CLineModelElementList::addToSpectrumAmplitudeOffset(const CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& modelfluxAxis )
{


    Log.LogDetail( "Elementlist: Adding n=%d ampOffsets", m_ampOffsetsIdxStart.size());
    for( UInt32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
        for( UInt32 k=m_ampOffsetsIdxStart[i]; k<=m_ampOffsetsIdxStop[i]; k++ )
        {
            modelfluxAxis[k] += m_ampOffsetsX0[i] + m_ampOffsetsX1[i]*spectralAxis[k]+ m_ampOffsetsX2[i]*spectralAxis[k]*spectralAxis[k];
        }

    }
    return true;
}

Int32 CLineModelElementList::prepareAmplitudeOffset(const CSpectrumFluxAxis& spcFlux)
{
    m_ampOffsetsX0.clear();
    m_ampOffsetsX1.clear();
    m_ampOffsetsX2.clear();
    m_ampOffsetsIdxStart.clear();
    m_ampOffsetsIdxStop.clear();

    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    if(validEltsIdx.size()<1)
    {
        return -1;
    }
    std::vector<UInt32> supportIdxes = getSupportIndexes( validEltsIdx );
    if(supportIdxes.size()<1)
    {
        return -1;
    }

    Int32 idxPrevious = supportIdxes[0];
    m_ampOffsetsIdxStart.push_back(supportIdxes[0]);
    for( UInt32 i=1; i<supportIdxes.size(); i++ )
    {
        UInt32 idxCurrent = supportIdxes[i];
        if(idxCurrent>idxPrevious+1)
        {
            m_ampOffsetsIdxStop.push_back(idxPrevious);
            m_ampOffsetsIdxStart.push_back(idxCurrent);
        }
        idxPrevious = idxCurrent;
    }
    m_ampOffsetsIdxStop.push_back(supportIdxes[supportIdxes.size()-1]);

    /*
    //estimate mean fluxNoCOnt in each support
    for( Int32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
        Float64 sum = 0.0;
        Int32 count = 0;

        for( Int32 k=m_ampOffsetsIdxStart[i]; k<=m_ampOffsetsIdxStop[i]; k++ )
        {
            sum +=spcFlux[k];
            count +=1;
        }
        Float64 mean = 0.0;
        if(count>0)
        {
            mean = sum/(Float64)count;
        }

        m_ampOffsetsA.push_back(mean);
    }
    //*/
    for( UInt32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
        m_ampOffsetsX0.push_back(0.0);
        m_ampOffsetsX1.push_back(0.0);
        m_ampOffsetsX2.push_back(0.0);
    }

    if(1)
    {
        for( UInt32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
        {
            Log.LogDebug( "Elementlist: i=%d, m_ampOffsetsIdxStart: %d, m_ampOffsetsIdxStop: %d", i, m_ampOffsetsIdxStart[i], m_ampOffsetsIdxStop[i] );
        }
    }
    return 0;
}


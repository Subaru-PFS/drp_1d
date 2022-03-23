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
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/line/catalogsTplShape.h"
#include "RedshiftLibrary/line/linetags.h"

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <string>
#include <fstream>
#include <iostream>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


bool CLineCatalogsTplShape::Init(Int32 enableISMCalzetti, 
                                std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti,
                                Float64 nsigmasupport)
{
    m_nsigmasupport = nsigmasupport;
    m_opt_dust_calzetti = enableISMCalzetti;
    m_ismCorrectionCalzetti = ismCorrectionCalzetti;

    return true;
}



/**
 * @brief CLineCatalogsTplShape::GetRestLinesList
 * @param index
 * WARNING: ismCoeff not applied on the restlines provided by that function.
 */
CLineCatalog::TLineVector CLineCatalogsTplShape::GetRestLinesList( Int32 index ) const
{
    Int32 typeFilter=-1;
    Int32 forceFilter=-1;

    CLineCatalog::TLineVector restLineList = m_lineRatioCatalogs[index].GetFilteredList( typeFilter, forceFilter);
    return restLineList;
}

Int32 CLineCatalogsTplShape::GetCatalogsCount() const
{
    return m_lineRatioCatalogs.size();
}

const TFloat64List&  CLineCatalogsTplShape::getCatalogsPriors()
{
  if(m_catalogsPriors.empty())
    {
      for(CLineRatioCatalog cat : m_lineRatioCatalogs) m_catalogsPriors.push_back(cat.getPrior());
    }
  return m_catalogsPriors;
}

std::string CLineCatalogsTplShape::GetCatalogName(Int32 idx) const
{
  return m_lineRatioCatalogs[idx].getName();
}

Float64 CLineCatalogsTplShape::GetIsmCoeff(Int32 idx) const
{
  return m_ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(idx));
} 

Int32 CLineCatalogsTplShape::GetIsmIndex(Int32 idx) const
{
  return m_lineRatioCatalogs[idx].getIsmIndex();
}

bool CLineCatalogsTplShape::GetCatalogVelocities(Int32 idx, Float64& elv, Float64& alv ) const
{
  //TODO generic velocity groups : there should not be hardcoded values, this should return a map
  elv = m_lineRatioCatalogs[idx].getVelocity("em_vel");
    alv = m_lineRatioCatalogs[idx].getVelocity("abs_vel");
    return true;
}


bool CLineCatalogsTplShape::InitLineCorrespondingAmplitudes(const CLineModelElementList &LineModelElementList)
{
    //first set all corresponding amplitudes to 0.0;
    for( Int32 iElts=0; iElts<LineModelElementList.size(); iElts++ )
    {
        Int32 nLines = LineModelElementList[iElts]->GetSize();
        TFloat64List thisCatLinesCorresp(nLines, 0.0); //is nLines cte among Elts?
        std::vector<TFloat64List> thisElementLinesCorresp(GetCatalogsCount(), thisCatLinesCorresp);

        m_LineCatalogLinesCorrespondingNominalAmp.push_back(thisElementLinesCorresp);

        //now set the non-zero amp correspondences
        for(Int32 iCatalog=0; iCatalog<GetCatalogsCount(); iCatalog++)
        {
            CLineCatalog::TLineVector currentCatalogLineList = m_lineRatioCatalogs[iCatalog].GetList();
            for(Int32 kL=0; kL<currentCatalogLineList.size(); kL++)
            {
                Float64 nominalAmp = currentCatalogLineList[kL].GetNominalAmplitude();
                Float64 restLambda = currentCatalogLineList[kL].GetPosition();
                Float64 dustCoeff = m_ismCorrectionCalzetti->GetDustCoeff( GetIsmIndex(iCatalog), restLambda);
                nominalAmp*=dustCoeff;
                //find line in the elementList
                Int32 nLines = LineModelElementList[iElts]->GetSize();
                for(Int32 j=0; j<nLines; j++){

                    if(LineModelElementList[iElts]->m_Lines[j].GetName() == currentCatalogLineList[kL].GetName())
                    {
                        m_LineCatalogLinesCorrespondingNominalAmp[iElts][iCatalog][j]=nominalAmp;
                    }
                }
            }
        }
    }

    //Now log the linesCorrespondingNominalAmp
    for( Int32 iElts=0; iElts<LineModelElementList.size(); iElts++ )
    {
        for(Int32 k=0; k<GetCatalogsCount(); k++)
        {
	  Log.LogDebug(Formatter()<<"log linesCorrespondingNominalAmp for "<<m_lineRatioCatalogs[k].getName()); 
            Int32 nLines = LineModelElementList[iElts]->GetSize();
            for(Int32 j=0; j<nLines; j++){
	        Float64 ebv = m_ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(k));
                Float64 nomAmp = m_LineCatalogLinesCorrespondingNominalAmp[iElts][k][j];
                std::string lineName = LineModelElementList[iElts]->m_Lines[j].GetName();
                Log.LogDebug("    CatalogsTplShape - linesCorrespondingNominalAmp iElt=%d, iCatalog=%d, iLine=%d with name=%s, ebv=%f: NominalAmpFound = %e", iElts, k, j, lineName.c_str(), ebv, nomAmp);
            }
        }
    }

    return 0;
}

const CLineCatalog& CLineCatalogsTplShape::GetCatalog(Int32 iCatalog) const
{
    return m_lineRatioCatalogs[iCatalog];
}

/**
 * \brief Calculates the best fit between the linemodel fitted amplitudes and the tplShaped catalogs: (for lm-rigidity=tplcorr)
 *
 **/
Float64 CLineCatalogsTplShape::GetBestFit( const CLineCatalog::TLineVector& restLineList, const TFloat64List &fittedAmplitudes, const TFloat64List & fittedErrors, TFloat64List &amplitudesCorrected, std::string& bestTplName) const
{
    Float64 coeffMin = -1;
    TInt32List mask;
    TFloat64List bestFitAmplitudes;
    //bestFitAmplitudes.resize(restLineList.size());
    for(Int32 iCatalogs=0; iCatalogs<m_lineRatioCatalogs.size(); iCatalogs++)
    {
        CLineCatalog::TLineVector currentCatalogLineList = m_lineRatioCatalogs[iCatalogs].GetList();

        //create the amplitude float vectors
        TFloat64List tplshapeAmplitudes;
        TFloat64List linemodelAmplitudes;
        TFloat64List linemodelErrors;

        mask.resize(fittedAmplitudes.size());

        for( Int32 iRestLine=0; iRestLine<restLineList.size(); iRestLine++ )
        {
            if(fittedAmplitudes[iRestLine]>=0 && fittedErrors[iRestLine]>0)
            {
                mask[iRestLine]=1;
                Int32 iTplshapeLineFound = -1;
                for(Int32 itplshapeLine=0; itplshapeLine<currentCatalogLineList.size(); itplshapeLine++)
                {
                    std::string tplshapeLineName =  currentCatalogLineList[itplshapeLine].GetName();
                    std::string restLineName = restLineList[iRestLine].GetName();
                    if ( restLineName==tplshapeLineName ){
                        iTplshapeLineFound = itplshapeLine;
                        break;
                    }
                }

                if ( iTplshapeLineFound < 0 )
                {
                    tplshapeAmplitudes.push_back(0.0);
                }else
                {
                    Float64 amp = currentCatalogLineList[iTplshapeLineFound].GetNominalAmplitude();
                    tplshapeAmplitudes.push_back(amp);
                }
                linemodelAmplitudes.push_back(fittedAmplitudes[iRestLine]);
                linemodelErrors.push_back(fittedErrors[iRestLine]);
            }else
            {
                mask[iRestLine]=0;
            }
        }

        if(linemodelAmplitudes.size()>1 && linemodelAmplitudes.size()==tplshapeAmplitudes.size())
        {
            TFloat64List ampsCorrected;
            ampsCorrected.resize(linemodelAmplitudes.size());
            Float64 fit = GetFit(linemodelAmplitudes, linemodelErrors, tplshapeAmplitudes, ampsCorrected);
            if(fit>0.0 && !std::isnan(fit) && (fit<coeffMin || coeffMin==-1))
            {
                coeffMin = fit;
                bestFitAmplitudes = ampsCorrected;
                bestTplName = GetCatalogName(iCatalogs);
            }
        }

    }
    //coeff min normalization
    //coeffMin = sqrt(coeffMin);
    //coeffMin /= 1.0;

    //fill the corrected amplitudes vector
    Int32 iTplAmps = 0;
    for( Int32 iRestLine=0; iRestLine<restLineList.size(); iRestLine++ )
    {
        if(mask[iRestLine]>0 && coeffMin>=0)
        {
            amplitudesCorrected[iRestLine]=bestFitAmplitudes[iTplAmps];
            iTplAmps++;
        }else
        {
            amplitudesCorrected[iRestLine]=-1.0;
        }

    }

    return coeffMin;
}

Float64 CLineCatalogsTplShape::GetFit( const TFloat64List &ampsLM, const TFloat64List &errLM, const TFloat64List &ampsTPL, TFloat64List &ampsCorrected )  const
{

    Float64 N = ampsLM.size();
//    // Normalize AmpsLM first
//    Float64 normalizeCoeff = 0.0;
//    for(Int32 k=0; k<N; k++)
//    {
//        normalizeCoeff+= ampsLM[k];
//    }
//    for(Int32 k=0; k<N; k++)
//    {
//        ampsLM[k] = ampsLM[k]/normalizeCoeff;
//    }
    Float64 normalizeCoeff = 1.0;

    //estimate fitting amplitude
    Float64 sumLM = 0.0;
    Float64 sumTPL = 0.0;
    Float64 sumCross= 0.0;
    Float64 sumTPL2 = 0.0;

    for(Int32 k=0; k<N; k++)
    {
        Float64 err2 = 1.0 / (errLM[k] * errLM[k]);

        sumCross+=ampsLM[k]*ampsTPL[k]*err2;
        sumTPL2+=ampsTPL[k]*ampsTPL[k]*err2;

        sumLM += ampsLM[k]*err2;
        sumTPL += ampsTPL[k]*err2;
    }
//    if ( sumLM==0 || sumTPL==0 )
//    {
//        return -1.0;
//    }
//    Float64 ampl = sumLM / sumTPL;
    if ( sumCross==0 || sumTPL2==0 )
    {
        return -1.0;
    }
    Float64 ampl = sumCross / sumTPL2;

    Float64 fit=0.0;
    Float64 diff;
    for(Int32 k=0; k<N; k++)
    {
        Float64 err2 = 1.0 / (errLM[k] * errLM[k]);
        diff = ampsLM[k]-ampl*ampsTPL[k];
        fit += diff*diff*err2;
    }

    //fill the amps_corrected vector
    for(Int32 k=0; k<N; k++)
    {
        ampsCorrected[k] = ampsTPL[k]*ampl*normalizeCoeff;
    }

    return fit;
}



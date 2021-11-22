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
#ifndef _ELEMENTS_H
#define _ELEMENTS_H

#include "RedshiftLibrary/linemodel/element.h"

namespace NSEpic
{
  static Int32 defaultIdx = -1;
  class CElements 
  {
  private:
    std::vector<std::shared_ptr<CLineModelElement> > m_Elements;
    
  public:
    std::vector<Float64> m_ampOffsetsX0;
    std::vector<Float64> m_ampOffsetsX1;
    std::vector<Float64> m_ampOffsetsX2;
    std::vector<Int32> m_ampOffsetsIdxStart;
    std::vector<Int32> m_ampOffsetsIdxStop;

    std::vector<Int32> m_elementsDisabledIndexes;

    
    CElements(const CRayCatalog::TRayVector &restRayList);

    const CRayCatalog::TRayVector &m_RestRayList;
    
    std::vector<UInt32> GetModelValidElementsIndexes();

    void SetElementAmplitude(Int32 j, Float64 a, Float64 snr);
    Float64 GetElementAmplitude(Int32 j);

    std::vector<UInt32> getOverlappingElements(UInt32 ind , const std::vector<UInt32> & excludedInd,Float64 redshift, Float64 overlapThres=0.1);

    Int32 GetModelValidElementsNDdl();
    Int32 GetModelNonZeroElementsNDdl();

     void SetSourcesizeDispersion(Float64 sizeArcsec);
    
    std::vector<std::vector<Int32>> GetModelVelfitGroups(Int32 lineType );

    Int32 FindElementIndex(Int32 LineCatalogIndex);
    Int32 FindElementIndex(std::string LineTagStr, Int32 linetype=-1, Int32& lineIdx=defaultIdx);

    Float64 getModelErrorUnderElement( UInt32 eltId,const CSpectrumFluxAxis& spcFluxAxis,const CSpectrumFluxAxis& modelFluxAxis);

    std::vector<UInt32> getSupportIndexes(const std::vector<UInt32> & EltsIdx);

    Int32 getIndexAmpOffset(UInt32 index);
    void setAmplitudeOffsetsCoeffsAt(UInt32 index,Float64 x0,Float64 x1,Float64 x2);
    Int32 prepareAmplitudeOffset(const CSpectrumFluxAxis &spcFlux);
    Bool addToSpectrumAmplitudeOffset(const CSpectrumSpectralAxis& spectralAxis,CSpectrumFluxAxis &modelfluxAxis);

    bool IsElementIndexInDisabledList(Int32 index);
    void SetElementIndexesDisabledAuto();
    void ResetElementIndexesDisabled();
    
    const std::shared_ptr<CLineModelElement> &operator[](UInt32 i) const {return m_Elements[i];}
    UInt32 size() const {return m_Elements.size();}
    void push_back(std::shared_ptr<CLineModelElement> elt){ m_Elements.push_back(elt);}
  };
}
#endif

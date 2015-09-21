#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/linemodel/singleline.h>

#include <epic/core/debug/assert.h>
#include <epic/core/log/log.h>

#include <algorithm>

using namespace NSEpic;

CLineModelElementList::CLineModelElementList( const CSpectrum& spectrum, const CRayCatalog::TRayVector& restRayList)
{
    LoadCatalog(restRayList);
    m_RestRayList = restRayList;
    m_SpectrumModel = new CSpectrum(spectrum);
    m_SpcFluxAxis = spectrum.GetFluxAxis();
}

CLineModelElementList::~CLineModelElementList()
{
}

const CSpectrum& CLineModelElementList::GetModelSpectrum() const
{
    return *m_SpectrumModel;
}

void CLineModelElementList::LoadCatalog(const CRayCatalog::TRayVector& restRayList)
{
    Float64 nominalWidthDefault = 3.18; //suited to PFS RJLcont simulations

    //Load OIII lines
    std::vector<Int32> OIIIaIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1)");
    if(OIIIaIdx.size()==1){
        addSingleLine(restRayList[OIIIaIdx[0]], OIIIaIdx[0], nominalWidthDefault);
    }
    std::vector<Int32> OIIIbIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1/3)");
    if(OIIIbIdx.size()==1){
        addSingleLine(restRayList[OIIIbIdx[0]], OIIIbIdx[0], nominalWidthDefault);
    }

    //Load NII lines
    std::vector<Int32> NIIdx = findLineIdxInCatalog( restRayList, "[NII]");
    if(NIIdx.size()==2){
        addSingleLine(restRayList[NIIdx[0]], NIIdx[0], nominalWidthDefault);
        addSingleLine(restRayList[NIIdx[1]], NIIdx[1], nominalWidthDefault);
    }

    //Load OII line
    std::vector<Int32> OIIdx = findLineIdxInCatalog( restRayList, "[OII]");
    if(OIIdx.size()==1){
        addSingleLine(restRayList[OIIdx[0]], OIIdx[0], 3.55);
    }

    //Load the rest of the single lines
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if ( FindElementIndex(iRestRay)==-1 )
        {
            addSingleLine(restRayList[iRestRay], iRestRay, nominalWidthDefault);
        }
    }
}

void CLineModelElementList::fit(Float64 redshift, CLineModelResult::SLineModelSolution& modelSolution)
{
    //initialize the model
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
        modelFluxAxis[i] = 0.0;
    }

    //fit the model amplitudes
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->fitAmplitude(spectralAxis, m_SpcFluxAxis, redshift);
    }

    //eventually apply rules
    //...

    //create spectrum model
    modelSolution = GetModelSolution();
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, redshift);
    }
}

std::vector<int> CLineModelElementList::findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag)
{
    std::vector<Int32> indexes;
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        std::string name = restRayList[iRestRay].GetName();
        std::size_t foundstra = name.find(strTag.c_str());
        if (foundstra!=std::string::npos){
            indexes.push_back(iRestRay);
        }
    }
    return indexes;
}

void CLineModelElementList::addSingleLine(const CRay &r, Int32 index, Float64 nominalWidth)
{
    //CSingleLine line = CSingleLine(r, nominalWidth);
    std::vector<Int32> a;
    a.push_back(index);
    //CSingleLine c(r, nominalWidth, a);
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CSingleLine(r, nominalWidth, a)));
    //m_Elements.push_back(new CSingleLine(r, nominalWidth, a));
    m_LineCatalogIndexes.push_back(index);
}

void CLineModelElementList::applyRules()
{
    //    //*
    //    Float64 doublettol=0.2;
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "[OIII](doublet-1)", "[OIII](doublet-1/3)", 0.334*(1.0+doublettol));
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "[OIII](doublet-1/3)", "[OIII](doublet-1)", (1.0+doublettol)/0.334);
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "Halpha", "Hbeta", 1.0/2.86);
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "Hbeta", "Hgamma", 0.47);
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "Hgamma", "Hdelta", 1.0);
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "[NII](doublet-1)", "[NII](doublet-1/2.95)", (1.0+doublettol)/2.95);
    //    Apply2LinesAmplitudeRule(restRayList, modelSolution.Amplitudes, modelSolution.OutsideLambdaRange, "[NII](doublet-1/2.95)", "[NII](doublet-1)", 2.95*(1.0+doublettol));
    //    //*/
}

CLineModelResult::SLineModelSolution CLineModelElementList::GetModelSolution()
{
    CLineModelResult::SLineModelSolution modelSolution;
    for( UInt32 iRestRay=0; iRestRay<m_RestRayList.size(); iRestRay++ )
    {
        Int32 eIdx = FindElementIndex(iRestRay);
        Int32 subeIdx = m_Elements[eIdx]->FindElementIndex(iRestRay);

        //modelSolution.fittingIndexRange.push_back( m_Elements[eIdx].);
        modelSolution.Amplitudes.push_back(m_Elements[eIdx]->GetFittedAmplitude(subeIdx));
        //modelSolution.Widths.push_back(-1.0);
        //modelSolution.OutsideLambdaRange.push_back(true);
    }

    return modelSolution;
}

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


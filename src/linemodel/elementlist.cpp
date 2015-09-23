#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/linemodel/singleline.h>
#include <epic/redshift/linemodel/multiline.h>

#include <epic/core/debug/assert.h>
#include <epic/core/log/log.h>

#include <algorithm>

using namespace NSEpic;

CLineModelElementList::CLineModelElementList( const CSpectrum& spectrum, const CRayCatalog::TRayVector& restRayList)
{
    //LoadCatalog(restRayList);
    LoadCatalogMultilineBalmer(restRayList);
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
    Float64 nominalWidthDefault = 3.3; //suited to PFS RJLcont simulations

    //Load OIII lines
    std::vector<Int32> OIIIaIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1)");
    std::vector<Int32> OIIIbIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1/3)");
    if(OIIIaIdx.size()==1 && OIIIbIdx.size()==1){
        addDoubleLine(restRayList[OIIIaIdx[0]], restRayList[OIIIbIdx[0]], OIIIaIdx[0], OIIIbIdx[0], nominalWidthDefault, 1.0, 1.0/3.0);
    }else{
        if(OIIIaIdx.size()==1){
            addSingleLine(restRayList[OIIIaIdx[0]], OIIIaIdx[0], nominalWidthDefault);
        }
        if(OIIIbIdx.size()==1){
            addSingleLine(restRayList[OIIIbIdx[0]], OIIIbIdx[0], nominalWidthDefault);
        }
    }

    //Load NII lines
    std::vector<Int32> NII1dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1)");
    std::vector<Int32> NII2dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1/2.95)");
    if(NII1dx.size()==1 && NII2dx.size()==1){
        addDoubleLine(restRayList[NII1dx[0]], restRayList[NII2dx[0]], NII1dx[0], NII2dx[0], nominalWidthDefault, 1.0, 1.0/2.95);
    }else{
        if(NII1dx.size()==1){
            addSingleLine(restRayList[NII1dx[0]], NII1dx[0], nominalWidthDefault);
        }
        if(NII2dx.size()==1){
            addSingleLine(restRayList[NII2dx[0]], NII2dx[0], nominalWidthDefault);
        }
    }

    //Load OII line doublet
    std::vector<Int32> OII1dx = findLineIdxInCatalog( restRayList, "[OII](doublet-1)");
    std::vector<Int32> OII2dx = findLineIdxInCatalog( restRayList, "[OII](doublet-1/3)");
    if(OII1dx.size()==1 && OII2dx.size()==1){
        addDoubleLine(restRayList[OII1dx[0]], restRayList[OII2dx[0]], OII1dx[0], OII2dx[0], 3.2, 1.0, 1.0/2.9);
    }else{
        if(OII1dx.size()==1){
            addSingleLine(restRayList[OII1dx[0]], OII1dx[0], nominalWidthDefault);
        }
        if(OII2dx.size()==1){
            addSingleLine(restRayList[OII2dx[0]], OII2dx[0], nominalWidthDefault);
        }
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

void CLineModelElementList::LoadCatalogMultilineBalmer(const CRayCatalog::TRayVector& restRayList)
{
    Float64 nominalWidthDefault = 3.3; //suited to PFS RJLcont simulations

    //Load OIII lines
    std::vector<Int32> OIIIaIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1)");
    std::vector<Int32> OIIIbIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1/3)");
    if(OIIIaIdx.size()==1 && OIIIbIdx.size()==1){
        addDoubleLine(restRayList[OIIIaIdx[0]], restRayList[OIIIbIdx[0]], OIIIaIdx[0], OIIIbIdx[0], nominalWidthDefault, 1.0, 1.0/3.0);
    }else{
        if(OIIIaIdx.size()==1){
            addSingleLine(restRayList[OIIIaIdx[0]], OIIIaIdx[0], nominalWidthDefault);
        }
        if(OIIIbIdx.size()==1){
            addSingleLine(restRayList[OIIIbIdx[0]], OIIIbIdx[0], nominalWidthDefault);
        }
    }

    //Load NII lines
    std::vector<Int32> NII1dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1)");
    std::vector<Int32> NII2dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1/2.95)");
    if(NII1dx.size()==1 && NII2dx.size()==1){
        addDoubleLine(restRayList[NII1dx[0]], restRayList[NII2dx[0]], NII1dx[0], NII2dx[0], nominalWidthDefault, 1.0, 1.0/2.95);
    }else{
        if(NII1dx.size()==1){
            addSingleLine(restRayList[NII1dx[0]], NII1dx[0], nominalWidthDefault);
        }
        if(NII2dx.size()==1){
            addSingleLine(restRayList[NII2dx[0]], NII2dx[0], nominalWidthDefault);
        }
    }

    //Load OII line doublet
    std::vector<Int32> OII1dx = findLineIdxInCatalog( restRayList, "[OII](doublet-1)");
    std::vector<Int32> OII2dx = findLineIdxInCatalog( restRayList, "[OII](doublet-1/3)");
    if(OII1dx.size()==1 && OII2dx.size()==1){
        addDoubleLine(restRayList[OII1dx[0]], restRayList[OII2dx[0]], OII1dx[0], OII2dx[0], 3.2, 1.0, 1.0/2.9);
    }else{
        if(OII1dx.size()==1){
            addSingleLine(restRayList[OII1dx[0]], OII1dx[0], nominalWidthDefault);
        }
        if(OII2dx.size()==1){
            addSingleLine(restRayList[OII2dx[0]], OII2dx[0], nominalWidthDefault);
        }
    }

    //Load Balmer multilines
    std::vector<Int32> Halphaidx = findLineIdxInCatalog( restRayList, "Halpha");
    std::vector<Int32> Hbetaidx = findLineIdxInCatalog( restRayList, "Hbeta");
    std::vector<Int32> Hgammaidx = findLineIdxInCatalog( restRayList, "Hgamma");
    std::vector<Int32> Hdeltaidx = findLineIdxInCatalog( restRayList, "Hdelta");
    if(Halphaidx.size()==1 && Hbetaidx.size()==1 && Hgammaidx.size()==1 && Hdeltaidx.size()==1){
        std::vector<CRay> lines;
        lines.push_back(restRayList[Halphaidx[0]]);
        lines.push_back(restRayList[Hbetaidx[0]]);
        lines.push_back(restRayList[Hgammaidx[0]]);
        lines.push_back(restRayList[Hdeltaidx[0]]);

        std::vector<Float64> amps;
        amps.push_back(1.0);
        amps.push_back(0.190);
        amps.push_back(0.071);
        amps.push_back(0.035);

        std::vector<Int32> inds;
        inds.push_back(Halphaidx[0]);
        inds.push_back(Hbetaidx[0]);
        inds.push_back(Hgammaidx[0]);
        inds.push_back(Hdeltaidx[0]);

        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, amps, nominalWidthDefault, inds)));
    }else{
        if(Halphaidx.size()==1){
            addSingleLine(restRayList[Halphaidx[0]], Halphaidx[0], nominalWidthDefault);
        }
        if(Hbetaidx.size()==1){
            addSingleLine(restRayList[Hbetaidx[0]], Hbetaidx[0], nominalWidthDefault);
        }
        if(Hgammaidx.size()==1){
            addSingleLine(restRayList[Hgammaidx[0]], Hgammaidx[0], nominalWidthDefault);
        }
        if(Hdeltaidx.size()==1){
            addSingleLine(restRayList[Hdeltaidx[0]], Hdeltaidx[0], nominalWidthDefault);
        }
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

void CLineModelElementList::LoadCatalogSingleLines(const CRayCatalog::TRayVector& restRayList)
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

    //eventually apply rules,
    // WARNING: no noise taken into account for now...
    //applyRules();

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
}

void CLineModelElementList::addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2)
{
    std::vector<CRay> lines;
    lines.push_back(r1);
    lines.push_back(r2);

    std::vector<Float64> amps;
    amps.push_back(a1);
    amps.push_back(a2);

    std::vector<Int32> a;
    a.push_back(index1);
    a.push_back(index2);
    //CSingleLine c(r, nominalWidth, a);
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, amps, nominalWidth, a)));
    //m_Elements.push_back(new CSingleLine(r, nominalWidth, a));
}

void CLineModelElementList::applyRules()
{
    //todo: check for the noise when applying amplitudes rules...

    Apply2SingleLinesAmplitudeRule("Halpha", "Hbeta", 1.0/2.86);
    Apply2SingleLinesAmplitudeRule("Hbeta", "Hgamma", 0.47);
    Apply2SingleLinesAmplitudeRule("Hgamma", "Hdelta", 1.0);
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



Void CLineModelElementList::Apply2SingleLinesAmplitudeRule(std::string lineA, std::string lineB, Float64 coeff )
{
    Int32 iA = FindElementIndex(lineA);
    if(m_Elements[iA]->GetSize()>1){
        iA=-1;
    }
    Int32 iB = FindElementIndex(lineB);
    if(m_Elements[iB]->GetSize()>1){
        iB=-1;
    }
    if(iA==-1 || iB==-1 || iA==iB){
        return;
    }

    Float64 ampA = m_Elements[iA]->GetFittedAmplitude(0);
    m_Elements[iB]->LimitFittedAmplitude(0, coeff*ampA);
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

Int32 CLineModelElementList::FindElementIndex(std::string LineTagStr)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->FindElementIndex(LineTagStr) !=-1){
            idx = iElts;
            break;
        }
    }
    return idx;
}


#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/linemodel/singleline.h>
#include <epic/redshift/linemodel/multiline.h>

#include <epic/redshift/gaussianfit/multigaussianfit.h>
#include <gsl/gsl_multifit.h>

#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/spectrum/template/template.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


#include <epic/redshift/continuum/irregularsamplingmedian.h>
#include <epic/redshift/spectrum/io/fitswriter.h>

#include <epic/core/debug/assert.h>
#include <epic/core/log/log.h>


#include <boost/format.hpp>
#include <boost/chrono/thread_clock.hpp>

#include <algorithm>

using namespace NSEpic;

CLineModelElementList::CLineModelElementList(const CSpectrum& spectrum, const CSpectrum &spectrumNoContinuum, const CRayCatalog::TRayVector& restRayList, Int32 widthType)
{
    m_LineWidthType = widthType;

    //PFS
    m_nominalWidthDefaultEmission = 3.4;//3.4; //suited to PFS RJLcont simulations
    m_nominalWidthDefaultAbsorption = m_nominalWidthDefaultEmission;

    //LoadCatalog(restRayList);
    LoadCatalog_tplExtendedBlue(restRayList);
    //LoadCatalogMultilineBalmer(restRayList);
    //LoadCatalogSingleLines(restRayList);

    //LogCatalogInfos();
    m_RestRayList = restRayList;
    m_SpectrumModel = new CSpectrum(spectrum);

    m_SpcFluxAxis = spectrum.GetFluxAxis();
    m_SpcNoContinuumFluxAxis = spectrumNoContinuum.GetFluxAxis();
    m_ContinuumFluxAxis.SetSize( m_SpcFluxAxis.GetSamplesCount() );
    m_precomputedFineGridContinuumFlux = NULL;
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
    //Load OIII multilines
    std::vector<Int32> OIIIaIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> OIIIbIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1/3)", CRay::nType_Emission);
    if(OIIIaIdx.size()==1 && OIIIbIdx.size()==1){
        addDoubleLine(restRayList[OIIIaIdx[0]], restRayList[OIIIbIdx[0]], OIIIaIdx[0], OIIIbIdx[0], m_nominalWidthDefaultEmission, 1.0, 1.0/3.0);
    }else{
        if(OIIIaIdx.size()==1){
            addSingleLine(restRayList[OIIIaIdx[0]], OIIIaIdx[0], m_nominalWidthDefaultEmission);
        }
        if(OIIIbIdx.size()==1){
            addSingleLine(restRayList[OIIIbIdx[0]], OIIIbIdx[0], m_nominalWidthDefaultEmission);
        }
    }

    //Load NII multilines
    std::vector<Int32> NII1dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> NII2dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1/2.95)", CRay::nType_Emission);
    if(NII1dx.size()==1 && NII2dx.size()==1){
        addDoubleLine(restRayList[NII1dx[0]], restRayList[NII2dx[0]], NII1dx[0], NII2dx[0], m_nominalWidthDefaultEmission, 1.0, 1.0/2.95);
    }else{
        if(NII1dx.size()==1){
            addSingleLine(restRayList[NII1dx[0]], NII1dx[0], m_nominalWidthDefaultEmission);
        }
        if(NII2dx.size()==1){
            addSingleLine(restRayList[NII2dx[0]], NII2dx[0], m_nominalWidthDefaultEmission);
        }
    }

    //Load OII line doublet
    std::vector<Int32> OII1dx = findLineIdxInCatalog( restRayList, "[OII]3729", CRay::nType_Emission);
    std::vector<Int32> OII2dx = findLineIdxInCatalog( restRayList, "[OII]3726", CRay::nType_Emission);
    if(OII1dx.size()==1 && OII2dx.size()==1){
        //addDoubleLine(restRayList[OII1dx[0]], restRayList[OII2dx[0]], OII1dx[0], OII2dx[0], m_nominalWidthDefaultEmission, 1.0, 2.9);
        addDoubleLine(restRayList[OII1dx[0]], restRayList[OII2dx[0]], OII1dx[0], OII2dx[0], m_nominalWidthDefaultEmission, 1.0, 1.0);
    }else{
        if(OII1dx.size()==1){
            addSingleLine(restRayList[OII1dx[0]], OII1dx[0], m_nominalWidthDefaultEmission);
        }
        if(OII2dx.size()==1){
            addSingleLine(restRayList[OII2dx[0]], OII2dx[0], m_nominalWidthDefaultEmission);
        }
    }


    //Load FeII 2586/2600 multilines
    std::vector<Int32> aIdx = findLineIdxInCatalog( restRayList, "FeII2600", CRay::nType_Absorption);
    std::vector<Int32> bIdx = findLineIdxInCatalog( restRayList, "FeII2586", CRay::nType_Absorption);
    std::vector<Int32> cIdx = findLineIdxInCatalog( restRayList, "FeII2382", CRay::nType_Absorption);
    std::vector<Int32> dIdx = findLineIdxInCatalog( restRayList, "FeII2374", CRay::nType_Absorption);
    std::vector<Int32> eIdx = findLineIdxInCatalog( restRayList, "FeII2344", CRay::nType_Absorption);
    std::vector<Int32> fIdx = findLineIdxInCatalog( restRayList, "FeII2260", CRay::nType_Absorption);
    std::vector<Int32> gIdx = findLineIdxInCatalog( restRayList, "FeII2249", CRay::nType_Absorption);
    std::vector<Int32> hIdx = findLineIdxInCatalog( restRayList, "FeII1608", CRay::nType_Absorption);
    std::vector<CRay> lines;
    std::vector<Float64> amps;
    std::vector<Int32> inds;
//    if(aIdx.size()==1){
//        lines.push_back(restRayList[aIdx[0]]);
//        amps.push_back(0.239);//2600
//        inds.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        lines.push_back(restRayList[bIdx[0]]);
//        amps.push_back(0.0691); //2586
//        inds.push_back(bIdx[0]);
//    }
//    if(cIdx.size()==1){
//        lines.push_back(restRayList[cIdx[0]]);
//        amps.push_back(0.32); //2382
//        inds.push_back(cIdx[0]);
//    }
//    if(dIdx.size()==1){
//        lines.push_back(restRayList[dIdx[0]]);
//        amps.push_back(0.0313); //2374
//        inds.push_back(dIdx[0]);
//    }
//    if(eIdx.size()==1){
//        lines.push_back(restRayList[eIdx[0]]);
//        amps.push_back(0.114); //2344
//        inds.push_back(eIdx[0]);
//    }
//    if(fIdx.size()==1){
//        lines.push_back(restRayList[fIdx[0]]);
//        amps.push_back(0.00244); //2260
//        inds.push_back(fIdx[0]);
//    }
//    if(gIdx.size()==1){
//        lines.push_back(restRayList[gIdx[0]]);
//        amps.push_back(0.001821);//2249
//        inds.push_back(gIdx[0]);
//    }
//    if(hIdx.size()==1){
//        lines.push_back(restRayList[hIdx[0]]);
//        amps.push_back(0.058);//1608
//        inds.push_back(hIdx[0]);
//    }
//    if(lines.size()>0){
//        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, amps, m_nominalWidthDefaultAbsorption, inds)));
//    }

    //Load MgII multilines
    aIdx = findLineIdxInCatalog( restRayList, "MgII2803", CRay::nType_Absorption);
    bIdx = findLineIdxInCatalog( restRayList, "MgII2796", CRay::nType_Absorption);
    std::vector<CRay> linesMgII;
    std::vector<Float64> ampsMgII;
    std::vector<Int32> indsMgII;
    if(aIdx.size()==1){
        linesMgII.push_back(restRayList[aIdx[0]]);
        ampsMgII.push_back(0.3054);//MgII2803
        indsMgII.push_back(aIdx[0]);
    }
    if(bIdx.size()==1){
        linesMgII.push_back(restRayList[bIdx[0]]);
        ampsMgII.push_back(0.6123); //MgII2796
        indsMgII.push_back(bIdx[0]);
    }
    if(lines.size()>0){
        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesMgII, m_LineWidthType, ampsMgII, m_nominalWidthDefaultAbsorption, indsMgII)));
    }

//    //Load SiIV multilines
//    aIdx = findLineIdxInCatalog( restRayList, "SiIV1402", CRay::nType_Absorption);
//    bIdx = findLineIdxInCatalog( restRayList, "SiIV1393", CRay::nType_Absorption);
//    std::vector<CRay> linesSiIV;
//    std::vector<Float64> ampsSiIV;
//    std::vector<Int32> indsSiIV;
//    if(aIdx.size()==1){
//        linesSiIV.push_back(restRayList[aIdx[0]]);
//        ampsSiIV.push_back(0.262);//SiIV1402
//        indsSiIV.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        linesSiIV.push_back(restRayList[bIdx[0]]);
//        ampsSiIV.push_back(0.528); //SiIV1393
//        indsSiIV.push_back(bIdx[0]);
//    }
//    if(lines.size()>0){
//        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesSiIV, m_LineWidthType, ampsSiIV, m_nominalWidthDefaultAbsorption, indsSiIV)));
//    }

//    //Load CIV multilines, emission possible = warning
//    aIdx = findLineIdxInCatalog( restRayList, "CIV1550", CRay::nType_Absorption);
//    bIdx = findLineIdxInCatalog( restRayList, "CIV1548", CRay::nType_Absorption);
//    std::vector<CRay> linesCIV;
//    std::vector<Float64> ampsCIV;
//    std::vector<Int32> indsCIV;
//    if(aIdx.size()==1){
//        linesCIV.push_back(restRayList[aIdx[0]]);
//        ampsCIV.push_back(0.09522);//CIV1550
//        indsCIV.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        linesCIV.push_back(restRayList[bIdx[0]]);
//        ampsCIV.push_back(0.1908); //CIV1548
//        indsCIV.push_back(bIdx[0]);
//    }
//    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesCIV, m_LineWidthType, ampsCIV, m_nominalWidthDefaultAbsorption, indsCIV)));

//    //Load CaII multilines
//    aIdx = findLineIdxInCatalog( restRayList, "CaII_H", CRay::nType_Absorption);
//    bIdx = findLineIdxInCatalog( restRayList, "CaII_K", CRay::nType_Absorption);
//    std::vector<CRay> linesCaII;
//    std::vector<Float64> ampsCaII;
//    std::vector<Int32> indsCaII;
//    if(aIdx.size()==1){
//        linesCaII.push_back(restRayList[aIdx[0]]);
//        ampsCaII.push_back(1.0);//CaII_H
//        indsCaII.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        linesCaII.push_back(restRayList[bIdx[0]]);
//        ampsCaII.push_back(0.6); //CaII_K
//        indsCaII.push_back(bIdx[0]);
//    }
//    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesCaII, m_LineWidthType, ampsCaII, m_nominalWidthDefaultAbsorption, indsCaII)));


    //Load the rest of the single lines
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if ( FindElementIndex(iRestRay)==-1 )
        {
            Float64 nominalwidth = m_nominalWidthDefaultEmission;
            if(restRayList[iRestRay].GetType() == CRay::nType_Emission){
                nominalwidth = m_nominalWidthDefaultEmission;
            }else{
                nominalwidth = m_nominalWidthDefaultAbsorption;
            }
            addSingleLine(restRayList[iRestRay], iRestRay, nominalwidth);
        }
    }
}


void CLineModelElementList::LoadCatalog_tplExtendedBlue(const CRayCatalog::TRayVector& restRayList)
{
    std::vector<CRay> lines;
    std::vector<Float64> amps;
    std::vector<Int32> inds;

    //Load OIII multilines
    std::vector<Int32> OIIIaIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> OIIIbIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1/3)", CRay::nType_Emission);
    if(OIIIaIdx.size()==1){
        lines.push_back(restRayList[OIIIaIdx[0]]);
        amps.push_back(0.84);//[OIII](doublet-1)
        inds.push_back(OIIIaIdx[0]);
    }
    if(OIIIbIdx.size()==1){
        lines.push_back(restRayList[OIIIbIdx[0]]);
        amps.push_back(0.28);//[OIII](doublet-1/3)
        inds.push_back(OIIIbIdx[0]);
    }
    std::vector<Int32> Hbetadx = findLineIdxInCatalog( restRayList, "Hbeta", CRay::nType_Emission);
    if(Hbetadx.size()==1){
        lines.push_back(restRayList[Hbetadx[0]]);
        amps.push_back(0.338); //Hbetadx
        inds.push_back(Hbetadx[0]);
    }
//    if(lines.size()>0){
//        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, amps, m_nominalWidthDefaultEmission, inds)));
//    }
//    lines.clear();
//    amps.clear();
//    inds.clear();
    //Load OII line doublet

//    std::vector<Int32> OIIdx = findLineIdxInCatalog( restRayList, "[OII]3727", CRay::nType_Emission);
//    if(OIIdx.size()==1){
//        lines.push_back(restRayList[OIIdx[0]]);
//        amps.push_back(1.0); //[OII]3727
//        inds.push_back(OIIdx[0]);
//    }

    std::vector<Int32> OII1dx = findLineIdxInCatalog( restRayList, "[OII]3729", CRay::nType_Emission);
    std::vector<Int32> OII2dx = findLineIdxInCatalog( restRayList, "[OII]3726", CRay::nType_Emission);
    if(OII1dx.size()==1){
        lines.push_back(restRayList[OII1dx[0]]);
        amps.push_back(1.0); //[OII]3729
        inds.push_back(OII1dx[0]);
    }
    if(OII2dx.size()==1){
        lines.push_back(restRayList[OII2dx[0]]);
        amps.push_back(1.0); //[OII]3726
        inds.push_back(OII2dx[0]);
    }
    if(lines.size()>0){
        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, amps, m_nominalWidthDefaultEmission, inds)));
    }

    //Load NII multilines
    std::vector<Int32> NII1dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> NII2dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1/2.95)", CRay::nType_Emission);
    if(NII1dx.size()==1 && NII2dx.size()==1){
        addDoubleLine(restRayList[NII1dx[0]], restRayList[NII2dx[0]], NII1dx[0], NII2dx[0], m_nominalWidthDefaultEmission, 1.0, 1.0/2.95);
    }else{
        if(NII1dx.size()==1){
            addSingleLine(restRayList[NII1dx[0]], NII1dx[0], m_nominalWidthDefaultEmission);
        }
        if(NII2dx.size()==1){
            addSingleLine(restRayList[NII2dx[0]], NII2dx[0], m_nominalWidthDefaultEmission);
        }
    }


    //Load FeII 2586/2600 multilines
    std::vector<Int32> aIdx = findLineIdxInCatalog( restRayList, "FeII2600", CRay::nType_Absorption);
    std::vector<Int32> bIdx = findLineIdxInCatalog( restRayList, "FeII2586", CRay::nType_Absorption);
    std::vector<Int32> cIdx = findLineIdxInCatalog( restRayList, "FeII2382", CRay::nType_Absorption);
    std::vector<Int32> dIdx = findLineIdxInCatalog( restRayList, "FeII2374", CRay::nType_Absorption);
    std::vector<Int32> eIdx = findLineIdxInCatalog( restRayList, "FeII2344", CRay::nType_Absorption);
    std::vector<Int32> fIdx = findLineIdxInCatalog( restRayList, "FeII2260", CRay::nType_Absorption);
    std::vector<Int32> gIdx = findLineIdxInCatalog( restRayList, "FeII2249", CRay::nType_Absorption);
    std::vector<Int32> hIdx = findLineIdxInCatalog( restRayList, "FeII1608", CRay::nType_Absorption);
    lines.clear();
    amps.clear();
    inds.clear();
//    if(aIdx.size()==1){
//        lines.push_back(restRayList[aIdx[0]]);
//        amps.push_back(0.239);//2600
//        inds.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        lines.push_back(restRayList[bIdx[0]]);
//        amps.push_back(0.0691); //2586
//        inds.push_back(bIdx[0]);
//    }
//    if(cIdx.size()==1){
//        lines.push_back(restRayList[cIdx[0]]);
//        amps.push_back(0.32); //2382
//        inds.push_back(cIdx[0]);
//    }
//    if(dIdx.size()==1){
//        lines.push_back(restRayList[dIdx[0]]);
//        amps.push_back(0.0313); //2374
//        inds.push_back(dIdx[0]);
//    }
//    if(eIdx.size()==1){
//        lines.push_back(restRayList[eIdx[0]]);
//        amps.push_back(0.114); //2344
//        inds.push_back(eIdx[0]);
//    }
//    if(fIdx.size()==1){
//        lines.push_back(restRayList[fIdx[0]]);
//        amps.push_back(0.00244); //2260
//        inds.push_back(fIdx[0]);
//    }
//    if(gIdx.size()==1){
//        lines.push_back(restRayList[gIdx[0]]);
//        amps.push_back(0.001821);//2249
//        inds.push_back(gIdx[0]);
//    }
//    if(hIdx.size()==1){
//        lines.push_back(restRayList[hIdx[0]]);
//        amps.push_back(0.058);//1608
//        inds.push_back(hIdx[0]);
//    }
//    if(lines.size()>0){
//        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, amps, m_nominalWidthDefaultAbsorption, inds)));
//    }

//    //Load MgII multilines
//    aIdx = findLineIdxInCatalog( restRayList, "MgII2803", CRay::nType_Absorption);
//    bIdx = findLineIdxInCatalog( restRayList, "MgII2796", CRay::nType_Absorption);
//    std::vector<CRay> linesMgII;
//    std::vector<Float64> ampsMgII;
//    std::vector<Int32> indsMgII;
//    if(aIdx.size()==1){
//        linesMgII.push_back(restRayList[aIdx[0]]);
//        ampsMgII.push_back(0.3054);//MgII2803
//        indsMgII.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        linesMgII.push_back(restRayList[bIdx[0]]);
//        ampsMgII.push_back(0.6123); //MgII2796
//        indsMgII.push_back(bIdx[0]);
//    }
//    if(lines.size()>0){
//        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesMgII, m_LineWidthType, ampsMgII, m_nominalWidthDefaultAbsorption, indsMgII)));
//    }

    //Load SiIV multilines
    aIdx = findLineIdxInCatalog( restRayList, "SiIV1402", CRay::nType_Absorption);
    bIdx = findLineIdxInCatalog( restRayList, "SiIV1393", CRay::nType_Absorption);
    std::vector<CRay> linesSiIV;
    std::vector<Float64> ampsSiIV;
    std::vector<Int32> indsSiIV;
    if(aIdx.size()==1){
        linesSiIV.push_back(restRayList[aIdx[0]]);
        ampsSiIV.push_back(0.262);//SiIV1402
        indsSiIV.push_back(aIdx[0]);
    }
    if(bIdx.size()==1){
        linesSiIV.push_back(restRayList[bIdx[0]]);
        ampsSiIV.push_back(0.528); //SiIV1393
        indsSiIV.push_back(bIdx[0]);
    }
    if(lines.size()>0){
        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesSiIV, m_LineWidthType, ampsSiIV, m_nominalWidthDefaultAbsorption, indsSiIV)));
    }

//    //Load CIV multilines, emission possible = warning
//    aIdx = findLineIdxInCatalog( restRayList, "CIV1550", CRay::nType_Absorption);
//    bIdx = findLineIdxInCatalog( restRayList, "CIV1548", CRay::nType_Absorption);
//    std::vector<CRay> linesCIV;
//    std::vector<Float64> ampsCIV;
//    std::vector<Int32> indsCIV;
//    if(aIdx.size()==1){
//        linesCIV.push_back(restRayList[aIdx[0]]);
//        ampsCIV.push_back(0.09522);//CIV1550
//        indsCIV.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        linesCIV.push_back(restRayList[bIdx[0]]);
//        ampsCIV.push_back(0.1908); //CIV1548
//        indsCIV.push_back(bIdx[0]);
//    }
//    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesCIV, m_LineWidthType, ampsCIV, m_nominalWidthDefaultAbsorption, indsCIV)));

//    //Load CaII multilines
//    aIdx = findLineIdxInCatalog( restRayList, "CaII_H", CRay::nType_Absorption);
//    bIdx = findLineIdxInCatalog( restRayList, "CaII_K", CRay::nType_Absorption);
//    std::vector<CRay> linesCaII;
//    std::vector<Float64> ampsCaII;
//    std::vector<Int32> indsCaII;
//    if(aIdx.size()==1){
//        linesCaII.push_back(restRayList[aIdx[0]]);
//        ampsCaII.push_back(1.0);//CaII_H
//        indsCaII.push_back(aIdx[0]);
//    }
//    if(bIdx.size()==1){
//        linesCaII.push_back(restRayList[bIdx[0]]);
//        ampsCaII.push_back(0.6); //CaII_K
//        indsCaII.push_back(bIdx[0]);
//    }
//    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(linesCaII, m_LineWidthType, ampsCaII, m_nominalWidthDefaultAbsorption, indsCaII)));


    //Load the rest of the single lines
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if ( FindElementIndex(iRestRay)==-1 )
        {
            Float64 nominalwidth = m_nominalWidthDefaultEmission;
            if(restRayList[iRestRay].GetType() == CRay::nType_Emission){
                nominalwidth = m_nominalWidthDefaultEmission;
            }else{
                nominalwidth = m_nominalWidthDefaultAbsorption;
            }
            addSingleLine(restRayList[iRestRay], iRestRay, nominalwidth);
        }
    }
}

void CLineModelElementList::LoadCatalogMultilineBalmer(const CRayCatalog::TRayVector& restRayList)
{
    //Load OIII multilines
    std::vector<Int32> OIIIaIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> OIIIbIdx = findLineIdxInCatalog( restRayList, "[OIII](doublet-1/3)", CRay::nType_Emission);
    if(OIIIaIdx.size()==1 && OIIIbIdx.size()==1){
        addDoubleLine(restRayList[OIIIaIdx[0]], restRayList[OIIIbIdx[0]], OIIIaIdx[0], OIIIbIdx[0], m_nominalWidthDefaultEmission, 1.0, 1.0/3.0);
    }else{
        if(OIIIaIdx.size()==1){
            addSingleLine(restRayList[OIIIaIdx[0]], OIIIaIdx[0], m_nominalWidthDefaultEmission);
        }
        if(OIIIbIdx.size()==1){
            addSingleLine(restRayList[OIIIbIdx[0]], OIIIbIdx[0], m_nominalWidthDefaultEmission);
        }
    }

    //Load NII multilines
    std::vector<Int32> NII1dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> NII2dx = findLineIdxInCatalog( restRayList, "[NII](doublet-1/2.95)", CRay::nType_Emission);
    if(NII1dx.size()==1 && NII2dx.size()==1){
        addDoubleLine(restRayList[NII1dx[0]], restRayList[NII2dx[0]], NII1dx[0], NII2dx[0], m_nominalWidthDefaultEmission, 1.0, 1.0/2.95);
    }else{
        if(NII1dx.size()==1){
            addSingleLine(restRayList[NII1dx[0]], NII1dx[0], m_nominalWidthDefaultEmission);
        }
        if(NII2dx.size()==1){
            addSingleLine(restRayList[NII2dx[0]], NII2dx[0], m_nominalWidthDefaultEmission);
        }
    }

    //Load OII line doublet
    std::vector<Int32> OII1dx = findLineIdxInCatalog( restRayList, "[OII](doublet-1)", CRay::nType_Emission);
    std::vector<Int32> OII2dx = findLineIdxInCatalog( restRayList, "[OII](doublet-1/3)", CRay::nType_Emission);
    if(OII1dx.size()==1 && OII2dx.size()==1){
        addDoubleLine(restRayList[OII1dx[0]], restRayList[OII2dx[0]], OII1dx[0], OII2dx[0], 3.2, 1.0, 1.0/2.9);
    }else{
        if(OII1dx.size()==1){
            addSingleLine(restRayList[OII1dx[0]], OII1dx[0], m_nominalWidthDefaultEmission);
        }
        if(OII2dx.size()==1){
            addSingleLine(restRayList[OII2dx[0]], OII2dx[0], m_nominalWidthDefaultEmission);
        }
    }

    //Load Balmer multilines
    std::vector<Int32> Halphaidx = findLineIdxInCatalog( restRayList, "Halpha", CRay::nType_Emission);
    std::vector<Int32> Hbetaidx = findLineIdxInCatalog( restRayList, "Hbeta", CRay::nType_Emission);
    std::vector<Int32> Hgammaidx = findLineIdxInCatalog( restRayList, "Hgamma", CRay::nType_Emission);
    std::vector<Int32> Hdeltaidx = findLineIdxInCatalog( restRayList, "Hdelta", CRay::nType_Emission);
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

        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, amps, m_nominalWidthDefaultEmission, inds)));
    }else{
        if(Halphaidx.size()==1){
            addSingleLine(restRayList[Halphaidx[0]], Halphaidx[0], m_nominalWidthDefaultEmission);
        }
        if(Hbetaidx.size()==1){
            addSingleLine(restRayList[Hbetaidx[0]], Hbetaidx[0], m_nominalWidthDefaultEmission);
        }
        if(Hgammaidx.size()==1){
            addSingleLine(restRayList[Hgammaidx[0]], Hgammaidx[0], m_nominalWidthDefaultEmission);
        }
        if(Hdeltaidx.size()==1){
            addSingleLine(restRayList[Hdeltaidx[0]], Hdeltaidx[0], m_nominalWidthDefaultEmission);
        }
    }

    //Load the rest of the single lines
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if ( FindElementIndex(iRestRay)==-1 )
        {
            addSingleLine(restRayList[iRestRay], iRestRay, m_nominalWidthDefaultEmission);
        }
    }
}

void CLineModelElementList::LogCatalogInfos()
{
    Log.LogInfo( "\n");
    Log.LogInfo( "LineModel Infos: %d elements", m_Elements.size());
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        Int32 nRays = m_Elements[iElts]->GetSize();
        if(nRays<1)
        {
            Log.LogInfo( "LineModel ctlg: elt %d (%s): no rays", iElts, m_Elements[iElts]->GetElementTypeTag().c_str());

        }
        for(UInt32 j=0; j<nRays; j++){
            std::string nominalAmpStr = "";
            if(nRays>1){
                nominalAmpStr = boost::str(boost::format("(nominal amp = %.4f)") % m_Elements[iElts]->GetNominalAmplitude(j));
            }
            Log.LogInfo( "LineModel ctlg: elt %d (%s): ray %d = %s %s", iElts, m_Elements[iElts]->GetElementTypeTag().c_str(), j, m_Elements[iElts]->GetRayName(j).c_str(), nominalAmpStr.c_str());
        }
    }
    Log.LogInfo( "\n");
}



void CLineModelElementList::LoadCatalogSingleLines(const CRayCatalog::TRayVector& restRayList)
{
    //Load the rest of the single lines
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if ( FindElementIndex(iRestRay)==-1 )
        {
            addSingleLine(restRayList[iRestRay], iRestRay, m_nominalWidthDefaultEmission);
        }
    }
}

void CLineModelElementList::LoadContinuum()
{
    //path p( templatePath );
    //path name = p.leaf();
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/emission/NEW_Im_extended.dat";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/galaxy/EW_SB2extended.dat";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/galaxy/BulgedataExtensionData.dat";

    std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission/NEW_Im_extended_blue_continuum.txt";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission/NEW_Im_extended_continuum.txt";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/galaxy/EW_SB2extended.txt";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/galaxy/BulgedataExtensionData.txt";

    //std::string templatePath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/47002690000013_flam.txt";


    //CRef<CTemplate> tmpl = new CTemplate();
    CTemplate tpl;
    CSpectrumIOGenericReader asciiReader;
    if( !asciiReader.Read( templatePath.c_str(), tpl ) ) {
        Log.LogError("Fail to read template: %s", templatePath);
        return;
    }

//    //save as a spectrum
//    CSpectrumIOFitsWriter writer;
//    Bool retVal1 = writer.Write( "tplexported.fits",  tpl);


//    //create a continuum tpl
//    if(1)
//    {
//        // Remove continuum
//        CContinuumIrregularSamplingMedian continuum;
//        CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;

//        Int32 retVal = continuum.RemoveContinuum( tpl, fluxAxisWithoutContinuumCalc );
//        CSpectrumFluxAxis fluxAxis = tpl.GetFluxAxis();
//        fluxAxis.Subtract(fluxAxisWithoutContinuumCalc);
//        CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
//        //*//debug:
//        // save continuum tpl and xmad,  flux data
//        FILE* f = fopen( "continuum_tpl_dbg.txt", "w+" );
//        for( Int32 t=0;t<fluxAxisWithoutContinuumCalc.GetSamplesCount();t++)
//        {
//            fprintf( f, "%f %f\n", t, tplSpectralAxis[t], fluxAxis[t]);//*1e12);
//        }
//        fclose( f );
//        //*/
//        return;
//    }


    //*/
    // Precalculate a fine grid template to be used for the 'closest value' rebin method
    Int32 n = tpl.GetSampleCount();
    CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
    CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
    //Float64 dLambdaTgt =  1.0 * ( spectrum.GetMeanResolution()*0.9 )/( 1+sortedRedshifts[sortedRedshifts.size()-1] );
    Float64 dLambdaTgt =  0.1; //should be sufficient for the continuum
    //Float64 lmin = tplSpectralAxis[0];
    Float64 lmin = 0;
    Float64 lmax = tplSpectralAxis[n-1];
    Int32 nTgt = (lmax-lmin)/dLambdaTgt + 2.0/dLambdaTgt;

    m_precomputedFineGridContinuumFlux = (Float64*)malloc(nTgt*sizeof(Float64));

    //inialise and allocate the gsl objects
    Float64* Ysrc = tplFluxAxis.GetSamples();
    Float64* Xsrc = tplSpectralAxis.GetSamples();

    //spline
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init (spline, Xsrc, Ysrc, n);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    Int32 k = 0;
    Float64 x = 0.0;
    for(k=0; k<nTgt; k++){
        x = lmin + k*dLambdaTgt;
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            m_precomputedFineGridContinuumFlux[k] = 0.0; //todo, make sure this is never used in the next steps...
        }else{
            m_precomputedFineGridContinuumFlux[k] = gsl_spline_eval (spline, x, accelerator);
        }
    }
    //*/

}

void CLineModelElementList::PrepareContinuum(Float64 z)
{    
    const CSpectrumSpectralAxis& targetSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const Float64* Xtgt = targetSpectralAxis.GetSamples();
    Float64* Yrebin = m_ContinuumFluxAxis.GetSamples();

    if(m_precomputedFineGridContinuumFlux == NULL){
        for ( Int32 i = 0; i<targetSpectralAxis.GetSamplesCount(); i++)
        {
            //Yrebin[i] = m_SpcNoContinuumFluxAxis[i]; //
            Yrebin[i] = 0.0;
        }
        return;
    }

    TFloat64Range currentRange(Xtgt[0], Xtgt[targetSpectralAxis.GetSamplesCount()-1]);

    // Move cursors up to lambda range start
    Int32 j = 0;
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] < currentRange.GetBegin() )
    {
        Yrebin[j] = 0.0;
        j++;
    }
    //* // Precomputed FINE GRID nearest sample,
    Int32 k = 0;
    Float64 dl = 0.1;
    Float64 Coeffk = 1.0/dl/(1+z);
    // For each sample in the target spectrum
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
    {
        k = (int)(Xtgt[j]*Coeffk+0.5);
        Yrebin[j] = m_precomputedFineGridContinuumFlux[k];
        //Yrebin[j] = 0.0;
        j++;

    }
    //*/
    while( j < targetSpectralAxis.GetSamplesCount() )
    {
        Yrebin[j] = 0.0;
        j++;
    }


    //fit the continuum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    //const Float64* flux = m_SpcFluxAxis.GetSamples();
    const Float64* flux = m_SpcNoContinuumFluxAxis.GetSamples();

    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = m_SpcFluxAxis.GetError();

    Float64 sumCross = 0.0;
    Float64 sumGauss = 0.0;
    Float64 err2 = 0.0;
    Int32 num = 0;

    Float64 x=0.0;
    Float64 y=0.0;
    Float64 yg=0.0;

    //A estimation
    for ( Int32 i = 0; i<targetSpectralAxis.GetSamplesCount(); i++)
    {
        y = flux[i];
        x = spectral[i];
        yg = Yrebin[i];

        num++;
        err2 = 1.0 / (error[i] * error[i]);
        sumCross += yg*y*err2;
        sumGauss += yg*yg*err2;
    }

    if ( num==0 || sumGauss==0 )
    {
        return;
    }

    Float64 A = std::max(0.0, sumCross / sumGauss);
    for ( Int32 i = 0; i<targetSpectralAxis.GetSamplesCount(); i++)
    {
        Yrebin[i] *=A;
    }
}

void CLineModelElementList::fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelResult::SLineModelSolution& modelSolution)
{
    m_Redshift = redshift;

    //prepare the continuum
    PrepareContinuum(redshift);

    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    CSpectrumFluxAxis spcFluxAxisNoContinuum(m_SpcFluxAxis);

    const Float64* error = m_SpcFluxAxis.GetError();
    Float64* errorNoContinuum = spcFluxAxisNoContinuum.GetError();

    spcFluxAxisNoContinuum.SetSize( modelFluxAxis.GetSamplesCount() );

    for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
        modelFluxAxis[i] = m_ContinuumFluxAxis[i];
        spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
        errorNoContinuum[i] = error[i];
    }

    //prepare the elements
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->prepareSupport(spectralAxis, redshift, lambdaRange);
    }


    //fit the amplitudes of each element independently
    if(0)
    {
        //fit the model amplitudes individually
        for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
        {
            m_Elements[iElts]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
        }
    }
    //else{
    //    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    //    {
    //        m_Elements[iElts]->SetFittedAmplitude(1e-19);
    //    }
    //}

    //fit the amplitude of all elements together with iterative solver: Nelder Mead Simplex
    if(0){
        //fit the amplitudes together
        fitAmplitudesSimplex();
    }
    //fit the amplitude of all elements together with linear solver: gsl_multifit_wlinear
    if(0){
        std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
        fitAmplitudesLinSolve(validEltsIdx);
    }

    //fit the amplitudes of each element independently, unless there is overlap
    if(1)
    {
        std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();

        std::vector<Int32> indexesFitted;
        for( UInt32 iValidElts=0; iValidElts<validEltsIdx.size(); iValidElts++ )
        {
            Int32 iElts = validEltsIdx[iValidElts];

            //skip if already fitted
            for(Int32 i=0; i<indexesFitted.size(); i++){
                if(iElts == indexesFitted[i]){
                    continue;
                }
            }

            //do the fit on the ovelapping elements
            std::vector<Int32> overlappingInds = getOverlappingElements(iElts);
            if(overlappingInds.size()<2){
                m_Elements[iElts]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
            }else{
                fitAmplitudesLinSolve(overlappingInds);
            }

            //update the already fitted list
            for(Int32 i=0; i<overlappingInds.size(); i++){
                indexesFitted.push_back(overlappingInds[i]);
            }

        }
    }


    //Apply rules,
    applyRules();

    refreshModel();

    //create spectrum model
    modelSolution = GetModelSolution();

}

void CLineModelElementList::refreshModel()
{

    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //create spectrum model
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis);
    }
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_Redshift);
    }

}

void CLineModelElementList::fitAmplitudesSimplex()
{
    CMultiGaussianFit fitter;

    Int32 status = fitter.Compute( *this );
//    if(status!=NSEpic::CMultiGaussianFit::nStatus_Success){
//        continue;
//    }

//    Float64 gaussAmp;
//    Float64 gaussPos;
//    Float64 gaussWidth;
//    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
//    Float64 gaussAmpErr;
//    Float64 gaussPosErr;
//    Float64 gaussWidthErr;
//    fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );


}

std::vector<Int32> CLineModelElementList::getSupportIndexes( std::vector<Int32> EltsIdx)
{
    std::vector<Int32> indexes;

    TInt32RangeList support;
    for( UInt32 i=0; i<EltsIdx.size(); i++ )
    {
        Int32 iElts = EltsIdx[i];

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


std::vector<Int32> CLineModelElementList::getOverlappingElements(  Int32 ind )
{
    TInt32RangeList refsupport = m_Elements[ind]->getSupport();
    Int32 linetype = m_RestRayList[m_Elements[ind]->m_LineCatalogIndexes[0]].GetType();


    std::vector<Int32> indexes;

    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetype){
            continue;
        }

        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            for( UInt32 iRefS=0; iRefS<refsupport.size(); iRefS++ )
            {
                Int32 x1, y1, x2, y2;
                x1 = refsupport[iRefS].GetBegin();
                x2 = refsupport[iRefS].GetEnd();
                y1 = s[iS].GetBegin();
                y2 = s[iS].GetEnd();

                if( std::max(x1,y1) < std::min(x2,y2) ){
                    indexes.push_back(iElts);
                    continue;
                }
            }
        }
    }

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    return indexes;
}

void CLineModelElementList::fitAmplitudesLinSolve( std::vector<Int32> EltsIdx )
{
    boost::chrono::thread_clock::time_point start_prep = boost::chrono::thread_clock::now();

    Int32 nddl = EltsIdx.size();
    if(nddl<1){
        return;
    }
    std::vector<Int32> xInds = getSupportIndexes( EltsIdx );

    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        SetElementAmplitude(EltsIdx[iddl], 1.0);
    }


    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    const Float64* spectral = spectralAxis.GetSamples();
    Float64* flux = m_SpcFluxAxis.GetSamples();
    const Float64* error = m_SpcFluxAxis.GetError();


    //Linear fit
    int i, n;
    Float64 fval;
    double xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = xInds.size();

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);

    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    //
    //todo: normalize, center...
    //
    Int32 idx = 0;
    for (i = 0; i < n; i++)
    {
        idx = xInds[i];
        xi = spectral[idx];
        yi = flux[idx];
        ei = error[idx];

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            fval =  m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, m_Redshift);
            gsl_matrix_set (X, i, iddl, fval);
        }

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    //
    boost::chrono::thread_clock::time_point stop_prep = boost::chrono::thread_clock::now();
    Float64 duration_prep = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep - start_prep).count();
    boost::chrono::thread_clock::time_point start_fit = boost::chrono::thread_clock::now();

    {
      gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
      gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
      gsl_multifit_linear_free (work);
    }

    //
    boost::chrono::thread_clock::time_point stop_fit = boost::chrono::thread_clock::now();
    Float64 duration_fit = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit - start_fit).count();
    //Log.LogInfo( "LineModel linear fit: prep = %.3f - fit = %.3f", duration_prep, duration_fit);

//#define C(i) (gsl_vector_get(c,(i)))
//#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        Float64 a = gsl_vector_get(c,iddl);
        SetElementAmplitude(EltsIdx[iddl], a);
    }
    //refreshModel();

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
}


Float64 CLineModelElementList::getLeastSquareMerit(const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 fit = 0.0;
    const Float64* error = spcFluxAxis.GetError();
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());

    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        // fit
        diff = (Yspc[j] - Ymodel[j]);
        fit += (diff*diff) / (error[j]*error[j]);
        //fit += (diff*diff) / (1e-16*1e-16);
        //fit += (diff*diff)/ (error[0]*error[0]);
        //fit += pow( Yspc[j] - Ymodel[j] , 2.0 );
    }
    //fit /= numDevs;

    return fit;
}

Float64 CLineModelElementList::getLeastSquareMeritUnderElements()
{
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 fit = 0;
    const Float64* error = spcFluxAxis.GetError();
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    TInt32RangeList support;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
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
            numDevs++;
            // fit
            diff = (Yspc[j] - Ymodel[j]);
            fit += (diff*diff) / (error[j]*error[j]);
            //fit += diff*diff;

        }
    }


    return fit;
}

std::vector<int> CLineModelElementList::findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag, Int32 type)
{
    std::vector<Int32> indexes;
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if(restRayList[iRestRay].GetType() != type){
            continue;
        }
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
    //CSingleLine line = CSingleLine(r, m_LineWidthType, nominalWidth);
    std::vector<Int32> a;
    a.push_back(index);
    //CSingleLine c(r, nominalWidth, a);
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CSingleLine(r, m_LineWidthType, nominalWidth, a)));
    //m_Elements.push_back(new CSingleLine(r, m_LineWidthType, nominalWidth, a));
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
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, amps, nominalWidth, a)));
    //m_Elements.push_back(new CSingleLine(r, nominalWidth, a));
}

void CLineModelElementList::applyRules()
{
    //todo: check for the noise when applying amplitudes rules...

    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Halpha", "Hbeta", 1.0/2.86);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hbeta", "Hgamma", 0.47);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hgamma", "Hdelta", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hdelta", "Hepsilon", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hepsilon", "H8", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "H8", "H9", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "H9", "H10", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "H10", "H11", 1.0);

    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HbetaA", "HgammaA", 0.47);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HgammaA", "HdeltaA", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HdeltaA", "HepsilonA", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HepsilonA", "H8A", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "H8A", "H9A", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "H9A", "H10A", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "H10A", "H11A", 1.0);

    ApplyStrongHigherWeakRule(CRay::nType_Emission);
    ApplyStrongHigherWeakRule(CRay::nType_Absorption);

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

Void CLineModelElementList::ApplyStrongHigherWeakRule( Int32 linetype )
{
    Float64 coeff = 0.5;
    Float64 erStrong=-1.0;
    Float64 maxiStrong = FindHighestStrongLineAmp(linetype, erStrong);


    for( UInt32 iRestRayWeak=0; iRestRayWeak<m_RestRayList.size(); iRestRayWeak++ ) //loop on the strong lines
    {
        if(m_RestRayList[iRestRayWeak].GetForce() != CRay::nForce_Weak){
            continue;
        }
        if(m_RestRayList[iRestRayWeak].GetType() != linetype){
            continue;
        }
        Int32 eIdxWeak = FindElementIndex(iRestRayWeak);
        Int32 subeIdxWeak = m_Elements[eIdxWeak]->FindElementIndex(iRestRayWeak);

        if(m_Elements[eIdxWeak]->IsOutsideLambdaRange(subeIdxWeak) == true){
            continue;
        }

        Float64 nSigma = 1.0;
        Float64 ampA = maxiStrong;
        Float64 erA = erStrong;

        Float64 ampB = m_Elements[eIdxWeak]->GetFittedAmplitude(subeIdxWeak);
        Float64 erB = m_Elements[eIdxWeak]->GetFittedAmplitudeErrorSigma(subeIdxWeak);

        Float64 maxB = (coeff*ampA) + coeff*(erA*nSigma);

        m_Elements[eIdxWeak]->LimitFittedAmplitude(subeIdxWeak, maxB);

    }

}

Float64 CLineModelElementList::FindHighestStrongLineAmp( Int32 linetype , Float64 &er)
{
    Float64 maxi = -1.0;
    for( UInt32 iRestRayStrong=0; iRestRayStrong<m_RestRayList.size(); iRestRayStrong++ ) //loop on the strong lines
    {
        if(m_RestRayList[iRestRayStrong].GetForce() != CRay::nForce_Strong){
            continue;
        }
        if(m_RestRayList[iRestRayStrong].GetType() != linetype){
            continue;
        }

        Int32 eIdxStrong = FindElementIndex(iRestRayStrong);
        Int32 subeIdxStrong = m_Elements[eIdxStrong]->FindElementIndex(iRestRayStrong);

        if(m_Elements[eIdxStrong]->IsOutsideLambdaRange(subeIdxStrong) == true){
            continue;
        }


        Float64 ampStrong = m_Elements[eIdxStrong]->GetFittedAmplitude(subeIdxStrong);
        if(maxi<ampStrong){
            maxi = ampStrong;
            er = m_Elements[eIdxStrong]->GetFittedAmplitudeErrorSigma(subeIdxStrong);
        }
    }
    return maxi;
}

Void CLineModelElementList::Apply2SingleLinesAmplitudeRule( Int32 linetype, std::string lineA, std::string lineB, Float64 coeff )
{
    Int32 iA = FindElementIndex(lineA, linetype);
    if(iA==-1){
        return;
    }
    if(m_Elements[iA]->GetSize()>1){
        iA=-1;
    }
    Int32 iB = FindElementIndex(lineB, linetype);
    if(iB==-1){
        return;
    }
    if(m_Elements[iB]->GetSize()>1){
        iB=-1;
    }
    if(iA==-1 || iB==-1 || iA==iB){
        return;
    }

    if(m_Elements[iA]->IsOutsideLambdaRange() == false){
        Float64 nSigma = 1.0;
        Float64 ampA = m_Elements[iA]->GetFittedAmplitude(0);
        Float64 erA = m_Elements[iA]->GetFittedAmplitudeErrorSigma(0);

        Float64 ampB = m_Elements[iB]->GetFittedAmplitude(0);
        Float64 erB = m_Elements[iB]->GetFittedAmplitudeErrorSigma(0);

        Float64 maxB = (coeff*ampA) + (erA*nSigma);

        m_Elements[iB]->LimitFittedAmplitude(0, maxB);
    }
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
        modelSolution.Errors.push_back(m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(subeIdx));
        //modelSolution.Widths.push_back(-1.0);
        //modelSolution.OutsideLambdaRange.push_back(true);
    }

    return modelSolution;
}

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

std::vector<Int32> CLineModelElementList::GetModelValidElementsIndexes()
{
    std::vector<Int32> nonZeroIndexes;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }

        nonZeroIndexes.push_back(iElts);
    }
    return nonZeroIndexes;
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

Int32 CLineModelElementList::FindElementIndex(std::string LineTagStr, Int32 linetype)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetype){
            continue;
        }

        if(m_Elements[iElts]->FindElementIndex(LineTagStr) !=-1){
            idx = iElts;
            break;
        }
    }
    return idx;
}


void CLineModelElementList::SetElementAmplitude(Int32 j, Float64 a)
{
    if(j>=0 && j<m_Elements.size())
    {
        m_Elements[j]->SetFittedAmplitude(a);
    }
    return;
}

Float64 CLineModelElementList::GetElementAmplitude(Int32 j)
{
    Float64 a=-1.0;
    if(j>=0 && j<m_Elements.size())
    {
        a = m_Elements[j]->GetElementAmplitude();
    }
    return a;
}


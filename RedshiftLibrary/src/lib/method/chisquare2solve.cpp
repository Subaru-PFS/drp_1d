#include <RedshiftLibrary/method/chisquare2solve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <float.h>
using namespace NSEpic;
using namespace std;

CMethodChisquare2Solve::CMethodChisquare2Solve( std::string calibrationPath )
{
    m_chiSquareOperator = new COperatorChiSquare2( calibrationPath );
}

CMethodChisquare2Solve::~CMethodChisquare2Solve()
{
    delete m_chiSquareOperator;
}


const std::string CMethodChisquare2Solve::GetDescription()
{
    std::string desc;

    desc = "Method chisquare2solve:\n";

    desc.append("\tparam: chisquare2solve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
    desc.append("\tparam: chisquare2solve.overlapThreshold = <float value>\n");
    desc.append("\tparam: chisquare2solve.interpolation = {""precomputedfinegrid"", ""lin""}\n");
    desc.append("\tparam: chisquare2solve.extinction = {""yes"", ""no""}\n");
    desc.append("\tparam: chisquare2solve.dustfit = {""yes"", ""no""}\n");


    return desc;

}


std::shared_ptr<CChisquare2SolveResult> CMethodChisquare2Solve::Compute(CDataStore& resultStore,
                                                                              const CSpectrum& spc,
                                                                              const CSpectrum& spcWithoutCont,
                                                                              const CTemplateCatalog& tplCatalog,
                                                                              const TStringList& tplCategoryList,
                                                                              const TFloat64Range& lambdaRange,
                                                                              const TFloat64List& redshifts,
                                                                              Float64 overlapThreshold,
                                                                              std::vector<CMask> maskList,
                                                                              std::string spcComponent,
                                                                              std::string opt_interp,
                                                                              std::string opt_extinction,
                                                                              std::string opt_dustFit)
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "chisquare2solve" );
    std::string scopeStr = "chisquare";

    Int32 _type;
    if(spcComponent=="raw"){
       _type = CChisquare2SolveResult::nType_raw;
       scopeStr = "chisquare";
    }else if(spcComponent=="nocontinuum"){
       _type = CChisquare2SolveResult::nType_noContinuum;
       scopeStr = "chisquare_nocontinuum";
    }else if(spcComponent=="continuum"){
        _type = CChisquare2SolveResult::nType_continuumOnly;
        scopeStr = "chisquare_continuum";
    }else if(spcComponent=="all"){
        _type = CChisquare2SolveResult::nType_all;
    }

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            Solve( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshifts, overlapThreshold, maskList, _type, opt_interp, opt_extinction, opt_dustFit);

            storeResult = true;
        }
    }


    if( storeResult )
    {
        std::shared_ptr< CChisquare2SolveResult>  ChisquareSolveResult = std::shared_ptr< CChisquare2SolveResult>( new CChisquare2SolveResult() );
        ChisquareSolveResult->m_type = _type;

        //std::string opt_combinePdf = "marg";
        std::string opt_combinePdf = "bestchi2";
        CombinePDF(resultStore, scopeStr, opt_combinePdf);

        return ChisquareSolveResult;
    }

    return NULL;
}

Bool CMethodChisquare2Solve::Solve(CDataStore& resultStore,
                                    const CSpectrum& spc,
                                    const CSpectrum& spcWithoutCont,
                                    const CTemplate& tpl,
                                    const CTemplate& tplWithoutCont,
                                    const TFloat64Range& lambdaRange,
                                    const TFloat64List& redshifts,
                                    Float64 overlapThreshold,
                                    std::vector<CMask> maskList,
                                    Int32 spctype,
                                    std::string opt_interp,
                                   std::string opt_extinction,
                                   std::string opt_dustFitting )
{
    CSpectrum _spc;
    CTemplate _tpl;
    std::string scopeStr = "chisquare";
    Int32 _ntype = 1;
    Int32 _spctype = spctype;
    Int32 _spctypetab[3] = {CChisquare2SolveResult::nType_raw, CChisquare2SolveResult::nType_noContinuum, CChisquare2SolveResult::nType_continuumOnly};


    Int32 enable_extinction = 0; //TODO: extinction should be deactivated for nocontinuum anyway ? TBD
    if(opt_extinction=="yes")
    {
        enable_extinction = 1;
    }
    Int32 enable_dustFitting = 0;
    if(opt_dustFitting=="yes")
    {
        enable_dustFitting = 1;
    }



    //case: nType_all
    if(spctype == CChisquare2SolveResult::nType_all){
        _ntype = 3;
    }

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == CChisquare2SolveResult::nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = spctype;
        }

        if(_spctype == CChisquare2SolveResult::nType_continuumOnly){
            // use continuum only
            _spc = spc;
            CSpectrumFluxAxis spcfluxAxis = _spc.GetFluxAxis();
            spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
            CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
            sfluxAxisPtr = spcfluxAxis;
            _tpl = tpl;
            CSpectrumFluxAxis tplfluxAxis = _tpl.GetFluxAxis();
            tplfluxAxis.Subtract(tplWithoutCont.GetFluxAxis());
            CSpectrumFluxAxis& tfluxAxisPtr = _tpl.GetFluxAxis();
            tfluxAxisPtr = tplfluxAxis;


            scopeStr = "chisquare_continuum";
        }else if(_spctype == CChisquare2SolveResult::nType_raw){
            // use full spectrum
            _spc = spc;
            _tpl = tpl;
            scopeStr = "chisquare";

        }else if(_spctype == CChisquare2SolveResult::nType_noContinuum){
            // use spectrum without continuum
            _spc = spc;
            CSpectrumFluxAxis spcfluxAxis = spcWithoutCont.GetFluxAxis();
            CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
            sfluxAxisPtr = spcfluxAxis;
            _tpl = tpl;
            CSpectrumFluxAxis tplfluxAxis = tplWithoutCont.GetFluxAxis();
            CSpectrumFluxAxis& tfluxAxisPtr = _tpl.GetFluxAxis();
            tfluxAxisPtr = tplfluxAxis;
            scopeStr = "chisquare_nocontinuum";
            //
            enable_dustFitting = 0;
        }

        // Compute merit function
        //CRef<CChisquareResult>  chisquareResult = (CChisquareResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
        auto  chisquareResult = std::dynamic_pointer_cast<CChisquareResult>( m_chiSquareOperator->Compute( _spc, _tpl, lambdaRange, redshifts, overlapThreshold, maskList, opt_interp, enable_extinction, enable_dustFitting ) );

        if( !chisquareResult )
        {
            //Log.LogInfo( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            resultStore.StoreScopedPerTemplateResult( tpl, scopeStr.c_str(), chisquareResult );

            //Save intermediate chisquare results
            if(m_opt_enableSaveIntermediateChisquareResults && chisquareResult->ChiSquareIntermediate.size()>0 && chisquareResult->ChiSquareIntermediate.size()==chisquareResult->Redshifts.size())
            {
                Int32 nISM = chisquareResult->ChiSquareIntermediate[0].size();
                if(chisquareResult->ChiSquareIntermediate[0].size()>0)
                {
                    Int32 nIGM = chisquareResult->ChiSquareIntermediate[0][0].size();

                    for(Int32 kism=0; kism<nISM; kism++)
                    {
                        for(Int32 kigm=0; kigm<nIGM; kigm++)
                        {
                            std::shared_ptr<CChisquareResult> result_chisquare_intermediate = std::shared_ptr<CChisquareResult>( new CChisquareResult() );
                            result_chisquare_intermediate->Init( chisquareResult->Redshifts.size(), 0, 0);
                            for(Int32 kz=0; kz<chisquareResult->Redshifts.size(); kz++)
                            {
                                result_chisquare_intermediate->Redshifts[kz] = chisquareResult->Redshifts[kz];
                                result_chisquare_intermediate->ChiSquare[kz] = chisquareResult->ChiSquareIntermediate[kz][kism][kigm];
                            }

                            std::string resname = (boost::format("%s_intermediate_ism%d_igm%d") % scopeStr.c_str() % kism % kigm).str();
                            resultStore.StoreScopedPerTemplateResult( tpl, resname.c_str(), result_chisquare_intermediate );
                        }
                    }
                }
            }
        }
    }

    return true;
}


Int32 CMethodChisquare2Solve::CombinePDF(CDataStore &store, std::string scopeStr, std::string opt_combine )
{
    Log.LogInfo("chisquare2solve: Pdfz computation");
    std::string scope = "chisquare2solve.";
    scope.append(scopeStr.c_str());

    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
    CPdfz pdfz;
    Float64 cstLog = -1;

    Int32 retPdfz=-1;
    std::vector<TFloat64List> priors;
    std::vector<TFloat64List> chiSquares;
    std::vector<Float64> redshifts;
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
        Int32 nISM = -1;
        Int32 nIGM = -1;
        if(meritResult->ChiSquareIntermediate.size()>0)
        {
            nISM = meritResult->ChiSquareIntermediate[0].size();
            if(meritResult->ChiSquareIntermediate[0].size()>0)
            {
                nIGM = meritResult->ChiSquareIntermediate[0][0].size();
            }
        }
        if(cstLog==-1)
        {
            cstLog = meritResult->CstLog;
        }else if ( cstLog != meritResult->CstLog)
        {
            Log.LogError("chisquare2solve: Found different cstLog values in results... val-1=%f != val-2=%f");
        }
        if(redshifts.size()==0)
        {
            redshifts = meritResult->Redshifts;
        }

        for(Int32 kism=0; kism<nISM; kism++)
        {
            for(Int32 kigm=0; kigm<nIGM; kigm++)
            {
                TFloat64List _prior;
                _prior = pdfz.GetConstantLogZPrior(meritResult->Redshifts.size());
                priors.push_back(_prior);

                //correct chi2 for ampl. marg. if necessary: todo add switch, currently deactivated
                TFloat64List logLikelihoodCorrected(meritResult->ChiSquareIntermediate.size(), DBL_MAX);
                for ( UInt32 kz=0; kz<meritResult->Redshifts.size(); kz++)
                {
                    logLikelihoodCorrected[kz] = meritResult->ChiSquareIntermediate[kz][kism][kigm];// + resultXXX->ScaleMargCorrectionTplshapes[][]?;
                }
                chiSquares.push_back(logLikelihoodCorrected);
            }
        }
    }

    if(opt_combine=="marg" && chiSquares.size()>0)
    {
        retPdfz = pdfz.Marginalize( redshifts, chiSquares, priors, cstLog, postmargZResult);
    }else if(opt_combine=="bestchi2")
    {
        retPdfz = pdfz.BestChi2( redshifts, chiSquares, priors, cstLog, postmargZResult);
    }else{
        Log.LogError("Chisquare2: Unable to parse pdf combination method option");
    }


    if(retPdfz!=0)
    {
        Log.LogError("chisquare2: Pdfz computation failed");
    }else{
        store.StoreGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
    }

    return retPdfz;
}




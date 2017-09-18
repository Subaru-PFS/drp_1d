#include <RedshiftLibrary/method/chisquarelogsolve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/operator/pdfLogresult.h>

#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <float.h>
using namespace NSEpic;
using namespace std;

CMethodChisquareLogSolve::CMethodChisquareLogSolve( std::string calibrationPath )
{
    bool enableLogRebin = true;
    m_chiSquareOperator = new COperatorChiSquareLogLambda( calibrationPath, enableLogRebin );
}

CMethodChisquareLogSolve::~CMethodChisquareLogSolve()
{
    delete m_chiSquareOperator;
}


const std::string CMethodChisquareLogSolve::GetDescription()
{
    std::string desc;

    desc = "Method chisquarelogsolve:\n";

    desc.append("\tparam: chisquarelogsolve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
    desc.append("\tparam: chisquarelogsolve.overlapThreshold = <float value>\n");
    desc.append("\tparam: chisquare.interpolation = {""precomputedfinegrid"", ""lin""}\n");
    desc.append("\tparam: chisquare.extinction = {""yes"", ""no""}\n");
    desc.append("\tparam: chisquare.dustfit = {""yes"", ""no""}\n");


    return desc;

}


std::shared_ptr<const CChisquareLogSolveResult> CMethodChisquareLogSolve::Compute(CDataStore& resultStore,
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

    CDataStore::CAutoScope resultScope( resultStore, "chisquarelogsolve" );
    std::string scopeStr = "chisquare";

    Int32 _type;
    if(spcComponent=="raw"){
       _type = CChisquareLogSolveResult::nType_raw;
       scopeStr = "chisquare";
    }else if(spcComponent=="nocontinuum"){
       _type = CChisquareLogSolveResult::nType_noContinuum;
       scopeStr = "chisquare_nocontinuum";
    }else if(spcComponent=="continuum"){
        _type = CChisquareLogSolveResult::nType_continuumOnly;
        scopeStr = "chisquare_continuum";
    }else if(spcComponent=="all"){
        _type = CChisquareLogSolveResult::nType_all;
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
        std::shared_ptr< CChisquareLogSolveResult>  ChisquareSolveResult = std::shared_ptr< CChisquareLogSolveResult>( new CChisquareLogSolveResult() );
        ChisquareSolveResult->m_type = _type;

        //std::string opt_combinePdf = "marg";
        std::string opt_combinePdf = "bestchi2";
        //std::string opt_combinePdf = "bestproba";
        CombinePDF(resultStore, scopeStr, opt_combinePdf);

        return ChisquareSolveResult;
    }

    return NULL;
}

Bool CMethodChisquareLogSolve::Solve(CDataStore& resultStore,
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
    Int32 _spctypetab[3] = {CChisquareLogSolveResult::nType_raw, CChisquareLogSolveResult::nType_noContinuum, CChisquareLogSolveResult::nType_continuumOnly};


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
    if(spctype == CChisquareLogSolveResult::nType_all){
        _ntype = 3;
    }

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == CChisquareLogSolveResult::nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = spctype;
        }

        if(_spctype == CChisquareLogSolveResult::nType_continuumOnly){
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
        }else if(_spctype == CChisquareLogSolveResult::nType_raw){
            // use full spectrum
            _spc = spc;
            _tpl = tpl;
            scopeStr = "chisquare";

        }else if(_spctype == CChisquareLogSolveResult::nType_noContinuum){
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
            bool enableSaveIntermediateChisquareResults = true;
            if(enableSaveIntermediateChisquareResults && chisquareResult->ChiSquareIntermediate.size()>0 && chisquareResult->ChiSquareIntermediate.size()==chisquareResult->Redshifts.size())
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


Int32 CMethodChisquareLogSolve::CombinePDF(CDataStore &store, std::string scopeStr, std::string opt_combine )
{
    Log.LogInfo("chisquarelogsolve: Pdfz computation");
    std::string scope = "chisquarelogsolve.";
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
            Log.LogError("chisquarelogsolve: Found different cstLog values in results... val-1=%f != val-2=%f");
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
        Log.LogError("Chisquarelog: Unable to parse pdf combination method option");
    }


    if(retPdfz!=0)
    {
        Log.LogError("Chisquarelog: Pdfz computation failed");
    }else{
        store.StoreGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
    }

    return retPdfz;
}


Int32 CMethodChisquareLogSolve::CombinePDF_deprecated(CDataStore &store, std::string scopeStr )
{
    std::string scope = "chisquarelogsolve.";
    scope.append(scopeStr.c_str());

    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    auto postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
    Bool initPostMarg = false;
    std::vector<UInt32> nSum;

    Float64 MaxiLogEvidence = -DBL_MAX;
    TFloat64List LogEvidences;
    Float64 sumModifiedEvidences = 0;
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );

        //Todo: Check if the status is OK ?
        //meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;

        std::vector<Float64> zPrior=pdfz.GetConstantLogZPrior(meritResult->Redshifts.size());
        Int32 retPdfz = pdfz.Compute(meritResult->ChiSquare, meritResult->Redshifts, meritResult->CstLog, zPrior, logProba, logEvidence);
        if(retPdfz!=0)
        {
            Log.LogError("chisquarelogsolve: Pdfz computation failed for tpl %s", (*it).first.c_str());
        }else{
            LogEvidences.push_back(logEvidence);
            if(MaxiLogEvidence<logEvidence){
                MaxiLogEvidence=logEvidence;
            }
        }
    }
    //Using computational trick to sum the evidences
    for(Int32 k=0; k<LogEvidences.size(); k++)
    {
        sumModifiedEvidences += exp(LogEvidences[k]-MaxiLogEvidence);
    }
    Float64 logSumEvidence = MaxiLogEvidence + log(sumModifiedEvidences);

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );

        //Todo: Check if the status is OK ?
        //meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        std::vector<Float64> zPrior=pdfz.GetConstantLogZPrior(meritResult->Redshifts.size());
        Int32 retPdfz = pdfz.Compute(meritResult->ChiSquare, meritResult->Redshifts, meritResult->CstLog, zPrior, logProba, logEvidence);
        if(retPdfz!=0)
        {
            Log.LogError("chisquarelogsolve: Pdfz computation failed for tpl %s", (*it).first.c_str());
        }else{
            if(!initPostMarg)
            {
                nSum.resize(meritResult->Redshifts.size());
                postmargZResult->countTPL = meritResult->Redshifts.size(); // assumed 1 model per z
                postmargZResult->Redshifts.resize(meritResult->Redshifts.size());
                postmargZResult->valProbaLog.resize(meritResult->Redshifts.size());
                for ( UInt32 k=0; k<meritResult->Redshifts.size(); k++)
                {
                    postmargZResult->Redshifts[k] = meritResult->Redshifts[k] ;
                    postmargZResult->valProbaLog[k] = log(0.0);
                    nSum[k] = 0;
                }
                initPostMarg = true;
            }else
            {
                //check if the redshift bins are the same
                for ( UInt32 k=0; k<meritResult->Redshifts.size(); k++)
                {
                    if(postmargZResult->Redshifts[k] != meritResult->Redshifts[k])
                    {
                        Log.LogError("chisquarelogsolve: Pdfz computation (z-bins comparison) failed for tpl %s", (*it).first.c_str());
                        break;
                    }
                }
            }
            for ( UInt32 k=0; k<meritResult->Redshifts.size(); k++)
            {
                if(meritResult->Status[k]== COperator::nStatus_OK)
                {
                    Float64 valProba = exp(postmargZResult->valProbaLog[k]);
                    Float64 valProbaAdd = exp(logProba[k]+logEvidence-logSumEvidence);
                    postmargZResult->valProbaLog[k] = log( valProba + valProbaAdd );
                    nSum[k]++;
                }
            }
        }


    }

    //THIS DOES NOT ALLOW Marginalization with coverage<100% for ALL templates
    for ( UInt32 k=0; k<postmargZResult->Redshifts.size(); k++)
    {
        if(nSum[k]!=meritResults.size())
        {
            postmargZResult->valProbaLog[k] = NAN;
            Log.LogError("chisquarelogsolve: Pdfz computation failed. Not all templates have 100% coverage for all redshifts!");
        }
    }
    store.StoreGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF

    return 0;
}




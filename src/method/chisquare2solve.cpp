#include <epic/redshift/method/chisquare2solve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/operator/chisquare2.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>


#include <epic/redshift/spectrum/io/fitswriter.h>

using namespace NSEpic;
using namespace std;

CMethodChisquare2Solve::CMethodChisquare2Solve()
{

}

CMethodChisquare2Solve::~CMethodChisquare2Solve()
{

}


const std::string CMethodChisquare2Solve::GetDescription()
{
    std::string desc;

    desc = "Method chisquare2solve:\n";

    desc.append("\tparam: chisquare2solve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
    desc.append("\tparam: chisquare2solve.overlapThreshold = <float value>\n");
    desc.append("\tparam: chisquare.interpolation = {""precomputedfinegrid"", ""lin""}\n");
    desc.append("\tparam: chisquare.extinction = {""yes"", ""no""}\n");


    return desc;

}


std::shared_ptr<const CChisquare2SolveResult> CMethodChisquare2Solve::Compute(  CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, std::string spcComponent, std::string opt_interp, std::string opt_extinction)
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "chisquare2solve" );

    Int32 _type;
    if(spcComponent=="raw"){
       _type = CChisquare2SolveResult::nType_raw;
    }else if(spcComponent=="nocontinuum"){
       _type = CChisquare2SolveResult::nType_noContinuum;
    }else if(spcComponent=="continuum"){
        _type = CChisquare2SolveResult::nType_continuumOnly;
    }else if(spcComponent=="all"){
        _type = CChisquare2SolveResult::nType_all;
    }


    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );

            /*
            //optionnally export the templates into fits files
            if(1){
                CSpectrum spctpl= CSpectrum(tpl);
                std::string tplName = tpl.GetName();
                tplName.append(".fits");

                CSpectrumIOFitsWriter writer;
                Bool retValueExportTpl = writer.Write( tplName.c_str(),  spctpl);
            }
            //*/

            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            Solve( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshifts, overlapThreshold, _type, opt_interp, opt_extinction);

            storeResult = true;
        }
    }


    if( storeResult )
    {
        std::shared_ptr< CChisquare2SolveResult>  ChisquareSolveResult = std::shared_ptr< CChisquare2SolveResult>( new CChisquare2SolveResult() );
        ChisquareSolveResult->m_type = _type;
        return ChisquareSolveResult;
    }

    return NULL;
}

Bool CMethodChisquare2Solve::Solve( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, Int32 spctype, std::string opt_interp, std::string opt_extinction )
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
        }

        // Compute merit function
        COperatorChiSquare2 chiSquare;
        //CRef<CChisquareResult>  chisquareResult = (CChisquareResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
        auto  chisquareResult = std::dynamic_pointer_cast<CChisquareResult>( chiSquare.Compute( _spc, _tpl, lambdaRange, redshifts, overlapThreshold, opt_interp, enable_extinction ) );

        if( !chisquareResult )
        {
            //Log.LogInfo( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            resultStore.StoreScopedPerTemplateResult( tpl, scopeStr.c_str(), chisquareResult );
        }
    }

    return true;
}

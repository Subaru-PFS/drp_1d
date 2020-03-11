#include <RedshiftLibrary/linemodel/templatesortho.h>
#include <RedshiftLibrary/linemodel/elementlist.h>



using namespace NSEpic;


CTemplatesOrthogonalization::CTemplatesOrthogonalization(const CTemplateCatalog& tplCatalog,
                                                         const TStringList& tplCategoryList,
                                                         const std::string calibrationPath,
                                                         const CRayCatalog::TRayVector& restRayList,
                                                         const std::string& opt_fittingmethod,
                                                         const std::string& opt_continuumcomponent,
                                                         const std::string& widthType,
                                                         const Float64 opt_nsigmasupport,
                                                         const Float64 resolution,
                                                         const Float64 velocityEmission,
                                                         const Float64 velocityAbsorption,
                                                         const std::string& opt_rules,
                                                         const std::string& opt_rigidity,
                                                         bool enableOrtho)
{

    m_enableOrtho = enableOrtho;

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );


            std::string rigidity = opt_rigidity.c_str();
            std::string rules = opt_rules.c_str();
            //temporary options override to be removed when full tpl ortho is implemented
            Bool enableOverride = true;
            if(enableOverride){
                rigidity = "rules";
                rules = "no";
            }
            std::string opt_fittingmethod2 = "hybrid";

            Log.LogDetail("    tplOrthogonalization: now processing tpl=%s", tpl.GetName().c_str() );
            Int32 ret = OrthogonalizeTemplate(tpl,
                                  calibrationPath,
                                  restRayList,
                                  opt_fittingmethod2,
                                  widthType,
                                  opt_nsigmasupport,
                                  resolution,
                                  velocityEmission,
                                  velocityAbsorption,
                                  rules,
                                  rigidity);
           if(ret!=0)
           {
               //do something...
           }
        }
    }

    std::shared_ptr<CTemplateCatalog> tplCatalogPtr = std::shared_ptr<CTemplateCatalog>( new CTemplateCatalog(m_tplCatalogOrthogonal) );
//    std::shared_ptr<CTemplate> tplOrtho = std::shared_ptr<CTemplate>( new CTemplate( inputTemplate.GetName().c_str(), inputTemplate.GetCategory() ) );

    m_tplOrthoStore.Add(tplCatalogPtr);

}

/**
 * \brief Empty destructor.
 **/
CTemplatesOrthogonalization::~CTemplatesOrthogonalization()
{

}

/**
 * @brief getOrthogonalTemplate
 * Computes as follows:
 * - process linemodel on the spectrum
 * - creates the output template from the subtraction of the input spectrum and the fitted linemodel
 * - add the newly created template to the tplCatalogOrtho member
 * @return
 */
Int32 CTemplatesOrthogonalization::OrthogonalizeTemplate(const CTemplate& inputTemplate,
                            const std::string opt_calibrationPath,
                            const CRayCatalog::TRayVector &restRayList,
                            const std::string &opt_fittingmethod,
                            const std::string &opt_lineWidthType,
                            const Float64 opt_nsigmasupport,
                            const Float64 opt_resolution,
                            const Float64 opt_velocityEmission,
                            const Float64 opt_velocityAbsorption,
                            const std::string &opt_rules,
                            const std::string &opt_rigidity)
{

    std::shared_ptr<CTemplate> tplOrtho = std::shared_ptr<CTemplate>( new CTemplate( inputTemplate.GetName().c_str(), inputTemplate.GetCategory() ) );
    *tplOrtho = inputTemplate; //todo: check if this is a true copy of the samples

    bool enableModelSubtraction = m_enableOrtho;
    if(enableModelSubtraction){

        std::string opt_continuumcomponent = "fromspectrum";
        CSpectrum spectrum = CSpectrum(inputTemplate);
        CSpectrum spectrumContinuumZero = CSpectrum(inputTemplate);
        CSpectrumFluxAxis& continuumZeroFluxAxis = spectrumContinuumZero.GetFluxAxis();
        for(UInt32 i=0; i<continuumZeroFluxAxis.GetSamplesCount(); i++){
            continuumZeroFluxAxis[i] = 0.0; //put zero as continuum here
        }

        CTemplateCatalog tplCatalogUnused;
        TStringList tplCategoryListUnused;


        //Compute linemodel on the template
        CLineModelElementList model( spectrum,
                                     spectrumContinuumZero,
                                     tplCatalogUnused,
                                     tplCatalogUnused,
                                     tplCategoryListUnused,
                                     opt_calibrationPath,
                                     restRayList,
                                     opt_fittingmethod,
                                     opt_continuumcomponent,
                                     opt_lineWidthType,
                                     opt_nsigmasupport,
                                     opt_resolution,
                                     opt_velocityEmission,
                                     opt_velocityAbsorption,
                                     opt_rules,
                                     opt_rigidity);

        Float64 redshift = 0.0;
        TLambdaRange lambdaRange = inputTemplate.GetLambdaRange();
        Float64 contreest_iterations = 0;
        Bool enableLogging=true;
        CLineModelSolution modelSolution;
        CContinuumModelSolution continuumModelSolution;
        model.fit( redshift,
                   lambdaRange,
                   modelSolution,
                   continuumModelSolution,
                   contreest_iterations,
                   enableLogging );

        //get mtm and dtm cumulative vector and store it
        std::vector<Float64> lbda;
        std::vector<Float64> mtmCumul;
        model.getMTransposeMCumulative(lambdaRange, lbda, mtmCumul);


        //Subtract the fitted model from the original template
        model.refreshModel();
        CSpectrum modelSpc = model.GetModelSpectrum();
        /*//debug:
        FILE* f = fopen( "templatesortho_fittedmodel_dbg.txt", "w+" );
        for( Int32 t=0;t<modelSpc.GetSampleCount();t++)
        {
            fprintf( f, "%f %e\n", modelSpc.GetSpectralAxis()[t], modelSpc.GetFluxAxis()[t]);
        }
        fclose( f );
        //*/

        CSpectrumFluxAxis& modelFluxAxis = modelSpc.GetFluxAxis();
        CSpectrumFluxAxis& continuumOrthoFluxAxis = tplOrtho->GetFluxAxis();
        for(UInt32 i=0; i<continuumOrthoFluxAxis.GetSamplesCount(); i++){
            continuumOrthoFluxAxis[i] -= modelFluxAxis[i];
        }
        /*//debug:
        FILE* f2 = fopen( "templatesortho_orthotemplate_dbg.txt", "w+" );
        for( Int32 t=0;t<modelSpc.GetSampleCount();t++)
        {
            fprintf( f2, "%f %e\n", modelSpc.GetSpectralAxis()[t], continuumOrthoFluxAxis[t]);
        }
        fclose( f2 );
        //*/

    }


    m_tplCatalogOrthogonal.Add(tplOrtho);
    return 0;
}

CTemplateCatalog CTemplatesOrthogonalization::getOrthogonalTplCatalog()
{
    return m_tplCatalogOrthogonal;
}

CTemplatesOrthoStore CTemplatesOrthogonalization::getOrthogonalTplStore()
{
    return m_tplOrthoStore;
}


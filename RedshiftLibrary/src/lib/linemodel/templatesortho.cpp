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
#include "RedshiftLibrary/linemodel/templatesortho.h"
#include "RedshiftLibrary/linemodel/elementlist.h"



using namespace NSEpic;


CTemplatesOrthogonalization::CTemplatesOrthogonalization(const CTemplateCatalog& tplCatalog,
                                                         const TStringList& tplCategoryList,
                                                         const std::string calibrationPath,
                                                         const CRayCatalog::TRayVector& restRayList,
                                                         const std::string& opt_fittingmethod,
                                                         const std::string& widthType,
                                                         const Float64 opt_nsigmasupport,
                                                         const Float64 velocityEmission,
                                                         const Float64 velocityAbsorption,
                                                         const std::string& opt_rules,
                                                         const std::string& opt_rigidity,
                                                         std::shared_ptr<const CLSF> lsf,
                                                         bool enableOrtho):
m_LSF(lsf)                                                        
{

    m_enableOrtho = enableOrtho;
    Bool currentsampling = tplCatalog.m_logsampling; 
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];
        for(Bool sampling:{0, 1})
        {
            tplCatalog.m_logsampling = sampling;
            m_tplCatalogOrthogonal.m_logsampling = sampling; 
            TTemplateConstRefList  TplList = tplCatalog.GetTemplateList(TStringList{category});
            for(auto tpl:TplList )
            {
                std::string rigidity = opt_rigidity.c_str();
                std::string rules = opt_rules.c_str();
                //temporary options override to be removed when full tpl ortho is implemented
                Bool enableOverride = true;
                if(enableOverride){
                    rigidity = "rules";
                    rules = "no";
                }

                Log.LogDetail("    tplOrthogonalization: now processing tpl=%s", tpl->GetName().c_str() );
                Int32 ret = OrthogonalizeTemplate(*tpl,
                                    calibrationPath,
                                    restRayList,
                                    opt_fittingmethod,
                                    widthType,
                                    opt_nsigmasupport,
                                    velocityEmission,
                                    velocityAbsorption,
                                    rules,
                                    rigidity);
            if(ret!=0)
            {
              Log.LogError("error in orthogonalization");
                //do something...
            }
            }
        }
    }
    tplCatalog.m_logsampling = currentsampling;
    m_tplCatalogOrthogonal.m_logsampling = currentsampling; 
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
                            const Float64 opt_velocityEmission,
                            const Float64 opt_velocityAbsorption,
                            const std::string &opt_rules,
                            const std::string &opt_rigidity)
{

    std::shared_ptr<CTemplate> tplOrtho = std::make_shared<CTemplate>(inputTemplate);

    bool enableModelSubtraction = m_enableOrtho;
    if(enableModelSubtraction){

        std::string opt_continuumcomponent = "fromspectrum";
        Float64 opt_continuum_neg_threshold=-INFINITY; // not relevant in the "fromspectrum" case
        CSpectrum spectrum = inputTemplate;
        spectrum.SetLSF(m_LSF);

        std::string saveContinuumEstimationMethod = spectrum.GetContinuumEstimationMethod();
        spectrum.SetContinuumEstimationMethod("zero");

        CTemplateCatalog tplCatalogUnused;
        TStringList tplCategoryListUnused;

        //Compute linemodel on the template
        CLineModelElementList model( spectrum,
                                     tplCatalogUnused,
                                     tplCategoryListUnused,
                                     opt_calibrationPath,
                                     restRayList,
                                     opt_fittingmethod,
                                     opt_continuumcomponent,
                                     opt_continuum_neg_threshold,
                                     opt_lineWidthType,
                                     opt_nsigmasupport,
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

        //Restore the continuum estimation method
        spectrum.SetContinuumEstimationMethod(saveContinuumEstimationMethod);

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

        const CSpectrumFluxAxis& modelFluxAxis = modelSpc.GetFluxAxis();
        CSpectrumFluxAxis continuumOrthoFluxAxis = tplOrtho->GetFluxAxis();
        for(UInt32 i=0; i<continuumOrthoFluxAxis.GetSamplesCount(); i++){
            continuumOrthoFluxAxis[i] -= modelFluxAxis[i];
        }
        tplOrtho->SetFluxAxis(std::move(continuumOrthoFluxAxis));
        
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


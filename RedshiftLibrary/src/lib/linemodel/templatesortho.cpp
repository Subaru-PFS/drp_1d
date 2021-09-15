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
#include "RedshiftLibrary/spectrum/logrebinning.h"



using namespace NSEpic;
void CTemplatesOrthogonalization::Orthogonalize(CInputContext& inputContext, 
                                                const std::string category, 
                                                const std::string calibrationPath,
                                                std::shared_ptr<const CLSF> lsf)
{
    //retrieve params from InputContext
    m_LSF = lsf;
    m_enableOrtho = inputContext.GetParameterStore()->EnableTemplateOrthogonalization(category);
    
    Float64 opt_nsigmasupport;
    Float64 velocityEmission, velocityAbsorption;

    std::string widthType;
    std::string opt_rules, opt_rigidity;
    std::string opt_lineforcefilter, opt_linetypefilter;

    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.linewidthtype", widthType, "velocitydriven" );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.nsigmasupport", opt_nsigmasupport, 8.0 );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.velocityemission", velocityEmission, 200.0 );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.velocityabsorption", velocityAbsorption, 300.0 );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.rules", opt_rules, "all" );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.rigidity", opt_rigidity, "rules" );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.linetypefilter", opt_linetypefilter, "no" );
    inputContext.GetParameterStore()->Get( category + "linemodelsolve.linemodel.lineforcefilter", opt_lineforcefilter, "no" );

    std::string rigidity = opt_rigidity.c_str();
    std::string rules = opt_rules.c_str();
    //temporary options override to be removed when full tpl ortho is implemented
    Bool enableOverride = true;
    if(enableOverride){
        rigidity = "rules";
        rules = "no";
    }

    Int32 typeFilter = -1;
    if (opt_linetypefilter == "A") typeFilter = CRay::nType_Absorption;
    else if (opt_linetypefilter == "E")typeFilter = CRay::nType_Emission;
    Int32 forceFilter = -1; // CRay::nForce_Strong;
    if (opt_lineforcefilter == "S")
    forceFilter = CRay::nForce_Strong;

    CRayCatalog::TRayVector restRayList = inputContext.GetRayCatalog(category)->GetFilteredList(typeFilter, forceFilter);
    // prepare continuum templates catalog
    std::string opt_fittingmethod="hybrid";

    //retrieve templateCatalog
    std::shared_ptr<CTemplateCatalog> tplCatalog = inputContext.GetTemplateCatalog();
    Bool currentsampling = tplCatalog->m_logsampling; 

    tplCatalog->m_logsampling = 1; tplCatalog->m_orthogonal = 0;//orig log
    Bool rebinning = tplCatalog->GetTemplateCount(category) > 0;//is rebinning is true, ortho-log tpls are cleared obligatory

    //check if LSF has changed, if yes reorthog all
    bool differentLSF = false;
    std::vector<Bool> samplingList {0,1};
    Float64 lambda = (inputContext.m_lambdaRange.GetBegin() + inputContext.m_lambdaRange.GetEnd())/2;
    if(std::isnan(tplCatalog->m_ortho_LSFWidth)){
        tplCatalog->m_ortho_LSFWidth = m_LSF->GetWidth(lambda); //first time orthogonalizing - do it for all
        differentLSF = true;
    }else{
        if(tplCatalog->m_ortho_LSFWidth != m_LSF->GetWidth(lambda)) // ortho setting changed - do it for all
        {   
            differentLSF = true;
            tplCatalog->ClearTemplateList(category, 1, 0);//clear orthog templates - non-rebinned
            tplCatalog->ClearTemplateList(category, 1, 1);//clear orthog templates - rebinned
            tplCatalog->m_ortho_LSFWidth = m_LSF->GetWidth(lambda);
        }else{
            if(rebinning) //ortho setting is the same but rebinning happened, only log-ortho templates are concerned
                samplingList = {1};
        }
    }
    Bool needOrthogonalization = rebinning | differentLSF;
    if(!needOrthogonalization) return;

    for(Bool sampling:samplingList)
    {
        tplCatalog->m_logsampling = sampling;
        //orthogonalize all templates
        tplCatalog->m_orthogonal = 0;
        const TTemplateConstRefList TplList = std::const_pointer_cast<const CTemplateCatalog>(tplCatalog)->GetTemplateList(TStringList{category});
        tplCatalog->m_orthogonal = 1;
        for (auto tpl : TplList)
        {   Log.LogDetail("    TplOrthogonalization: now processing tpl=%s", tpl->GetName().c_str() );
            std::shared_ptr<CTemplate> _orthoTpl = OrthogonalizeTemplate(*tpl,
                                                                        calibrationPath,
                                                                        restRayList,
                                                                        opt_fittingmethod,
                                                                        widthType,
                                                                        opt_nsigmasupport,
                                                                        velocityEmission,
                                                                        velocityAbsorption,
                                                                        rules,
                                                                        rigidity);
            tplCatalog->Add(_orthoTpl);
        }
    }    
    tplCatalog->m_logsampling = currentsampling;
    tplCatalog->m_orthogonal = 0; 
    return;
}

/**
 * @brief getOrthogonalTemplate
 * Computes as follows:
 * - process linemodel on the spectrum
 * - creates the output template from the subtraction of the input spectrum and the fitted linemodel
 * - add the newly created template to the tplCatalogOrtho member
 * @return
 */
std::shared_ptr<CTemplate> CTemplatesOrthogonalization::OrthogonalizeTemplate(const CTemplate& inputTemplate,
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
    return tplOrtho;
}


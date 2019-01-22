#ifndef _REDSHIFT_OPERATOR_ZWEIMODELSOLVE_
#define _REDSHIFT_OPERATOR_ZWEIMODELSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/linemodel.h>
#include <RedshiftLibrary/linemodel/zweimodelresult.h>


namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CZweiModelSolve
{

public:

    CZweiModelSolve(std::string calibrationPath="");
    ~CZweiModelSolve();

    const std::string GetDescription();
    Bool PopulateParameters( CDataStore& dataStore );


    std::shared_ptr<CLineModelSolveResult> Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog& restraycatalog,
                                           const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog& restraycatalog,
                                 const TFloat64Range& lambdaRange, const TFloat64List& redshifts);

private:

    Int32 getValueFromRefFile( const char* filePath, std::string spcid, Int32 colID, Float64& zref, Int32 reverseInclusion );
    Int32 getVelocitiesFromRefFile( const char* filePath, std::string spcid, Float64& elv, Float64& alv );
    Int32 getContaminantNameFromFile(const char* filePath, std::string spcid, Int32 colID,
                                     std::string& contSpcFileName,
                                     std::string& contErrorFileName,
                                     Float64 &offsetLambdaContaminant,
                                     Int32 reverseInclusion );


    Int32 CombinePDF(CDataStore &store,
                     std::shared_ptr<const CLineModelResult> result,
                     std::string opt_rigidity,
                     std::string opt_combine,
                     Float64 opt_stronglinesprior);
    Int32 SaveContinuumPDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result);


    std::string m_opt_linetypefilter;
    std::string m_opt_lineforcefilter;
    std::string m_opt_fittingmethod;
    std::string m_opt_continuumcomponent;
    std::string m_opt_rigidity;
    std::string m_opt_lineWidthType;
    Float64 m_opt_resolution;
    Float64 m_opt_velocity_emission;
    Float64 m_opt_velocity_absorption;
    std::string m_opt_velocityfit;
    Float64 m_opt_em_velocity_fit_min;
    Float64 m_opt_em_velocity_fit_max;
    Float64 m_opt_abs_velocity_fit_min;
    Float64 m_opt_abs_velocity_fit_max;
    std::string m_opt_continuumreest;
    std::string m_opt_rules;

    //options for rigidity=tplshape
    std::string m_opt_tplratio_reldirpath;
    std::string m_opt_tplratio_ismfit="no";
    std::string m_opt_offsets_reldirpath;

    Float64 m_opt_extremacount;
    Float64 m_opt_firstpass_largegridstep;
    std::string m_opt_firstpass_largegridsampling;
    std::string m_opt_pdfcombination;
    Float64 m_opt_stronglinesprior;
    Float64 m_opt_euclidNHaEmittersPriorStrength;
    std::string m_opt_saveintermediateresults;
    bool m_opt_enableSaveChisquareTplshapeResults;

    std::string m_calibrationPath;
};


}

#endif

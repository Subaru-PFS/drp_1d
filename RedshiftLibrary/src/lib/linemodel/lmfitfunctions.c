
#include <RedshiftLibrary/linemodel/elementlist.h>

namespace NSEpic
{
class CLineModelElementList;


struct lmfitdata {
    size_t n;
    Float64 * y;
    CLineModelElementList* linemodel;
    std::vector<Int32> linemodel_elts_indexes;
    std::vector<Int32> linemodel_samples_indexes;
    Float64 normFactor;
    Int32 lineType;
    CTemplate tpl;
};

int
lmfit_f (const gsl_vector * x, void *data,
         gsl_vector * f)
{
    size_t n = ((struct lmfitdata *)data)->n;
    Float64 *y = ((struct lmfitdata *)data)->y;
    //std::shared_ptr<CLineModelElementList> linemodel = ((struct lmfitdata *)data)->linemodel;
    CLineModelElementList* linemodel = ((struct lmfitdata *)data)->linemodel;
    std::vector<Int32> elts_indexes = ((struct lmfitdata *)data)->linemodel_elts_indexes;
    std::vector<Int32> samples_indexes = ((struct lmfitdata *)data)->linemodel_samples_indexes;
    Float64 normFactor = ((struct lmfitdata *)data)->normFactor;
    Int32 lineType = ((struct lmfitdata *)data)->lineType;
    CTemplate tpl = ((struct lmfitdata *)data)->tpl;

    Int32 idxTplAmp = elts_indexes.size()+1;
    Float64 tplAmp= gsl_vector_get (x, idxTplAmp);
    linemodel.setFitContinuum_tplAmplitude( tplAmp)

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
    {
        Float64 amp = gsl_vector_get (x, iElt)/normFactor;
        linemodel->SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
    }

    Int32 idxVelocity = elts_indexes.size();
    Float64 velocity = gsl_vector_get (x, idxVelocity);
    if(lineType==CRay::nType_Emission)
    {
        linemodel->SetVelocityEmission(velocity);
    }else
    {
        linemodel->SetVelocityAbsorption(velocity);
    }
    linemodel->refreshModelUnderElements(elts_indexes);


    const CSpectrumFluxAxis& tplFluxAxis = tpl->GetFluxAxis();
    for (Int32 i = 0; i < n; i++)
    {

        Float64 Yi = linemodel->getModelFluxVal(samples_indexes[i])*normFactor;
        Float64 err = Yi - y[i]
        // TODO : add condition for fullmodel
        if(1){
          err += tplAmp * tplFluxAxis[samples_indexes[i]]* normFactor
        }
        gsl_vector_set (f, i, err);
    }

    return GSL_SUCCESS;
}

int
lmfit_df (const gsl_vector * x, void *data,
          gsl_matrix * J)
{
    size_t n = ((struct lmfitdata *)data)->n;
    //std::shared_ptr<CLineModelElementList> linemodel = ((struct lmfitdata *)data)->linemodel;
    CLineModelElementList* linemodel = ((struct lmfitdata *)data)->linemodel;
    std::vector<Int32> elts_indexes = ((struct lmfitdata *)data)->linemodel_elts_indexes;
    std::vector<Int32> samples_indexes = ((struct lmfitdata *)data)->linemodel_samples_indexes;
    Float64 normFactor = ((struct lmfitdata *)data)->normFactor;
    Int32 lineType = ((struct lmfitdata *)data)->lineType;

    Int32 idxTplAmp = elts_indexes.size()+1;
    Float64 tplAmp= gsl_vector_get (x, idxTplAmp);
    linemodel.setFitContinuum_tplAmplitude( tplAmp)

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
    {
        Float64 amp = gsl_vector_get (x, iElt)/normFactor;
        linemodel->SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
    }

    Int32 idxVelocity = elts_indexes.size();
    Float64 velocity = gsl_vector_get (x, idxVelocity);
    if(lineType==CRay::nType_Emission)
    {
        linemodel->SetVelocityEmission(velocity);
    }else
    {
        linemodel->SetVelocityAbsorption(velocity);
    }
    linemodel->refreshModelDerivSigmaUnderElements(elts_indexes);
    const CSpectrumFluxAxis& tplFluxAxis = tpl->GetFluxAxis();
    for (Int32 i = 0; i < n; i++)
    {
        for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
        {
            //TODO : demander confirmation pour l'ajout du normfactor
            Float64 dval = linemodel->getModelFluxDerivEltVal(elts_indexes[iElt], samples_indexes[i])*normFactor;
            gsl_matrix_set (J, i, iElt, dval);
        }
        Float64 dval = linemodel->getModelFluxDerivSigmaVal(samples_indexes[i])*normFactor;
        gsl_matrix_set (J, i, idxVelocity, dval);

        // TODO : add condition for fullmodel
        if(1){

          gsl_matrix_set (J, i, idxTplAmp, tplFluxAxis[samples_indexes[i]]);
        }
    }

    return GSL_SUCCESS;
}

}

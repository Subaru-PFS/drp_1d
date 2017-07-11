
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/lmfitcontroller.h>

#include <RedshiftLibrary/log/log.h>
namespace NSEpic
{
class CLineModelElementList;


struct lmfitdata {
    size_t n;
    Float64 * y;
    CLineModelElementList* linemodel;
    std::vector<Int32> linemodel_samples_indexes;
    const Float64* observeGridContinuumFlux;
    CLmfitController* controller;
};

int
lmfit_f (const gsl_vector * x, void *data,
         gsl_vector * f)
{
    size_t n = ((struct lmfitdata *)data)->n;
    Float64 *y = ((struct lmfitdata *)data)->y;
    //std::shared_ptr<CLineModelElementList> linemodel = ((struct lmfitdata *)data)->linemodel;
    CLineModelElementList* linemodel = ((struct lmfitdata *)data)->linemodel;
    CLmfitController* controller = ((struct lmfitdata *)data)->controller;
    std::vector<Int32> elts_indexes = controller->getFilteredIdx();
    std::vector<Int32> samples_indexes = ((struct lmfitdata *)data)->linemodel_samples_indexes;
    Float64 normFactor = controller->getNormFactor();
    //const Float64* observeGridContinuumFlux = ((struct lmfitdata *)data)->observeGridContinuumFlux;

    if(controller->isRedshiftFitted()){
      Float64 redshift = gsl_vector_get(x, controller-> getIndRedshift());
      linemodel->setRedshift(redshift , false);
      if(!controller->isContinuumFitted()){
        linemodel->setFitContinuum_tplAmplitude(linemodel->getFitContinuum_tplAmplitude());
      }
    }

    if(controller->isContinuumFitted()){
      Int32 idxContinuumTplAmp = controller->getIndContinuumAmp();
      Float64 continuumTplAmp= controller->continuumAmp_LmToModel(gsl_vector_get (x, idxContinuumTplAmp));///normFactor;
      linemodel->setFitContinuum_tplAmplitude( continuumTplAmp);
    }

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
    {
        Float64 amp = controller->lineAmp_LmToModel(gsl_vector_get (x, iElt));///normFactor;
        linemodel->SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
    }

    if(controller->isEmissionVelocityFitted()){
      Float64 velocity = controller->emiVel_LmToModel(gsl_vector_get (x, controller->getIndEmissionVel()));
      // Log.LogInfo("velocity emission val : %f", velocity);
      linemodel->SetVelocityEmission(velocity);
    }
    if(controller->isAbsorptionVelocityFitted()){
      Float64 velocity = controller->absVel_LmToModel(gsl_vector_get (x, controller->getIndAbsorptionVel()));
      linemodel->SetVelocityAbsorption(velocity);
    }

    //
    if(controller->isContinuumFitted()){
      linemodel->refreshModelInitAllGrid();
    }else{
      linemodel->refreshModelUnderElements(elts_indexes);
    }


    for (Int32 i = 0; i < n; i++)
    {

        Float64 Yi = linemodel->getModelFluxVal(samples_indexes[i])*normFactor;


        gsl_vector_set (f, i, Yi - y[i]);
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
    CLmfitController* controller = ((struct lmfitdata *)data)->controller;
    std::vector<Int32> elts_indexes = controller->getFilteredIdx();
    std::vector<Int32> samples_indexes = ((struct lmfitdata *)data)->linemodel_samples_indexes;
    Float64 normFactor = controller->getNormFactor();
    const Float64* observeGridContinuumFlux = ((struct lmfitdata *)data)->observeGridContinuumFlux;


    if(controller->isRedshiftFitted()){
      Float64 redshift = gsl_vector_get(x, controller-> getIndRedshift());
      //Log.LogInfo("redshift value %f", redshift);
      linemodel->setRedshift(redshift ,false );
      if(!controller->isContinuumFitted()){
        linemodel->setFitContinuum_tplAmplitude(linemodel->getFitContinuum_tplAmplitude());
      }
    }

    if(controller->isContinuumFitted()){
      Int32 idxContinuumTplAmp = controller->getIndContinuumAmp();
      Float64 continuumTplAmp= controller->continuumAmp_LmToModel(gsl_vector_get (x, idxContinuumTplAmp));///normFactor;
      linemodel->setFitContinuum_tplAmplitude( continuumTplAmp);
    }

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
    {
        Float64 amp = controller->lineAmp_LmToModel(gsl_vector_get (x, iElt));///normFactor;
        linemodel->SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
    }

    if(controller->isEmissionVelocityFitted()){
      Float64 velocity = controller->emiVel_LmToModel(gsl_vector_get (x, controller->getIndEmissionVel()));
      // Log.LogInfo("velocity emission val : %f", velocity);
      linemodel->SetVelocityEmission(velocity);
    }
    if(controller->isAbsorptionVelocityFitted()){
      Float64 velocity = controller->absVel_LmToModel(gsl_vector_get (x, controller->getIndAbsorptionVel()));
      linemodel->SetVelocityAbsorption(velocity);
    }

    if(controller->isAbsorptionVelocityFitted() && controller->isEmissionVelocityFitted()){
      linemodel->refreshModelDerivVelUnderElements(elts_indexes);
    }else if(controller->isAbsorptionVelocityFitted()){
      linemodel->refreshModelDerivVelAbsorptionUnderElements(elts_indexes);
    }else if(controller->isEmissionVelocityFitted()){
      linemodel->refreshModelDerivVelEmissionUnderElements(elts_indexes);
    }
    //linemodel->refreshModel();
    for (Int32 i = 0; i < n; i++)
    {
        for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
        {

            Float64 dval = linemodel->getModelFluxDerivEltVal(elts_indexes[iElt], samples_indexes[i])*normFactor*2*gsl_vector_get (x, iElt);
            gsl_matrix_set (J, i, iElt, dval);
        }

        if(controller->isEmissionVelocityFitted()){
            Float64 normEmiFactor = controller->getNormEmiFactor();
          Float64 dval = linemodel->getModelFluxDerivVelEmissionVal(samples_indexes[i])*normFactor/normEmiFactor*2*gsl_vector_get (x, controller->getIndEmissionVel());
          gsl_matrix_set (J, i, controller->getIndEmissionVel(), dval);
        }

        if(controller->isAbsorptionVelocityFitted()){
          Float64 normAbsFactor = controller->getNormAbsFactor();
          Float64 dval = linemodel->getModelFluxDerivVelAbsorptionVal(samples_indexes[i])*normFactor/normAbsFactor*2*gsl_vector_get (x, controller->getIndAbsorptionVel());
          gsl_matrix_set (J, i, controller->getIndAbsorptionVel(), dval);
        }

        if(controller->isContinuumFitted()){
          Float64 dval =  observeGridContinuumFlux[samples_indexes[i]];
          for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
          {
              dval += linemodel->getModelFluxDerivContinuumAmpEltVal(elts_indexes[iElt], samples_indexes[i]);
          }
          dval = dval *normFactor * 2 * gsl_vector_get (x,controller->getIndContinuumAmp());
          gsl_matrix_set (J, i, controller->getIndContinuumAmp(),dval);
        }

        if(controller->isRedshiftFitted()){
          Float64 dvalContinuum = linemodel->getModelFluxDerivZContinuumVal(samples_indexes[i] );
          Float64 dval = dvalContinuum;
          for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
          {
            for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
            {
              dval += linemodel-> getModelFluxDerivZEltVal(elts_indexes[iElt], samples_indexes[i], dvalContinuum);
            }
          }
          dval = dval * normFactor;
          gsl_matrix_set (J, i, controller->getIndRedshift(),dval);
        }
    }

    return GSL_SUCCESS;
}

}

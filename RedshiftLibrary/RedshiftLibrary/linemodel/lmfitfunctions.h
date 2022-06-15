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
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/linemodel/lmfitcontroller.h"

#include "RedshiftLibrary/log/log.h"

namespace NSEpic {
/**
 This class is use by the solver lmfit, it is compose of two callback functions
 the first one (lmfit_f) return a vector : for each lamdbda the difference
between model and flux the second one (lmfit_df) return a jacobian matrix : for
each lambda the derivatives amplitudes between model and flux (emission and
absorbtion lines)
**/
class CLineModelFitting;

struct lmfitdata {
  size_t n;
  Float64 *y;
  CLineModelFitting *linemodel;
  TInt32List linemodel_samples_indexes;
  const Float64 *observeGridContinuumFlux;
  CLmfitController *controller;
};

int lmfit_f(const gsl_vector *x, void *data, gsl_vector *f) {
  size_t n = ((struct lmfitdata *)data)->n;
  Float64 *y = ((struct lmfitdata *)data)->y;
  // std::shared_ptr<CLineModelFitting> linemodel = ((struct lmfitdata
  // *)data)->linemodel;
  CLineModelFitting *linemodel = ((struct lmfitdata *)data)->linemodel;
  CLmfitController *controller = ((struct lmfitdata *)data)->controller;
  TInt32List elts_indexes = controller->getFilteredIdx();
  TInt32List samples_indexes =
      ((struct lmfitdata *)data)->linemodel_samples_indexes;
  Float64 normFactor = controller->getNormFactor();
  // const Float64* observeGridContinuumFlux = ((struct lmfitdata
  // *)data)->observeGridContinuumFlux;

  if (controller->isRedshiftFitted()) {
    Float64 redshift = gsl_vector_get(x, controller->getIndRedshift());
    linemodel->setRedshift(redshift, !controller->isNoContinuum());
    if (!controller->isNoContinuum() && !controller->isContinuumFitted()) {
      // a continuum is fit by linemodel , but not by lmfit
      TFloat64List polyCoeffs;
      linemodel->setFitContinuum_tplAmplitude(
          linemodel->getFitContinuum_tplAmplitude(),
          linemodel->getFitContinuum_tplAmplitudeError(), polyCoeffs);
    }
  }

  if (controller->isContinuumFitted()) {
    Int32 idxContinuumTplAmp = controller->getIndContinuumAmp();
    Float64 continuumTplAmp = controller->continuumAmp_LmToModel(
        gsl_vector_get(x, idxContinuumTplAmp)); /// normFactor;
    Float64 continuumTplAmpErr =
        -1.0; // controller->continuumAmpErr_LmToModel(gsl_vector_get
              // (x, idxContinuumTplAmp));
    TFloat64List polyCoeffs;
    linemodel->setFitContinuum_tplAmplitude(continuumTplAmp, continuumTplAmpErr,
                                            polyCoeffs);
  }

  for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++) {
    Float64 amp =
        controller->lineAmp_LmToModel(gsl_vector_get(x, iElt)); /// normFactor;
    linemodel->m_Elements.SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
  }

  if (controller->isEmissionVelocityFitted()) {
    Float64 velocity = controller->emiVel_LmToModel(
        gsl_vector_get(x, controller->getIndEmissionVel()));
    // Log.LogInfo("velocity emission val : %f", velocity);
    linemodel->SetVelocityEmission(velocity);
  }
  if (controller->isAbsorptionVelocityFitted()) {
    Float64 velocity = controller->absVel_LmToModel(
        gsl_vector_get(x, controller->getIndAbsorptionVel()));
    linemodel->SetVelocityAbsorption(velocity);
  }

  //
  if (controller->isContinuumFitted()) {
    linemodel->refreshModelInitAllGrid();
  } else {
    linemodel->refreshModelUnderElements(elts_indexes);
  }

  for (Int32 i = 0; i < n; i++) {
    Float64 Yi = linemodel->getModelFluxVal(samples_indexes[i]) * normFactor;
    gsl_vector_set(f, i, Yi - y[i]);
  }

  return GSL_SUCCESS;
}

int lmfit_df(const gsl_vector *x, void *data, gsl_matrix *J) {
  size_t n = ((struct lmfitdata *)data)->n;
  // std::shared_ptr<CLineModelFitting> linemodel = ((struct lmfitdata
  // *)data)->linemodel;
  CLineModelFitting *linemodel = ((struct lmfitdata *)data)->linemodel;
  CLmfitController *controller = ((struct lmfitdata *)data)->controller;
  TInt32List elts_indexes = controller->getFilteredIdx();
  TInt32List samples_indexes =
      ((struct lmfitdata *)data)->linemodel_samples_indexes;
  Float64 normFactor = controller->getNormFactor();
  const Float64 *observeGridContinuumFlux =
      ((struct lmfitdata *)data)->observeGridContinuumFlux;

  if (controller->isRedshiftFitted()) {
    Float64 redshift = gsl_vector_get(x, controller->getIndRedshift());
    // Log.LogInfo("redshift value %f", redshift);
    linemodel->setRedshift(redshift, false);
    if (!controller->isNoContinuum() && !controller->isContinuumFitted()) {
      // a continuum is fit by linemodel , but not by lmfit
      TFloat64List polyCoeffs;
      linemodel->setFitContinuum_tplAmplitude(
          linemodel->getFitContinuum_tplAmplitude(),
          linemodel->getFitContinuum_tplAmplitudeError(), polyCoeffs);
    }
  }

  if (controller->isContinuumFitted()) {
    Int32 idxContinuumTplAmp = controller->getIndContinuumAmp();
    Float64 continuumTplAmp = controller->continuumAmp_LmToModel(
        gsl_vector_get(x, idxContinuumTplAmp)); /// normFactor;
    Float64 continuumTplAmpErr =
        -1.0; // controller->continuumAmpErr_LmToModel(gsl_vector_get
              // (x, idxContinuumTplAmp));
    TFloat64List polyCoeffs;
    linemodel->setFitContinuum_tplAmplitude(continuumTplAmp, continuumTplAmpErr,
                                            polyCoeffs);
  }

  for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++) {
    Float64 amp =
        controller->lineAmp_LmToModel(gsl_vector_get(x, iElt)); /// normFactor;
    linemodel->m_Elements.SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
  }

  if (controller->isEmissionVelocityFitted()) {
    Float64 velocity = controller->emiVel_LmToModel(
        gsl_vector_get(x, controller->getIndEmissionVel()));
    // Log.LogInfo("velocity emission val : %f", velocity);
    linemodel->SetVelocityEmission(velocity);
  }
  if (controller->isAbsorptionVelocityFitted()) {
    Float64 velocity = controller->absVel_LmToModel(
        gsl_vector_get(x, controller->getIndAbsorptionVel()));
    linemodel->SetVelocityAbsorption(velocity);
  }

  if (controller->isAbsorptionVelocityFitted() &&
      controller->isEmissionVelocityFitted()) {
    linemodel->refreshModelDerivVelUnderElements(elts_indexes);
  } else if (controller->isAbsorptionVelocityFitted()) {
    linemodel->refreshModelDerivVelAbsorptionUnderElements(elts_indexes);
  } else if (controller->isEmissionVelocityFitted()) {
    linemodel->refreshModelDerivVelEmissionUnderElements(elts_indexes);
  }
  // linemodel->refreshModel();
  Float64 normAmpLine = controller->getNormAmpLine();
  for (Int32 i = 0; i < n; i++) {

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++) {
      Float64 dval = linemodel->getModelFluxDerivEltVal(elts_indexes[iElt],
                                                        samples_indexes[i]) *
                     normFactor / normAmpLine * 2 * gsl_vector_get(x, iElt);
      gsl_matrix_set(J, i, iElt, dval);
    }

    if (controller->isEmissionVelocityFitted()) {
      Float64 normEmiFactor = controller->getNormEmiFactor();
      Float64 dval =
          linemodel->getModelFluxDerivVelEmissionVal(samples_indexes[i]) *
          normFactor / normEmiFactor * 2 *
          gsl_vector_get(x, controller->getIndEmissionVel());
      // Log.LogInfo("Deriv sigma Emi for i %d = %f",i, dval);
      gsl_matrix_set(J, i, controller->getIndEmissionVel(), dval);
    }

    if (controller->isAbsorptionVelocityFitted()) {
      Float64 normAbsFactor = controller->getNormAbsFactor();
      Float64 dval =
          linemodel->getModelFluxDerivVelAbsorptionVal(samples_indexes[i]) *
          normFactor / normAbsFactor * 2 *
          gsl_vector_get(x, controller->getIndAbsorptionVel());
      // Log.LogInfo("Deriv sigma Abs for i %d = %f",i, dval);
      gsl_matrix_set(J, i, controller->getIndAbsorptionVel(), dval);
    }

    if (controller->isContinuumFitted()) {
      Float64 dval = observeGridContinuumFlux[samples_indexes[i]];
      for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++) {
        dval += linemodel->getModelFluxDerivContinuumAmpEltVal(
            elts_indexes[iElt], samples_indexes[i]);
      }
      dval = dval * normFactor * 2 *
             gsl_vector_get(x, controller->getIndContinuumAmp());
      gsl_matrix_set(J, i, controller->getIndContinuumAmp(), dval);
    }

    if (controller->isRedshiftFitted()) {
      if (!controller->isNoContinuum()) {
        Float64 dvalContinuum =
            linemodel->getModelFluxDerivZContinuumVal(samples_indexes[i]);
        Float64 dval = dvalContinuum;
        for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++) {
          dval += linemodel->getModelFluxDerivZEltVal(
              elts_indexes[iElt], samples_indexes[i], dvalContinuum);
        }
        dval = dval * normFactor;
        gsl_matrix_set(J, i, controller->getIndRedshift(), dval);
      } else {
        Float64 dval = 0.;
        for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++) {
          dval += linemodel->getModelFluxDerivZEltValNoContinuum(
              elts_indexes[iElt], samples_indexes[i]);
        }
        dval = dval * normFactor;
        gsl_matrix_set(J, i, controller->getIndRedshift(), dval);
      }
    }
  }
  /*
  // export for debug
  FILE* fspc = fopen( "model_flux.txt", "w+" );
  for (Int32 i = 0; i < n; i++)
  {
      fprintf( fspc, "%d %f\n", samples_indexes[i],
  linemodel->getModelFluxVal(samples_indexes[i])*normFactor);
  }
  fclose( fspc );
  //*/

  /*
  // export for debug
  fspc = fopen( "model_derivsigmaAbs.txt", "w+" );
  for (Int32 i = 0; i < n; i++)
  {
      Float64 normAbsFactor = controller->getNormAbsFactor();
      fprintf( fspc, "%d %f\n",
  samples_indexes[i],linemodel->getModelFluxDerivVelAbsorptionVal(samples_indexes[i])*normFactor/normAbsFactor*2*gsl_vector_get
  (x, controller->getIndAbsorptionVel()));
  }
  fclose( fspc );
  //*/

  return GSL_SUCCESS;
}

} // namespace NSEpic

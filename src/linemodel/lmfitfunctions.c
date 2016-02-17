
#include <epic/redshift/linemodel/elementlist.h>

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

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
    {
        Float64 amp = gsl_vector_get (x, iElt)/normFactor;
        linemodel->SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
    }

    Int32 idxVelocityEmission = elts_indexes.size();
    Float64 velocityEmission = gsl_vector_get (x, idxVelocityEmission);
    linemodel->SetVelocityEmission(velocityEmission);
    linemodel->refreshModelUnderElements(elts_indexes);

    for (Int32 i = 0; i < n; i++)
    {
        Float64 Yi = linemodel->getModelFluxVal(samples_indexes[i])*normFactor;
        gsl_vector_set (f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}

//int
//lmfit_f (const gsl_vector * x, void *data,
//         gsl_vector * f)
//{
//    size_t n = ((struct lmfitdata *)data)->n;
//    double *y = ((struct lmfitdata *)data)->y;

//    double A = gsl_vector_get (x, 0);
//    double lambda = gsl_vector_get (x, 1);
//    double b = gsl_vector_get (x, 2);

//    size_t i;

//    for (i = 0; i < n; i++)
//    {
//        /* Model Yi = A * exp(-lambda * i) + b */
//        double t = i;
//        double Yi = A * exp (-lambda * t) + b;
//        gsl_vector_set (f, i, Yi - y[i]);
//    }

//    return GSL_SUCCESS;
//}

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

    for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
    {
        Float64 amp = gsl_vector_get (x, iElt)/normFactor;
        linemodel->SetElementAmplitude(elts_indexes[iElt], amp, 0.0);
    }

    Int32 idxVelocityEmission = elts_indexes.size();
    Float64 velocityEmission = gsl_vector_get (x, idxVelocityEmission);
    linemodel->SetVelocityEmission(velocityEmission);
    linemodel->refreshModelDerivSigmaUnderElements(elts_indexes);

    for (Int32 i = 0; i < n; i++)
    {
        for (Int32 iElt = 0; iElt < elts_indexes.size(); iElt++)
        {
            Float64 dval = linemodel->getModelFluxDerivEltVal(elts_indexes[iElt], samples_indexes[i]);
            gsl_matrix_set (J, i, iElt, dval);
        }
        Float64 dval = linemodel->getModelFluxDerivSigmaVal(samples_indexes[i])*normFactor;
        gsl_matrix_set (J, i, idxVelocityEmission, dval);
    }

    return GSL_SUCCESS;
}

//int
//lmfit_df (const gsl_vector * x, void *data,
//          gsl_matrix * J)
//{
//    size_t n = ((struct lmfitdata *)data)->n;

//    double A = gsl_vector_get (x, 0);
//    double lambda = gsl_vector_get (x, 1);

//    size_t i;

//    for (i = 0; i < n; i++)
//    {
//        /* Jacobian matrix J(i,j) = dfi / dxj, */
//        /* where fi = (Yi - yi)/sigma[i],      */
//        /*       Yi = A * exp(-lambda * i) + b  */
//        /* and the xj are the parameters (A,lambda,b) */
//        double t = i;
//        double e = exp(-lambda * t);
//        gsl_matrix_set (J, i, 0, e);
//        gsl_matrix_set (J, i, 1, -t * A * e);
//        gsl_matrix_set (J, i, 2, 1.0);
//    }
//    return GSL_SUCCESS;
//}

}

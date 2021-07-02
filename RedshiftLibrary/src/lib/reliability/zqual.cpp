#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/processflow/result.h"

#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/datastore.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/extremum/extremum.h"

#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/reliability/pdfzFeatureResult.h"
#include "RedshiftLibrary/reliability/pdfzPredictResult.h"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <stdio.h>
#include <vector>

#include "RedshiftLibrary/reliability/zclassifierstore.h"
#include "RedshiftLibrary/reliability/zqual.h"
#include "RedshiftLibrary/reliability/zqualresult.h"

using namespace std;
using namespace NSEpic;

CQualz::CQualz()
{
}

CQualz::~CQualz() {}

std::shared_ptr<const CQualzResult>
CQualz::Compute(CDataStore &resultStore, CClassifierStore &classifierStore,
                const TFloat64Range &redshiftRange, Float64 &redshiftStep)
{
    Bool storeResult = false;
    CDataStore::CAutoScope resultScope(resultStore, "zReliability/result");

    Log.LogDetail("  ZClassifier: Solving with zstep = %f", redshiftStep);
    storeResult =
        Solve(resultStore, classifierStore, redshiftRange, redshiftStep);

    if (storeResult)
    {
        std::shared_ptr<const CQualzResult> SolveResult =
            (std::shared_ptr<const CQualzResult>)new CQualzResult();
        return SolveResult;
    }

    return NULL;
}

Bool CQualz::Solve(CDataStore &resultStore, CClassifierStore &classifierStore,
                   const TFloat64Range &redshiftRange, Float64 &redshiftStep)
{
    boost::posix_time::ptime startTime;

    /* **********************************************************************
     *              Z-FEATURES IN POSTERIOR PDF
     * ********************************************************************** */
    {
        Log.LogDetail("  ZClassifier: Extracting Features from PDF");
        startTime = boost::posix_time::microsec_clock::local_time();
        ExtractFeaturesPDF(resultStore, redshiftRange, redshiftStep);
        T0 = boost::posix_time::microsec_clock::local_time() - startTime;
    }

    {
        /* **********************************************************************
         *              Z-PROJECTION
         * **********************************************************************
         */
        {
            Log.LogDetail("  ZClassifier: Projecting");
            startTime = boost::posix_time::microsec_clock::local_time();
            ProjectPDF(classifierStore);
            T1 = boost::posix_time::microsec_clock::local_time() - startTime;
        }

        /* **********************************************************************
         *              STORE & DISPLAY
         * **********************************************************************
         */
        {
            Log.LogDetail("  ZClassifier: Displaying prediction");
            // DisplayPrediction();
            Log.LogDetail("  ZClassifier: PREDICTED CLASS zREL : %s", m_predLabel.c_str());
            Log.LogDetail("  ZClassifier: \t*associated posterior class prediction = : %f", m_predProba);

            Log.LogDetail("  ZClassifier: Storing results");
            auto zQual_vect = std::shared_ptr<CPdfzPredictResult>(new CPdfzPredictResult());
            zQual_vect->m_predLabel = m_predLabel;
            zQual_vect->m_score = m_score;
            zQual_vect->m_posterior = m_posterior;
            zQual_vect->m_Labels = classifierStore.GetLabels();
            zQual_vect->m_predProba = m_predProba;

            resultStore.StoreScopedGlobalResult("zpredict", zQual_vect);
        }
    }

    return true;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute descriptors of the redshift posterior PDF P (
 * z | D, I )
 *  ------------------------------------------------------------------------------------------------------
 */
Bool CQualz::ExtractFeaturesPDF(CDataStore &resultStore,
                                const TFloat64Range &redshiftRange,
                                Float64 &redshiftStep)
{
    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results = resultStore.GetGlobalResult(scope_res.c_str());
    auto logzpdf1d =
        std::dynamic_pointer_cast<const CPdfMargZLogResult>(results.lock());

    if (!logzpdf1d)
    {
        Log.LogError("ExtractFeaturesPDF: no results retrieved from scope: %s",
                     scope_res.c_str());
        throw runtime_error("ExtractFeaturesPDF: no results retrieved from scope");
    }

    Log.LogDebug(
        "  ZClassifier: ExtractFeaturesPDF - Checking PDF for Nan vector");
    nanVector = CheckPDF(logzpdf1d->valProbaLog);

    Log.LogDebug("  ZClassifier: ExtractFeaturesPDF - more PDF checks");
    if (logzpdf1d->Redshifts.size() != logzpdf1d->valProbaLog.size())
    {
        throw runtime_error("ZClassifier: ExtractFeaturesPDF: logzpdf1d->Redshifts.size() != logzpdf1d->valProbaLog.size()");
    }

    SPoint dist_peaks;
    Int32 significant_peaks = 0;
    Float64 lbins = log(redshiftStep);
    Float64 zmap, pzmap, cr_dz, cr_nbz, cr_cumpz, r2_cumpz;
    Float64 zmean, zdispersion, zskewness, zskewness2, zkurtosis, loc_z,
        loc_fracz;

    if (nanVector)
    {
        Log.LogError("ZClassifier: ExtractFeaturesPDF: valProbaLog is a nanVector");
        throw runtime_error("ZClassifier: ExtractFeaturesPDF: nanVector");
    } else
    {
        Log.LogDebug("  ZClassifier: ExtractFeaturesPDF - estimate zpdf");
        Int32 nb_pz = logzpdf1d->Redshifts.size();
        TFloat64List zpdf1d;
        zpdf1d.resize(nb_pz);
        zpdf1d[0] = exp(logzpdf1d->valProbaLog[0] + lbins);
        for (Int32 i = 1; i < nb_pz; i++)
        {
            zpdf1d[i] = exp(logzpdf1d->valProbaLog[i] + lbins);
        }
        /* **********************************************************************
         *              FIND EXTREMA IN ZPDF
         * **********************************************************************
         */
        Log.LogDebug("  ZClassifier: ExtractFeaturesPDF - find extrema");
        TPointList extremumList;
        Int32 extremumCount = 1000;
        Float64 radius = 0.005;

        CExtremum extremum( extremumCount, 2*radius); //no peak sep, no cut
        extremum.Find(logzpdf1d->Redshifts, logzpdf1d->valProbaLog,
                      extremumList); // already sorted list in logPZ in output
        // extremum.Find( logzpdf1d->Redshifts, zpdf1d, extremumList ); //
        // overflow or underflow of exponential values...
        // std::cout<<"extremumList.size() = " <<
        // extremumList.size()<<std::endl;

        Bool noExtremum = false;
        if (extremumList.size() == 0)
        {
            dist_peaks.X = NAN; // distance in z between the two first best peaks in zPDF
            dist_peaks.Y = NAN; // distance in pz
            significant_peaks = -1;
            noExtremum = true;
        } else if (extremumList.size() == 1)
        {
            dist_peaks.X = 0; // distance in z between the two first best peaks in zPDF
            dist_peaks.Y = 1; // distance in pz
            significant_peaks = 1;
        } else
        {
            // distance in z between the two first best peaks in zPDF
            dist_peaks.X = abs(extremumList[0].X - extremumList[1].X);
            // distance in pz
            dist_peaks.Y = (exp(extremumList[0].Y + lbins) - exp(extremumList[1].Y + lbins));
        }

        if (noExtremum)
        {
            Log.LogError("ZClassifier: ExtractFeaturesPDF: %d extrema found", extremumList.size());
            throw runtime_error("ZClassifier: ExtractFeaturesPDF: no extremum found");
        } else
        {

            /* **********************************************************************
             *              MAXIMUM A POSTERIORI Z-ESTIMATE
             * ***********************************************************************
             */
            Log.LogDebug("  ZClassifier: ExtractFeaturesPDF - find map");
            Int32 indice_map = 0;
            zmap = logzpdf1d->Redshifts[indice_map];
            pzmap = exp(logzpdf1d->valProbaLog[indice_map] + lbins);
            SPoint map_estimate;
            map_estimate.X = zmap;
            map_estimate.Y = pzmap;

            for (Int32 i = 1; i < nb_pz; i++)
            {
                if (zpdf1d[i] >= pzmap)
                {
                    indice_map = i;
                    pzmap = zpdf1d[indice_map];
                    zmap = logzpdf1d->Redshifts[indice_map];
                }
            }
            Int32 N2 = (redshiftRange.GetEnd() - redshiftRange.GetBegin()) /
                           redshiftStep +
                       1;
            Bool condi = (N2 == nb_pz);
            if (condi)
            {
                loc_z = 0;
            } else
            {
                loc_z = 1; // = 1 if constant prior
            }
            loc_fracz = 100 * nb_pz / (N2);

            significant_peaks = GetNPeaksKM(extremumList, zmap, lbins);

            /* **********************************************************************
             *              GLOBAL STATISTICS:
             *
             *              MEAN mu                         = sum [ p(z_i) * z_i
             * ] VARIANCE sigma2         = sum [ p(z_i) * (z_i - mu)^2 ]
             *              DISPERSION sigma        = sqrt [ sigma2 ]
             *
             *              CENTERED MOMENTS, (rank j)
             *                      m_j =  sum [ p(z_i) * (   ( z_i - mu ) /
             * sigma  )^j   ]
             *
             *              SKEWNESS s      = sum [ p(z_i) * (  ( z_i - mu ) /
             * sigma  )^3   ] = m_3 / sigma^3 = [ E[z^3] - 3*mu*sigma2 - mu^3 ]
             * / sigma.^3
             *
             *              KURTOSIS k   = sum [ p(z_i) * (  ( z_i - mu ) /
             * sigma  )^4   ] = m_4 / sigma2^2
             *
             * **********************************************************************
             */
            Float64 m3, zcenter, zmoment1, zmoment2, zmoment3, zmoment4;
            m3 = 0;
            zmoment1 = 0;
            for (Int32 i = 0; i < nb_pz; i++)
            {
                zmoment1 += zpdf1d[i] * logzpdf1d->Redshifts[i]; // = E[z]
                m3 += zpdf1d[i] * logzpdf1d->Redshifts[i] *
                      logzpdf1d->Redshifts[i] *
                      logzpdf1d->Redshifts[i]; // = E[z^3]
            }

            zmoment2 = 0;
            zmoment3 = 0;
            zmoment4 = 0;
            for (Int32 i = 0; i < nb_pz; i++)
            {
                zcenter = (logzpdf1d->Redshifts[i] - zmoment1);
                zmoment2 +=
                    zpdf1d[i] * zcenter * zcenter; // = E [ (z - E[z] ) ^2 ]
                zmoment3 += zpdf1d[i] * zcenter * zcenter *
                            zcenter; // = E [ (z - E[z] ) ^3 ]
                zmoment4 += zpdf1d[i] * zcenter * zcenter * zcenter *
                            zcenter; // = E [ (z - E[z] ) ^4 ]
            }
            Float64 mu = zmoment1, mu3 = zmoment1 * zmoment1 * zmoment1;
            Float64 sigma = sqrt(zmoment2), sigma2 = zmoment2,
                    sigma3 = sigma2 * sigma, sigma4 = sigma2 * sigma2;

            zmean = zmoment1;
            zdispersion = sqrt(zmoment2);
            zskewness = (m3 - 3 * mu * sigma2 - mu3) / sigma3;
            zskewness2 = zmoment3 / sigma3;
            zkurtosis = zmoment4 / sigma4;

            /* **********************************************************************
             *              CREDIBILITY REGION (95%) characteristics
             * ***********************************************************************
             */
            Log.LogDebug(
                "  ZClassifier: ExtractFeaturesPDF - estimate CR95 features");
            cr_dz = 0;
            cr_nbz = 1;
            cr_cumpz = pzmap;
            Float64 cr_zleft = zmap, cr_zright = zmap;
            Int32 indice_left = indice_map, indice_right = indice_map;

            Int32 isearchMethod = 2; // 2 is closer to Matlab's behaviour, 0 is
                                     // the initial C++ implementation
            Bool verboseCRIndexSearch = false;
            if (isearchMethod == 2)
            {
                Float64 criter = 1e-3;
                Int32 maxiborne = indice_map;
                Int32 maxiborneLimit = int(nb_pz - indice_map - 1);
                while (maxiborne < maxiborneLimit &&
                       logzpdf1d->Redshifts[maxiborne] < zmap + criter)
                {
                    maxiborne++;
                }
                if (maxiborne > indice_map)
                {
                    maxiborne -= 1;
                }
                Int32 margin = 10;
                maxiborne = maxiborne + margin;
                Bool stop = false;
                for (Int32 kr = 0; kr < maxiborne - indice_map; kr++)
                {
                    for (Int32 kl = 0; kl <= kr + 3; kl++)
                    {
                        indice_left = indice_map - kl;
                        indice_right = indice_map + kr;
                        cr_zleft = logzpdf1d->Redshifts[indice_left];
                        cr_zright = logzpdf1d->Redshifts[indice_right];
                        cr_cumpz = 0.0;
                        cr_nbz = 0;
                        for (Int32 k = indice_left; k <= indice_right; k++)
                        {
                            cr_nbz++;
                            cr_cumpz += zpdf1d[k];
                        }
                        if (cr_cumpz >= 0.95)
                        {
                            stop = true;
                        }
                        if (verboseCRIndexSearch)
                        {
                            std::cout
                                << " \t Idx search: "
                                << "\n"
                                << " \t indice_left: " << indice_left << "\n"
                                << " \t indice_right: " << indice_right << "\n"
                                << " \t cr_cumpz: " << cr_cumpz << "\n"
                                << " \t cr_nbz: " << cr_nbz << "\n"
                                << " \n";
                        }
                        if (stop)
                        {
                            break;
                        }
                    }
                    if (stop)
                    {
                        break;
                    }
                }
            } else if (isearchMethod == 1)
            {
                // right first, with criter 2e-3
                Float64 criter = 1e-3;
                Int32 maxiborne = indice_map;
                Int32 maxiborneLimit = int(nb_pz - indice_map);
                while (logzpdf1d->Redshifts[maxiborne] < zmap + criter &&
                       maxiborne < maxiborneLimit)
                {
                    maxiborne++;
                }
                Int32 margin = 10;
                while (cr_cumpz <= 0.95 && indice_right < (maxiborne + margin))
                {
                    indice_right++;
                    if (indice_right < nb_pz - 1)
                    {
                        cr_nbz++;
                        cr_cumpz += zpdf1d[indice_right];
                        cr_zright = logzpdf1d->Redshifts[indice_right];
                        if (cr_cumpz <= 0.95)
                        {
                            if (indice_left > 0)
                            {
                                indice_left--;
                                cr_nbz++;
                                cr_cumpz += zpdf1d[indice_left];
                                cr_zleft = logzpdf1d->Redshifts[indice_left];
                            }
                        }
                    }
                    if (verboseCRIndexSearch)
                    {
                        std::cout << " \t Idx search: "
                                  << "\n"
                                  << " \t indice_right: " << indice_right
                                  << "\n"
                                  << " \t indice_left: " << indice_left << "\n"
                                  << " \t cr_cumpz: " << cr_cumpz << "\n"
                                  << " \t cr_nbz: " << cr_nbz << "\n"
                                  << " \n";
                    }
                }
            } else
            {
                // Original Method : left first
                Bool cond_stop = false;
                Int32 crit_stop = round((nb_pz * 0.5) * 0.5);
                while (cr_cumpz <= 0.95 && indice_left >= 0 && !cond_stop)
                {
                    indice_left--;
                    if (indice_left >= (max)(0, indice_map - crit_stop))
                    {
                        cr_nbz++;
                        cr_cumpz += zpdf1d[indice_left];
                        cr_zleft = logzpdf1d->Redshifts[indice_left];
                        if (cr_cumpz <= 0.95)
                        {
                            indice_right++;
                            if (indice_right < nb_pz - 1)
                            {
                                cr_nbz++;
                                cr_cumpz += zpdf1d[indice_right];
                                cr_zright = logzpdf1d->Redshifts[indice_right];
                            }
                        }
                    }
                }
            }
            cr_dz = cr_zright -
                    cr_zleft; //+1*redshiftStep; //do not add redshiftStep to
                              //reproduce Matlab results best
            Int32 cr_idx_zleft = indice_left;
            Int32 cr_idx_zright = indice_right;
            /* **********************************************************************
             *              REGION R2 ( dz = 2 * espsilon ) characteristics
             * **********************************************************************
             */
            Log.LogDebug(
                "  ZClassifier: ExtractFeaturesPDF - estimate CR_2 features");
            Float64 zshift = 1e-3;
            Int32 methodR2Indexes = 1; // 0= initial method, 1= 20180417 method
                                       // that reproduces best Matlab
            if (methodR2Indexes == 0)
            {
                Int32 epsilon = round(zshift / redshiftStep);
                indice_left = (max)(0, indice_map - epsilon);
                indice_right = (min)(indice_map + epsilon, nb_pz - 1);
            } else if (methodR2Indexes == 1)
            {
                indice_left = indice_map;
                Int32 miniLimit = int(0);
                while (logzpdf1d->Redshifts[indice_left] >
                           zmap - zshift + redshiftStep &&
                       indice_left > miniLimit)
                {
                    indice_left--;
                }
                indice_right = indice_map;
                Int32 maxiLimit = int(nb_pz - 1);
                while (indice_right < maxiLimit &&
                       logzpdf1d->Redshifts[indice_right] < zmap + zshift)
                {
                    indice_right++;
                }
            }
            r2_cumpz = 0.0;
            for (int i = indice_left; i <= indice_right; i++)
            {
                r2_cumpz += zpdf1d[i];
            }
            Int32 cr2_idx_zleft = indice_left;
            Int32 cr2_idx_zright = indice_right;

            Bool showDescriptors = false;
            if (showDescriptors)
            {
                if (0)
                {
                    zdispersion = 0.000632706;
                    cr_cumpz = 0.95296;
                    cr_dz = 0.0015;
                    cr_nbz = 16.0000e+00;
                    r2_cumpz = 0.992187;
                }
                if (0)
                {
                    cr_cumpz = 0.959401987918265;
                    cr_dz = 0.000699999999999923;
                    cr_nbz = 8.0000e+00;
                }
                std::setprecision(20);
                std::cout << "-------------------------------------------------"
                             "-----------------"

                          << "CHECK VERSUS MATLAB FEATURES"
                          << " \t \n"
                          << " \t [STAT] dispersion = " << zdispersion << "\n"
                          << " \t [MAP] pzmap = " << pzmap << "\n"

                          << " \t [CR 95%] nz = " << cr_nbz << "\n"
                          << " \t [CR 95%] dz = " << cr_dz << "\n"
                          << " \t [CR 95%] pz = " << cr_cumpz << "\n"
                          << " \t -1\n"
                          << " \t [R2 1E-3] pz = " << r2_cumpz << "\n"
                          << " \t [KM] dist_pz = " << dist_peaks.Y << "\n"
                          << " \t [KM] #modes = " << significant_peaks << "\n"
                          << "-------------------------------------------------"
                             "-----------------"
                          << "\n"
                          << " \t   zmap = " << zmap << "\n"
                          << " \t   indice_map = " << indice_map << "\n"
                          << " \t   cr_zleft = " << cr_zleft << "\n"
                          << " \t   cr_zright = " << cr_zright << "\n"
                          << " \t   cr_idx_zleft = " << cr_idx_zleft << "\n"
                          << " \t   cr_idx_zright = " << cr_idx_zright << "\n"
                          << " \t   zmoment2 = " << zmoment2 << "\n"
                          << " \t   cr2_idx_zleft = " << cr2_idx_zleft << "\n"
                          << " \t   cr2_idx_zright = " << cr2_idx_zright << "\n"
                          << "\n";
            }

            /* **********************************************************************
             *              QUICK DISPLAY ON CONSOLE
             * **********************************************************************
             */
            if (0)
            {
                std::setprecision(20);
                std::cout << "-------------------------------------------------"
                             "-----------------"
                          << "\n"
                          << " \t [MAP] i = " << indice_map << "\n"
                          << " \t [MAP] zmap = " << zmap << "\n"
                          << " \t [MAP] pzmap = " << pzmap << "\n"
                          << " \t [MAP] loc_z = " << loc_z << "\n"
                          << " \t [MAP] loc_fracz = " << loc_fracz << "\n"
                          << "-------------------------------------------"
                          << "\n"
                          << " \t [STAT] mean = " << zmean << "\n"
                          << " \t [STAT] dispersion = " << zdispersion << "\n"
                          << " \t [STAT] skewness = " << zskewness << "\n"
                          << " \t [STAT] skewness2 = " << zskewness2 << "\n"
                          << " \t [STAT] kurtosis = " << zkurtosis << "\n"
                          << "-------------------------------------------"
                          << "\n"
                          << " \t [CR 95%] nz = " << cr_nbz << "\n"
                          << " \t [CR 95%] dz = " << cr_dz << "\n"
                          << " \t [CR 95%] pz = " << cr_cumpz << "\n"
                          << "-------------------------------------------"
                          << "\n"
                          //<< " \t [R2 1E-3] dz = " << r2_nbz << "\n"
                          //<< " \t [R2 1E-3] z_left = " << r2_zleft << "\n"
                          //<< " \t [R2 1E-3] z_right = " << r2_zright
                          << " \t [R2 1E-3] pz = " << r2_cumpz << "\n"
                          << "-------------------------------------------"
                          << "\n"
                          << " \t [KM] #modes = " << significant_peaks << "\n"
                          << " \t [KM] dist_z = " << dist_peaks.X << "\n"
                          << " \t [KM] dist_pz = " << dist_peaks.Y << "\n"
                          << " \t [KM] 1st Peak (z=" << extremumList[0].X
                          << " ; pz= " << exp(extremumList[0].Y + lbins)
                          << ") \n"
                          << " \t [KM] 2nd Peak (z=" << extremumList[1].X
                          << " ; pz= " << exp(extremumList[1].Y + lbins) << ") "
                          << std::endl;
            }
        }
    }

    /* **********************************************************************
     *              STORE FEATURE Vector
     * ********************************************************************** */
    Log.LogDebug("  ZClassifier: ExtractFeaturesPDF - store features");
    std::string desc_1 =
        "Localized_prior"; // boolean (0 if constant prior; 1 otherwise)
    std::string desc_2 =
        "%z candidats"; // if localized prior on z-candidates, fraction in the
                        // initial zRange (in %)
    std::string desc_3 = "std_p(z|D,I) ";
    std::string desc_4 = "skewness_p(z|D,I)";
    std::string desc_5 = "kurtosis_p(z|D,I)";
    std::string desc_6 = "p(zmap_1d|D,I)";
    std::string desc_7 =
        "#z_in_R1"; // Region around zMAP [R1], credibility region (CR) 95% ;
    std::string desc_8 = "dz_in_R1";
    std::string desc_9 = "sum_p(z_in_R1|D,I)";
    std::string desc_10 =
        "sum_p(z_in_R2|D,I)"; // Region around zMAP [R2], crit dz=0.001;
    std::string desc_11 = "dist_(2peaks)_in_dz";
    std::string desc_12 = "dist_(2peaks_in_proba";
    std::string desc_13 = "#km_peaks_1d";

    /*
    if ( false) { //test==true) {//nanVector==1 ) {
            std::cout << " \t [STAT] test = " << test << std::endl;
            loc_z =0;//NAN;
            loc_fracz =0;// NAN;
            zdispersion =0;//NAN;
            zskewness =0;//NAN;
            zkurtosis =0;//NAN;
            pzmap =0;//NAN;
            cr_nbz =0;//NAN;
            cr_dz =0;//NAN;
            cr_cumpz =0;//NAN;
            r2_cumpz =0;//NAN;
            dist_peaks.X =0;//NAN;
            dist_peaks.Y =0;//NAN;
            significant_peaks =0;//NAN;
    }*/

    if (m_doTEST)
    { // .m file
        /*
    'Std P(z|D,I)'
    'P(z_{MAP 1d}|D,I)'
    '#z\in CR_{95%}'
    '\Deltaz\in CR_{95%}'
    'Sum P(z\in CR_{95%} | D, I)'
    '\Deltaz\in Crit_{.1%}'
    'Sum P(z\in Crit_{.1%} | D, I)'
    '#dist Peaks_{pdf}'
    '#km Peak_{pdf} ~MAPz1d'
        */
        zdispersion = 5.8595e-05;
        pzmap = 6.2780e-01;
        cr_nbz = 3.0000e+00;
        cr_dz = 2.0000e-04;
        cr_cumpz = 9.8888e-01;
        // 2.0000e-03 ;
        r2_cumpz = 1.0000e+00;
        dist_peaks.Y = 6.2780e-01;
        significant_peaks = 1.0000e+00;
    }

    /*
            std::string desc_3 = "std_p(z|D,I) "; //2 std::string desc_6 =
       "p(zmap_1d|D,I)";                                          //5
            std::string desc_7 = "#z_in_R1"; //6 std::string desc_8 =
       "dz_in_R1"; //7 std::string desc_9 = "sum_p(z_in_R1|D,I)"; //8
            std::string desc_10 = "sum_p(z_in_R2|D,I)"; //9 std::string desc_12
       = "dist_(2peaks_in_proba";                  //11 std::string desc_13 =
       "#km_peaks_1d";                                           //12

            TInt64List selected = { 2, 5, 6, 7, 8, -1, 9, 11, 12 };
     */

    m_zfeatures.resize(13);
    m_zfeatures[0] = loc_z;
    m_zfeatures[1] = loc_fracz;
    m_zfeatures[2] = zdispersion;
    m_zfeatures[3] = zskewness;
    m_zfeatures[4] = zkurtosis;
    m_zfeatures[5] = pzmap;
    m_zfeatures[6] = cr_nbz;
    m_zfeatures[7] = cr_dz;
    m_zfeatures[8] = cr_cumpz;
    m_zfeatures[9] = r2_cumpz;
    m_zfeatures[10] = dist_peaks.X;
    m_zfeatures[11] = dist_peaks.Y;
    m_zfeatures[12] = significant_peaks;

    id_descriptors.resize(13);
    id_descriptors[0] = desc_1;
    id_descriptors[1] = desc_2;
    id_descriptors[2] = desc_3;
    id_descriptors[3] = desc_4;
    id_descriptors[4] = desc_5;
    id_descriptors[5] = desc_6;
    id_descriptors[6] = desc_7;
    id_descriptors[7] = desc_8;
    id_descriptors[8] = desc_9;
    id_descriptors[9] = desc_10;
    id_descriptors[10] = desc_11;
    id_descriptors[11] = desc_12;
    id_descriptors[12] = desc_13;

    auto zfeature_vect =
        std::shared_ptr<CPdfzFeatureResult>(new CPdfzFeatureResult());

    zfeature_vect->mapzfeatures[desc_1] = loc_z;
    zfeature_vect->mapzfeatures[desc_2] = loc_fracz;
    zfeature_vect->mapzfeatures[desc_3] = zdispersion;
    zfeature_vect->mapzfeatures[desc_4] = zskewness;
    zfeature_vect->mapzfeatures[desc_5] = zkurtosis;
    zfeature_vect->mapzfeatures[desc_6] = pzmap;
    zfeature_vect->mapzfeatures[desc_7] = cr_nbz;
    zfeature_vect->mapzfeatures[desc_8] = cr_dz;
    zfeature_vect->mapzfeatures[desc_9] = cr_cumpz;
    zfeature_vect->mapzfeatures[desc_10] = r2_cumpz;
    zfeature_vect->mapzfeatures[desc_11] = dist_peaks.X;
    zfeature_vect->mapzfeatures[desc_12] = dist_peaks.Y;
    zfeature_vect->mapzfeatures[desc_13] = significant_peaks;

    resultStore.StoreScopedGlobalResult("zfeatures", zfeature_vect);

    return true;
}

void CQualz::DisplayPrediction()
{
    std::cout << "-------------------------------------------------------------"
                 "--------------------------------------------\n"
              << "       PREDICTED CLASS zREL : " << m_predLabel << "\n"
              << " \t* [zQual] associated posterior class prediction = "
              << std::setprecision(10) << m_predProba << std::endl;

    /*
    std::cout<<std::setprecision(10)
    <<"\t* [Class C1] : \t PostProba = "<<m_posterior->data[0] <<" ;  Score = "
    <<m_score->data[0] <<"\n"
    <<"\t* [Class C2] : \t PostProba = "<<m_posterior->data[1] <<" ;  Score = "
    <<m_score->data[1] <<"\n"
    <<"\t* [Class C3] : \t PostProba = "<<m_posterior->data[2] <<" ;  Score = "
    <<m_score->data[2] <<"\n"
    <<"\t* [Class C4] : \t PostProba = "<<m_posterior->data[3] <<" ;  Score = "
    <<m_score->data[3] <<"\n"
    <<"\t* [Class C5] : \t PostProba = "<<m_posterior->data[4] <<" ;  Score = "
    <<m_score->data[4] //<<"\n"
    << std::endl;
    */

    if (disp_time)
    {
        std::cout << "\t* [zQual] extract features (in sec): "
                  << to_simple_string(T0) << "\n"
                  << "\t* [zQual] z-project (in sec): " << to_simple_string(T1)
                  << std::endl;
    }
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute # of significant peaks in the posterior zPDF
 * P ( z | D, I )
 *  ------------------------------------------------------------------------------------------------------
 */
Int32 CQualz::GetNPeaksKM(TPointList &peaks, Float64 &zmap_estimate,
                          Float64 &lbins)
{
    gsl_matrix *data;          // Input data = modes/peaks(zPDF)
    Int32 nClusters = 2;       //# Clusters
    Int32 nObs = peaks.size(); //# Observations
    Int32 nDesc = 1; //# Descriptor(s) (here only the value of the peaks(zPDF)
    Int32 maxIter = 100;       // 500;     // Maximum of allowed iterations
    std::string method = "L2"; // Distance method

    gsl_matrix *centroids; // [Output] Centroids of clusters
    gsl_vector_int *index; // [Output] Index values (from 0, to nClusters-1 )

    Int32 i, j, h, km;
    Float64 *weights;
    gsl_vector *temp1, *temp2, *temp3;
    Bool cond_converge = false;

    /* **********************************************************************
     *              STEP 0 - Allocate memory
     * ********************************************************************** */
    data = gsl_matrix_alloc(nObs, nDesc);
    for (i = 0; i < nObs; i++)
    {
        for (j = 0; j < nDesc; j++)
        {
            gsl_matrix_set(data, i, j, exp(peaks[i].Y + lbins));
        }
    }
    temp1 = gsl_vector_alloc(nDesc);
    temp2 = gsl_vector_alloc(nDesc);
    temp3 = gsl_vector_alloc(nObs);
    index = gsl_vector_int_alloc(nObs);
    centroids = gsl_matrix_alloc(nClusters, nDesc);
    weights = (Float64 *)malloc(sizeof(Float64) * nObs);
    gsl_vector_int *index_init = gsl_vector_int_alloc(nObs);

    /* **********************************************************************
     *              STEP 1 - Initialize centroïds
     * ********************************************************************** */
    Float64 minipz = 0, maxipz = 0;
    minipz = gsl_matrix_min(data);
    maxipz = gsl_matrix_max(data);
    gsl_matrix_minmax(data, &minipz, &maxipz);
    for (i = 0; i < centroids->size1; i++)
    {
        for (j = 0; j < centroids->size2; j++)
        {
            gsl_matrix_set(centroids, i, j,
                           minipz + (((Float64)rand()) / RAND_MAX) *
                                        (maxipz - minipz));
        }
    }
    /* **********************************************************************
     *              STEP 2 - Optimal search
     * ********************************************************************** */
    Int32 iter = 0;
    while (iter <= maxIter && !cond_converge)
    {
        iter++;
        //  STEP 2.1 - Update index
        Float64 dist;
        Float64 mindist;
        Int32 optimIndex = 0;
        for (i = 0; i < nObs; i++)
        {
            mindist = INFINITY;
            gsl_matrix_get_row(temp1, data, i);
            for (j = 0; j < nClusters; j++)
            {
                gsl_matrix_get_row(temp2, centroids, j);
                dist = GetDistanceKM(temp1, temp2, method);
                if (dist < mindist)
                {
                    mindist = dist;
                    optimIndex = j;
                }
            }
            gsl_vector_int_set(index, i, optimIndex);
        }

        // STEP 2.2 - Update centroïds
        Float64 weighted_mean;
        for (i = 0; i < nClusters; i++)
        {
            for (j = 0; j < nObs; j++)
            {
                if (gsl_vector_int_get(index, j) == i)
                {
                    weights[j] = 1.0;
                } else
                {
                    weights[j] = 0.0;
                }
            }

            for (h = 0; h < nDesc; h++)
            {
                gsl_matrix_get_col(temp3, data, h);
                weighted_mean =
                    gsl_stats_wmean(weights, 1, temp3->data, 1, temp3->size);
                gsl_matrix_set(centroids, i, h, weighted_mean);
            }
        }

        // STEP 2.3 - Convergence criteria
        Int32 counter = 0;
        for (i = 0; i < nObs; i++)
        {
            if (index->data[i] == index_init->data[i])
            {
                counter++;
            }
            index_init->data[i] = index->data[i];
        }
        cond_converge = (index->size - counter) ==
                        0; // the algorithm is stuck on the best solution
    }                      // end_while

    /* **********************************************************************
     *              STEP 3 - # Elements in zMap cluster
     * ********************************************************************** */
    TFloat64List nbElements;
    nbElements.resize(nClusters);
    Int32 ic, nelements, zmap_cluster = -1;
    for (ic = 0; ic < nClusters; ic++)
    {
        nelements = 0;
        for (i = 0; i < nObs; i++)
        {
            if (gsl_vector_int_get(index, i) == ic)
            {
                nelements++;
                if (peaks[i].X == zmap_estimate)
                {
                    zmap_cluster = ic; // cluster in where the zmap is located
                    // zmap_position = i;    // very likely to be the first in
                    // this list
                }
            }
        }
        nbElements[ic] = nelements;
        // std::cout<<"nelements in ["<<ic<<"] =" << nelements <<std::endl;
    }

    km = nbElements[zmap_cluster]; // # significant peaks

    /* **********************************************************************
     *              STEP 4 - Free memory
     * ********************************************************************** */
    free(weights);

    gsl_matrix_free(data);
    gsl_vector_int_free(index);
    gsl_matrix_free(centroids);

    gsl_vector_free(temp1);
    gsl_vector_free(temp2);
    gsl_vector_free(temp3);

    return km;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute the distance between two vectors: L1, L2,
 * COSINE, XC
 *  ------------------------------------------------------------------------------------------------------
 */
Float64 CQualz::GetDistanceKM(gsl_vector *vect1, gsl_vector *vect2,
                              std::string method)
{
    Int32 i;
    Float64 distance = 0.0;

    /* CITY-BLOCK DISTANCE */
    if (method == "L1")
    {
        Float64 diff = 0.0;
        for (i = 0; i < vect1->size; i++)
        {
            diff = vect1->data[i] - vect2->data[i];
            distance += abs(diff);
        }
    }

    /* EUCLIDEAN DISTANCE  */
    if (method == "L2")
    {
        Float64 diff = 0.0;
        for (i = 0; i < vect1->size; i++)
        {
            diff = vect1->data[i] - vect2->data[i];
            distance += diff * diff;
        }
    }

    /* COSINE DISTANCE */
    if (method == "cosine")
    {
        Float64 xcorr = 0.0, normx = 0.0, normc = 0.0;
        for (i = 0; i < vect1->size; i++)
        {
            xcorr += vect1->data[i] * vect2->data[i];
            normx += vect1->data[i] * vect1->data[i];
            normc += vect2->data[i] * vect2->data[i];
        }
        distance = 1 - xcorr / sqrt(normx * normc);
    }

    /* CORRELATION */
    if (method == "correlation")
    {
        Float64 center_c = 0.0, center_x = 0.0; // centered (~mean)
        for (i = 0; i < vect1->size; i++)
        {
            center_x += vect1->data[i];
            center_c += vect2->data[i];
        }
        center_x = center_x / vect1->size;
        center_c = center_c /
                   vect2->size; // vectors temp1 and temp2 have the same size

        Float64 xcorr_c = 0.0, normx_c = 0.0, normc_c = 0.0;
        for (i = 0; i < vect1->size; i++)
        {
            xcorr_c +=
                (vect1->data[i] - center_x) * (vect2->data[i] - center_c);
            normx_c += (vect1->data[i] * vect1->data[i] - center_x);
            normc_c += (vect2->data[i] * vect2->data[i] - center_c);
        }
        distance = 1 - xcorr_c / sqrt(normx_c * normc_c);
    }

    return distance;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Check if the full zPDF == NAN
 *                                      This can occur when wrong input flux &
 * noise spectra are processed
 *  ------------------------------------------------------------------------------------------------------
 */
Bool CQualz::CheckPDF(const TFloat64List &zpdf)
{
    Int32 counter = 0;
    for (Int32 i = 1; i < zpdf.size(); i++)
    {
        if (std::isnan(zpdf[i]))
        {
            counter++;
        }
    }
    if (counter == zpdf.size() - 1)
    {
        return true;
    } else
    {
        return false;
    }
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Project the feature vector into a trainde mapping
 *                                              and compute score & posterior
 * class probabilities in addiction to the predicted label
 *  ------------------------------------------------------------------------------------------------------
 */
void CQualz::ProjectPDF(CClassifierStore &classifierStore)
{
    if (nanVector)
    {
        m_predLabel = "";
    } else
    {

        m_predLabel = " ";
        m_score = gsl_vector_alloc(classifierStore.GetNbClasses());
        m_posterior = gsl_vector_alloc(classifierStore.GetNbClasses());

        GetScorePred(classifierStore); // store vector m_score, for all learners
        GetPosteriorPred(
            classifierStore); // store vector m_posterior, for all classes
        GetLabelPred(
            classifierStore); // decode score knowing the coding strategy (ECOC)
                              // to  derive the predicted label & its associated
                              // posterior proba
    }
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Get Predicted zReliability Class
 *                                      by combining the (binary) predictions of
 * each learner + the additional step of 'DECODE'  in ECOC
 *  ------------------------------------------------------------------------------------------------------
 */
void CQualz::GetLabelPred(CClassifierStore &classifierStore)
{
    Int32 opt_predL = 1;

    Int32 ind_min = 0;
    Float64 dist, delta_min = INFINITY;
    Int32 ind_max = 0;
    Float64 delta_max = -INFINITY;

    if (opt_predL == 1)
    {
        // OPTION 1 :   argmin LOSS FUNCTION ( QUADRATIC )
        // distance (class i ) = mean( (1 - getTimes( CODING_MATRIX , 2*SCORE -
        // 1 ) ).^2 ) / 2;
        gsl_vector *delta;
        delta = gsl_vector_alloc(classifierStore.GetNbClasses());

        gsl_vector *sc2;
        sc2 = gsl_vector_alloc(m_score->size);
        for (Int32 id_lrn = 0; id_lrn < classifierStore.GetNbLearners();
             id_lrn++)
        {
            sc2->data[id_lrn] = 2 * m_score->data[id_lrn] - 1;
        }
        gsl_matrix *prodM = GetTimesKL(classifierStore.GetCodingMatrix(), sc2);
        gsl_matrix_scale(prodM, -1);
        gsl_matrix_add_constant(prodM, 1);
        gsl_matrix *prodM_ = gsl_matrix_alloc(prodM->size1, prodM->size2);
        gsl_matrix_memcpy(prodM_, prodM);
        gsl_matrix_mul_elements(prodM_, prodM);

        for (Int32 id_class = 0; id_class < classifierStore.GetNbClasses();
             id_class++)
        {
            dist = 0.0;
            /*
            Float64 s = 0, m_binaryPred = 0;
            for ( Int32 id_lrn=0; id_lrn<classifierStore.GetNbLearners();
            id_lrn++) { double m = gsl_matrix_get(
            classifierStore.GetCodingMatrix(), id_class, id_lrn ); m_binaryPred
            = gsl_vector_alloc( classifierStore.GetNbLearners() ); if (
            m_score->data[id_lrn] <=0) { m_binaryPred = -1; } else {
                            m_binaryPred = +1;
                    }
                    dist = ( m - m_binaryPred );                    // DECODE
            (ECOC) s += dist*dist;
            }
            delta->data[id_class] = sqrt(s);                // minimizing the
            euclidean function
             */
            for (Int32 id_lrn = 0; id_lrn < classifierStore.GetNbLearners();
                 id_lrn++)
            {
                double m = gsl_matrix_get(prodM_, id_class, id_lrn);
                dist += m;
            }
            delta->data[id_class] = dist / 2;
            if (delta->data[id_class] < delta_min)
            {
                ind_min = id_class;
                delta_min = delta->data[id_class];
            }
        }
        m_idpredLabel = ind_min + 1; // the class_id start from 1 -> K classes

        m_predProba = m_posterior->data[ind_min];
        m_predLabel = classifierStore.GetLabel(ind_min);

        gsl_vector_free(delta);
        gsl_vector_free(sc2);
        gsl_matrix_free(prodM);
        gsl_matrix_free(prodM_);

    } else if (opt_predL == 2)
    {
        // OPTION 2:            argmax  posterior probabilities
        for (Int32 id_class = 0; id_class < classifierStore.GetNbClasses();
             id_class++)
        {
            if (m_posterior->data[id_class] > delta_max)
            {
                ind_max = id_class;
                delta_max = m_posterior->data[id_class];
            }
        }
        m_idpredLabel = ind_max + 1; // the class_id start from 1 -> K classes

        m_predProba = m_posterior->data[ind_max];
        m_predLabel = classifierStore.GetLabel(ind_max);
    }
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Get the score for each learner
 *  ------------------------------------------------------------------------------------------------------
 */
void CQualz::GetScorePred(CClassifierStore &classifierStore)
{
    Int32 M, P;

    auto mapLearners = classifierStore.GetLearners();
    CClassifierStore::MapLearners::const_iterator it = mapLearners.begin();
    auto learner = ((it->second));

    for (it = mapLearners.begin(); it != mapLearners.end(); ++it)
    {
        Int32 learner_id = (it->first);
        learner = (it->second);

        if (disp_details)
        {
            std::cout << "#####################################################"
                         "#####################\n"
                      << "LEARNER" << learner_id << std::endl;
        }
        Log.LogDetail("  ZClassifier: LEARNER #%d", learner_id);

        double learner_sc;
        Float64 learner_score;

        auto mu = learner->m_SVmu;
        auto sigma = learner->m_SVsigma;

        GetXc(mu, sigma);

        if (classifierStore.GetOptionClassifier() == 1)
        {
            // OPTION 1:                Compute: (x_c * sc')*sv_L *alpha  + bias

            M = learner->m_SVectors->size1; // nb_sVectors
            P = learner->m_SVectors->size2; // nb_descriptors
            if (disp_details)
            {
                std::cout << "nbCLASSIF = " << M << "\t"
                          << "nbDESCRIPTORS = " << P << std::endl;
            }

            gsl_vector *prod_xsvt;

            gsl_matrix *sv_transpose =
                gsl_matrix_alloc(P, M); // Size (learner->m_SVectors) = [M x P]
            gsl_matrix_transpose_memcpy(sv_transpose, learner->m_SVectors);

            gsl_vector *prod_svAlpha = gsl_vector_alloc(M);
            gsl_vector_memcpy(prod_svAlpha,
                              learner->m_SValpha); // Copy alpha into temp

            prod_xsvt = GetProductKL(
                m_Xc, sv_transpose); // Vector c [1 x M] = PRODUCT( Vector a [1
                                     // x P] , Matrix b [size P x M]). Element
                                     // c[j] = SUM_i( a[i] * b[i,j] )
            gsl_vector_mul(
                prod_svAlpha,
                learner->m_SVectorLabels); // Vector c[1:K , j] = ( Vector a
                                           // [1xK] .* Vector b [1xK] ). Element
                                           // c[i,j] = a[i] * b[i]

            // gsl_blas_ddot( prod_xsvt, prod_svAlpha, &learner_sc);         //
            // Float c =  ( Vector a [1xP] * Vector b [1xP] ). Element c =
            // SUM_j(a[j] * b[j] )
            learner_sc = 0;
            for (Int32 lig = 0; lig < prod_xsvt->size; lig++)
            {
                learner_sc += (gsl_vector_get(prod_xsvt, lig)) *
                              (gsl_vector_get(prod_svAlpha, lig));
            }
            learner_sc += learner->m_SVbias;

            // SCORE TRANSFORM
            GetSigmoid(learner_sc, learner->m_SVsigmoiid, learner_score);

            // FREE ALLOC_MEMORY
            gsl_vector_free(prod_xsvt);
            gsl_matrix_free(sv_transpose);
            gsl_vector_free(prod_svAlpha);

        } else if (classifierStore.GetOptionClassifier() == 2)
        {
            // OPTION 2:            Compute: (x_c * beta) + bias

            double learner_sc2 = 0;
            for (Int32 lig = 0; lig < learner->m_SVbeta->size; lig++)
            {
                learner_sc2 += (gsl_vector_get(m_Xc, lig)) *
                               (gsl_vector_get(learner->m_SVbeta, lig));
            }
            // gsl_vector* prod_svBeta = gsl_vector_alloc(M);
            // gsl_vector_memcpy( prod_svBeta, learner->m_SVbeta);
            // gsl_blas_ddot( m_Xc, prod_svBeta, &learner_sc2);
            learner_sc2 += learner->m_SVbias;
            GetSigmoid(learner_sc2, learner->m_SVsigmoiid, learner_score);
            learner_sc = learner_sc2;
        } else
        {
            Log.LogError("Unable to find a correct classifier option: found %d",
                         classifierStore.GetOptionClassifier());
        }

        if (disp_details)
        {
            std::cout << "SCORE LEARNER _= " << learner_id << "    "
                      << learner_sc << std::endl;
        }
        // UPDATE SCORE FOR EACH LEARNER
        m_score->data[learner_id - 1] = learner_score;
    }
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Get the posterior class probabilities (combine coding
 * strategy ECOC and individual responses of learners)
 *  ------------------------------------------------------------------------------------------------------
 */
void CQualz::GetPosteriorPred(CClassifierStore &classifierStore)
{
    Int32 K = classifierStore.GetNbClasses();
    Int32 L = classifierStore.GetNbLearners();

    Float64 value = 1 / ((Float64)K);

    gsl_vector *r = gsl_vector_alloc(K);
    gsl_vector_memcpy(r, m_score);

    gsl_vector *w = gsl_vector_alloc(L);
    gsl_vector_memcpy(w, classifierStore.GetLearnersWeight());

    Float64 dist_iter1 = 0, dist_iter2 = 0;
    gsl_vector *p_null = gsl_vector_alloc(1);
    p_null->data[0] = NAN;
    gsl_vector *p_init = gsl_vector_alloc(K);
    for (Int32 i = 0; i < K; i++)
    {
        p_init->data[i] = value;
    }

    // COMPUTE KL MIN_DISTANCE
    gsl_vector *p_iter1 =
        GetArgminKL(classifierStore, r, w, p_init, dist_iter1);
    gsl_vector *p_iter2 =
        GetArgminKL(classifierStore, r, w, p_null, dist_iter2);

    // UPDATE POSTERIOR CLASS PROBABILITIES
    if (dist_iter2 < dist_iter1)
    {
        gsl_vector_memcpy(m_posterior, p_iter2);
    } else
    {
        gsl_vector_memcpy(m_posterior, p_iter1);
    }

    // RESCALE IF NEEDED
    Float64 sump = 0;
    Float64 threshold = 100 * EPS_TOL;
    for (Int32 i = 0; i < m_posterior->size; i++)
    {
        if (m_posterior->data[i] < threshold)
        {
            m_posterior->data[i] = threshold;
        }
        sump += m_posterior->data[i]; // sump = gsl_blas_dasum( p ); // sum of
                                      // absolute elements . the probability
                                      // vector p is postive anyway
    }
    gsl_vector_scale(m_posterior, ((Float64)1 / sump));

    // FREE ALLOC_MEMORY
    gsl_vector_free(r);
    gsl_vector_free(w);
    gsl_vector_free(p_init);
    gsl_vector_free(p_null);
    gsl_vector_free(p_iter1);
    gsl_vector_free(p_iter2);
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute Kullback-Leibler minimization to compute
 * posterior class probabilities (here, numfits =2)
 *  ------------------------------------------------------------------------------------------------------
 */
gsl_vector *CQualz::GetArgminKL(CClassifierStore &classifierStore,
                                gsl_vector *r, gsl_vector *w_learner,
                                gsl_vector *p0, Float64 &distance)
{
    double sump;
    Bool opt_row = true;
    Float64 threshold = 100 * EPS_TOL;

    Int32 K = classifierStore.GetCodingMatrix()->size1; // nbClasses
    Int32 L = classifierStore.GetCodingMatrix()->size2; // nbLearners

    gsl_vector *p;

    gsl_vector *tempPos;
    gsl_vector *tempNeg;
    gsl_vector *tempPosNeg;

    /*
    gsl_matrix* Mpos = gsl_matrix_alloc(K,L);
    gsl_matrix* Mneg = gsl_matrix_alloc(K,L);
    gsl_matrix_memcpy( Mpos, classifierStore.GetCodingMatrixPos());
    gsl_matrix_memcpy( Mneg, classifierStore.GetCodingMatrixNeg());
     */

    p = gsl_vector_alloc(K);
    gsl_vector *r_estim = gsl_vector_alloc(r->size);

    // UPDATE P
    if (std::isnan(p0->data[0]))
    {
        gsl_matrix *M = gsl_matrix_alloc(K, L);
        gsl_matrix_memcpy(M, classifierStore.GetCodingMatrixPos());
        gsl_matrix_add(
            M, classifierStore.GetCodingMatrixNeg()); // M = Mpos + Mneg .
        for (Int32 i = 0; i < M->size1; i++)
        {
            for (Int32 j = 0; j < M->size2; j++)
            {
                if (gsl_matrix_get(M, i, j) == -1)
                {
                    gsl_matrix_set(M, i, j, 0);
                }
            }
        }
        gsl_matrix_transpose(M);

        Bool doquit = GetLSQnonNegKL(M, r, p);
        gsl_matrix_free(M);

        if (doquit)
        {
            // UPDATE R_ESTIM
            tempPos = GetSumKL(
                GetTimesKL(classifierStore.GetCodingMatrixPos(), p), opt_row);
            tempNeg = GetSumKL(
                GetTimesKL(classifierStore.GetCodingMatrixNeg(), p), opt_row);

            tempPosNeg = gsl_vector_alloc(tempPos->size);
            gsl_vector_memcpy(tempPosNeg, tempPos);
            gsl_vector_sub(tempPosNeg, tempNeg);
            gsl_vector_div(
                tempPos, tempPosNeg); // r_estim = tempPos / (tempPos - tempNeg)
            gsl_vector_memcpy(r_estim, tempPos); // Store "r_estim"

            // UPDATE DISTANCE
            distance = GetDistanceKL(r, r_estim, w_learner); // Store "distance"

            // FREE ALLOC_MEMORY
            gsl_vector_free(tempPos);
            gsl_vector_free(tempNeg);
            gsl_vector_free(tempPosNeg);

            return p;
        }

        // UPDATE THE Vector P
        sump = 0;
        for (Int32 i = 0; i < p->size; i++)
        {
            if (p->data[i] < threshold)
            {
                p->data[i] = threshold;
            }
            sump += p->data[i]; // sump = gsl_blas_dasum( p ); // sum of
                                // absolute elements . the probability vector p
                                // is postive anyway
        }

        gsl_vector_scale(p, ((Float64)1 / sump));
        for (Int32 i = 0; i < p->size; i++)
        {
            if (p->data[i] > 1.0)
            {
                p->data[i] = 1.0;
            }
        }
    } else
    {
        gsl_vector_memcpy(p, p0);
    }

    // UPDATE R_ESTIM
    tempPos =
        GetSumKL(GetTimesKL(classifierStore.GetCodingMatrixPos(), p), opt_row);
    tempNeg =
        GetSumKL(GetTimesKL(classifierStore.GetCodingMatrixNeg(), p), opt_row);
    tempPosNeg = gsl_vector_alloc(tempPos->size);
    gsl_vector_memcpy(tempPosNeg, tempPos);
    gsl_vector_sub(tempPosNeg, tempNeg);
    gsl_vector_div(tempPos,
                   tempPosNeg); // r_estim = tempPos ./ (tempPos - tempNeg)
    gsl_vector_memcpy(r_estim, tempPos);

    // UPDATE DISTANCE
    distance = GetDistanceKL(r, r_estim, w_learner);

    Int32 iter = 0, maxIter = 1e3;
    Float64 delta = INFINITY;
    Float64 dist_updated;

    while ((iter <= maxIter) && (delta > 1e-6))
    {
        iter++;

        // UPDATE  P
        p = GetNumDenKL(classifierStore, r, r_estim, p);

        sump = 0;
        for (Int32 i = 0; i < p->size; i++)
        {
            if (p->data[i] < threshold)
            {
                p->data[i] = threshold;
            }
            sump += p->data[i]; // sump = gsl_blas_dasum( p ); // sum of
                                // absolute elements . the probability vector p
                                // is postive anyway
        }
        gsl_vector_scale(p, ((Float64)1 / sump));

        // UPDATE R_ESTIM
        tempPos = GetSumKL(GetTimesKL(classifierStore.GetCodingMatrixPos(), p),
                           opt_row);
        tempNeg = GetSumKL(GetTimesKL(classifierStore.GetCodingMatrixNeg(), p),
                           opt_row);
        tempPosNeg = gsl_vector_alloc(tempPos->size);
        gsl_vector_memcpy(tempPosNeg, tempPos);
        gsl_vector_sub(tempPosNeg, tempNeg);
        gsl_vector_div(tempPos,
                       tempPosNeg); // r_estim = tempPos ./ (tempPos - tempNeg)
        gsl_vector_memcpy(r_estim, tempPos); // Store "r_estim"

        // UPDATE DISTANCE
        dist_updated = GetDistanceKL(r, r_estim, w_learner); // Store "distance"
        delta = distance - dist_updated;
        distance = dist_updated;
    }

    // FINAL ITERATION
    if (0)
    {
        std::cout << " >>>>>>>>>> ITER FINAL  <<<<<<<<<<<< :  " << iter
                  << std::endl;
        std::cout << "\t >>   [" << iter << "] delta = " << delta
                  << "\n"
                     "\t >>   ["
                  << iter << "] distance = " << distance
                  << "\n"
                     "\t >>   ["
                  << iter << "]  proba p = " << p->data[0] << "\t" << p->data[1]
                  << "\t" << p->data[2] << "\t" << p->data[3] << "\t"
                  << p->data[4] << "\n"
                  << std::endl;
    }

    // FREE ALLOC_MEMORY
    gsl_vector_free(tempPos);
    gsl_vector_free(tempNeg);
    gsl_vector_free(tempPosNeg);

    return p;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute c.*p = d             -> ( p= c\d) with
 * positive constraints on p
 *  ------------------------------------------------------------------------------------------------------
 */
Bool CQualz::GetLSQnonNegKL(gsl_matrix *c, gsl_vector *d, gsl_vector *p)
{
    Int32 s, K = c->size2;
    gsl_permutation *permut = gsl_permutation_alloc(
        K); // OVA coding strategy (the coding matrix is square)
    gsl_linalg_LU_decomp(c, permut, &s);
    gsl_linalg_LU_solve(c, permut, d, p);

    // FREE ALLOC_MEMORY
    gsl_permutation_free(permut);

    // UPDATE P
    Bool doquit = false;
    Int32 counter_p0 = 0, counter_p1 = 0, counter_p1id = 0;
    for (Int32 i = 0; i < p->size; i++)
    {
        if (p->data[i] == 0)
        {
            counter_p0++;
        }
        if (p->data[i] > 0)
        {
            counter_p1++;
            counter_p1id = i;
        }
    }
    if (counter_p0 == p->size)
    {
        for (Int32 i = 0; i < p->size; i++)
        {
            p->data[i] = 1 / K;
        }
        doquit = true;
    }
    if (counter_p1 == 1)
    {
        p->data[counter_p1id] = 1;
        doquit = true;
    }
    return doquit;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Intermediate step to minimize KL distance in function
 * "GetArgminKL(..)"
 *  ------------------------------------------------------------------------------------------------------
 */
gsl_vector *CQualz::GetNumDenKL(CClassifierStore &classifierStore,
                                gsl_vector *r, gsl_vector *r_estim,
                                gsl_vector *pold)
{
    Int32 K = classifierStore.GetCodingMatrix()->size1; // nbClasses
    Int32 L = classifierStore.GetCodingMatrix()->size2; // nbLearners

    gsl_vector *w = gsl_vector_alloc(L);
    gsl_vector_memcpy(w, classifierStore.GetLearnersWeight());

    // COMPUTE:
    //                                      SUM_[ @times( Mpos, w.*r ) - @times(
    //                                      Mneg, w.*r ) ]      and      SUM_[
    //                                      @times( Mpos, w.*rEstim ) - @times(
    //                                      Mneg, w.*rEstim ) ]
    gsl_vector *wr = gsl_vector_alloc(L);
    gsl_vector *wr1 = gsl_vector_alloc(L);
    gsl_vector *wr_estim = gsl_vector_alloc(L);
    gsl_vector *wr_estim1 = gsl_vector_alloc(L);

    gsl_vector_memcpy(wr, r);
    gsl_vector_memcpy(wr1, r);
    gsl_vector_memcpy(wr_estim, r_estim);
    gsl_vector_memcpy(wr_estim1, r_estim);

    gsl_vector_scale(wr1, -1);
    gsl_vector_add_constant(wr1, 1); //  = ( 1 - r )
    gsl_vector_scale(wr_estim1, -1);
    gsl_vector_add_constant(wr_estim1, 1); //= ( 1 - rEstim )

    gsl_vector_mul(wr, w);        //  = W.*r
    gsl_vector_mul(wr1, w);       //  = W.*( 1 - r )
    gsl_vector_mul(wr_estim, w);  //  = W.*rEstim
    gsl_vector_mul(wr_estim1, w); //  = W.*( 1 - rEstim )

    gsl_matrix *tempPos = GetTimesKL(classifierStore.GetCodingMatrixPos(),
                                     wr); // = @times( Mpos, w.*r )
    gsl_matrix *tempNeg = GetTimesKL(classifierStore.GetCodingMatrixNeg(),
                                     wr1); // = @times( Mneg, w.*(1-r) )
    gsl_matrix *tempPos_ = GetTimesKL(classifierStore.GetCodingMatrixPos(),
                                      wr_estim); // = @times( Mpos, w.*rEstim )
    gsl_matrix *tempNeg_ =
        GetTimesKL(classifierStore.GetCodingMatrixNeg(),
                   wr_estim1); // = @times( Mneg, w.*(1-rEstim) )

    gsl_matrix_sub(tempPos, tempNeg);   // result stored in tempPos
    gsl_matrix_sub(tempPos_, tempNeg_); // result stored in tempPos_

    // COMPUTE NUM & DEN
    gsl_vector *p = gsl_vector_alloc(K);
    gsl_vector *cond1 = gsl_vector_alloc(K);
    gsl_vector *cond2 = gsl_vector_alloc(K);
    gsl_vector *numer = gsl_vector_alloc(K);
    gsl_vector *denom = gsl_vector_alloc(K);

    Float64 s1, s2;
    Int32 cond_counter = 0;
    for (Int32 idrow = 0; idrow < K; idrow++)
    { // a sum in each row of the matrices
        s1 = 0.0;
        s2 = 0.0;
        for (Int32 idcol = 0; idcol < L; idcol++)
        {
            s1 += gsl_matrix_get(tempPos, idrow, idcol);
            s2 += gsl_matrix_get(tempPos_, idrow, idcol);
        }
        numer->data[idrow] = s1;
        denom->data[idrow] = s2;

        cond1->data[idrow] = (s2 <= 0) & (s1 > 0);
        cond2->data[idrow] = (s2 > 0);
        if (cond1->data[idrow] > 0)
        {
            cond_counter++;
        }
    }

    // UPDATE P
    if (cond_counter > 0)
    { // at least one number is !NAN & !0 in this vector
        for (Int32 idrow = 0; idrow < K; idrow++)
        {
            if (cond1->data[idrow])
            {
                p->data[idrow] = 1;
            } else
            {
                p->data[idrow] = 0;
            }
        }
    } else
    {
        for (Int32 idrow = 0; idrow < K; idrow++)
        {
            Bool conditio = (cond2->data[idrow]);
            if (conditio)
            {
                Float64 value = ((Float64)numer->data[idrow]) /
                                ((Float64)denom->data[idrow]);
                p->data[idrow] = pold->data[idrow] * value;
            } else
            {
                p->data[idrow] = 0;
            }
        }
    }

    // FREE ALLOC_MEMORY
    gsl_vector_free(w);
    gsl_vector_free(wr);
    gsl_vector_free(wr1);
    gsl_vector_free(wr_estim);
    gsl_vector_free(wr_estim1);

    gsl_matrix_free(tempPos);
    gsl_matrix_free(tempNeg);
    gsl_matrix_free(tempPos_);
    gsl_matrix_free(tempNeg_);

    gsl_vector_free(cond1);
    gsl_vector_free(cond2);
    gsl_vector_free(numer);
    gsl_vector_free(denom);

    return p;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute the distance KL
 *  ------------------------------------------------------------------------------------------------------
 */
Float64 CQualz::GetDistanceKL(gsl_vector *r, gsl_vector *r_estim,
                              gsl_vector *w_learner)
{
    Float64 threshold = 100 * EPS_TOL;
    Float64 s1 = 0.0, s2 = 0.0, distance = 0.0;

    for (Int32 i = 0; i < r->size; i++)
    {
        if (r->data[i] > threshold)
        {
            s1 += w_learner->data[i] * r->data[i] *
                  log(r->data[i] / r_estim->data[i]);
        }
        if ((1 - r->data[i]) > threshold)
        {
            s2 += w_learner->data[i] * (1 - r->data[i]) *
                  log((1 - r->data[i]) / (1 - r_estim->data[i]));
        }
    }
    distance = s1 + s2;

    return distance;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute the norm L1 of a gsl_matrix
 *
------------------------------------------------------------------------------------------------------
* Float64 CQualz::GetNormKL ( gsl_matrix* m )
{
        gsl_vector* temp = gsl_vector_alloc( m->size1 );
        Float64 s, maxi = 0.0;

        for ( Int32 j = 0; j < m->size2 ; j++) {
                gsl_matrix_get_col ( temp, m, j );

                //s = gsl_blas_dasum ( temp );          // absolute sum \sum
|x_i| of the elements of the vector x. s=0; for ( Int32 i = 0; i< temp->size;
i++) { s+=temp->data[i]; //sump = gsl_blas_dasum( p ); // sum of absolute
elements . the probability vector p is postive anyway
                }
                maxi = max( maxi, s );
        }
        Float64 normL1 = maxi;
        gsl_vector_free( temp );

        return normL1;
}
*/

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Center the feature vector with SV parameters (mean &
 * dispersion)
 *  ------------------------------------------------------------------------------------------------------
 */
void CQualz::GetXc(gsl_vector *mu, gsl_vector *sigma)
{
    /*
    std::string desc_1 = "Localized_prior";                                 //0
    std::string desc_2 = "%z candidats";                                    //1
    std::string desc_3 = "std_p(z|D,I) "; //2 std::string desc_4 =
    "skewness_p(z|D,I)";                       //3 std::string desc_5 =
    "kurtosis_p(z|D,I)";                               //4 std::string desc_6 =
    "p(zmap_1d|D,I)";                                          //5 std::string
    desc_7 = "#z_in_R1"; //6 std::string desc_8 = "dz_in_R1"; //7 std::string
    desc_9 = "sum_p(z_in_R1|D,I)";                                      //8
    std::string desc_10 = "sum_p(z_in_R2|D,I)";                             //9
    std::string desc_11 = "dist_(2peaks)_in_dz";            //10
    std::string desc_12 = "dist_(2peaks_in_proba";                  //11
    std::string desc_13 = "#km_peaks_1d"; //12
    */
    TInt64List selected = {2, 5, 6, 7, 8, -1, 9, 11, 12};
    Float64 xc;
    m_Xc = gsl_vector_alloc(mu->size);
    for (Int32 i = 0; i < mu->size; i++)
    {
        if (selected[i] == -1)
        {
            xc = 2e-3; // unused descriptor for recent tests (but still in the
                       // old training set). No effect, bcse cste in the current
                       // config
        } else
        {
            xc = m_zfeatures[selected[i]];
        }
        m_Xc->data[i] = (xc - mu->data[i]) / sigma->data[i];
    }
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute sigmoid function
 *  ------------------------------------------------------------------------------------------------------
 */
void CQualz::GetSigmoid(Float64 &sc, TFloat64List &paramsSig, Float64 &result)
{
    result = 1 + exp(paramsSig[0] * sc + paramsSig[1]);
    result = 1 / result;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute the product of a vector Y and matrix SV
 *                                              Output Result [1 x M] = PRODUCT(
 * Y [1 x P] , SV [size P x M] ). Element Result[j] = SUM_i( Y[i] * SV[i,j] )
 *  ------------------------------------------------------------------------------------------------------
 */
gsl_vector *CQualz::GetProductKL(gsl_vector *y, gsl_matrix *sv)
{
    gsl_vector *result = gsl_vector_alloc(y->size);
    // gsl_vector* temp = gsl_vector_alloc( sv->size1 );
    for (Int32 col = 0; col < sv->size2; col++)
    {
        double prod = 0;
        for (Int32 lig = 0; lig < sv->size1; lig++)
        {
            prod += (gsl_vector_get(y, lig)) * (gsl_matrix_get(sv, lig, col));
        }
        // gsl_matrix_get_col ( temp, sv, col );
        // gsl_blas_ddot( y, temp, &prod ); // DOUBLE c =  ( Vector a [1xP] *
        // Vector b [1xP] ). Element c = SUM_j( a[j] * b[j] )
        result->data[col] = prod;
    }
    // gsl_vector_free(temp);

    return result;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute @times function(vector Y , matrix M)
 *                                              Output Result [P x M] = PRODUCT(
 * V [P x 1] , M [size P x M] ). Column Result[:,j] = ( Y[i] * M[i,j] )_i
 *  ------------------------------------------------------------------------------------------------------
 */
gsl_matrix *CQualz::GetTimesKL(const gsl_matrix *m, gsl_vector *v)
{
    gsl_matrix *result = gsl_matrix_alloc(m->size1, m->size2);
    gsl_vector *temp = gsl_vector_alloc(m->size1);
    for (Int32 j = 0; j < m->size2; j++)
    {
        gsl_matrix_get_col(temp, m, j);
        gsl_vector_mul(temp,
                       v); // vector c[1:K , j] = ( Vector a [1xK] .* Vector b
                           // [1xK] ). Element c[i,j] = a[i] * b[i]
        gsl_matrix_set_row(result, j, temp);
    }
    gsl_vector_free(temp);
    return result;
}

/*  ------------------------------------------------------------------------------------------------------
 *      [ Function ]
 *                      >> Compute the sum of a matrix (on each column, or on
 * each row)
 *  ------------------------------------------------------------------------------------------------------
 */
gsl_vector *CQualz::GetSumKL(gsl_matrix *m, Bool opt_row)
{
    Float64 s = 0.0;

    gsl_vector *result;

    if (opt_row)
    {
        result = gsl_vector_alloc(m->size1);
        for (Int32 i = 0; i < m->size1; i++)
        {
            s = 0;
            for (Int32 j = 0; j < m->size2; j++)
            {
                s += gsl_matrix_get(m, i, j);
            }
            result->data[i] = s;
        }
    } else
    {
        result = gsl_vector_alloc(m->size2);
        for (Int32 j = 0; j < m->size2; j++)
        {
            s = 0;
            for (Int32 i = 0; i < m->size1; i++)
            {
                s += gsl_matrix_get(m, i, j);
            }
            result->data[j] = s;
        }
    }

    return result;
}

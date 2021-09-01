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
#ifndef _REDSHIFT_RAY_CATALOGSTPLSHAPE_
#define _REDSHIFT_RAY_CATALOGSTPLSHAPE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/ray.h"
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{
    class CLineModelElementList;

/**
 * \ingroup Redshift
 */
class CRayCatalogsTplShape
{

public:

    Bool Init(std::string calibrationPath, 
              std::string opt_tplratioCatRelPath, 
              Int32 enableISMCalzetti,
              std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti,
              Float64 nsigmasupport);

    Bool Load( const char* dirPath );
    bool LoadVelocities( const char* filepath, Int32 k );
    bool LoadPrior( const char* filepath, Int32 k );

    //Bool AreCatalogsAligned( const CRayCatalog::TRayVector& restRayList, Int32 typeFilter, Int32 forceFilter  );
    Float64 GetBestFit(const CRayCatalog::TRayVector& restRayList, std::vector<Float64> fittedAmplitudes, std::vector<Float64> fittedErrors, std::vector<Float64> &amplitudesCorrected , std::string &bestTplName);
    CRayCatalog::TRayVector GetRestLinesList( const Int32 index );
    Int32 GetCatalogsCount();
    const std::vector<Float64> & getCatalogsPriors();
    std::string GetCatalogName(Int32 idx);
    Int32 GetIsmIndex(Int32 idx);
    Float64 GetIsmCoeff(Int32 idx);

    Bool GetCatalogVelocities(Int32 idx, Float64& elv, Float64& alv );
    Bool SetMultilineNominalAmplitudes(CLineModelElementList& LineModelElementList, Int32 iLine);
    Bool SetLyaProfile(CLineModelElementList &LineModelElementList, Int32 iCatalog, bool forceLyaFitting, const Float64 nsigmasupport);
    Bool InitLineCorrespondingAmplitudes(CLineModelElementList &LineModelElementList);
    Bool SetMultilineNominalAmplitudesFast(CLineModelElementList &LineModelElementList, Int32 iCatalog);

private:
    Float64 GetFit(std::vector<Float64> ampsLM, std::vector<Float64> errLM, std::vector<Float64> ampsTPL , std::vector<Float64> &ampsCorrected);

    std::string tplshapedcatalog_relpath;

    std::vector<std::string> m_RayCatalogNames;
    std::vector<CRayCatalog> m_RayCatalogList;
    std::vector<std::vector<std::vector<Float64>>> m_RayCatalogLinesCorrespondingNominalAmp;
    std::vector<Float64> m_ELvelocities;
    std::vector<Float64> m_ABSvelocities;
    std::vector<Float64> m_Priors;
    std::vector<Int32> m_IsmIndexes;

    std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
    Int32 m_opt_dust_calzetti;
    Float64 m_nsigmasupport;
};


}

#endif

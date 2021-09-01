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
#ifndef _REDSHIFT_STATISTICS_PRIORHELPERCONTINUUM_
#define _REDSHIFT_STATISTICS_PRIORHELPERCONTINUUM_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <boost/format.hpp>
#include <cfloat>
#include <vector>
#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CPriorHelperContinuum
{

public:

    struct SPriorTZE{
        Float64 logpriorTZE;
        Float64 A_sigma;
        Float64 A_mean;
    };
    typedef std::vector<SPriorTZE> TPriorEList;
    typedef std::vector<TPriorEList> TPriorZEList;
    typedef std::vector<TPriorZEList> TPriorTZEList;

    CPriorHelperContinuum();
    ~CPriorHelperContinuum();

    bool Init( std::string priorDirPath );
    bool LoadFileEZ(const char* filePath , std::vector<std::vector<Float64>>& data);

    bool SetSize(UInt32 size);
    bool SetTNameData(UInt32 k, std::string tname);
    bool SetEZTData(UInt32 k, std::vector<std::vector<Float64>> ezt_data);
    bool SetAGaussmeanData(UInt32 k, std::vector<std::vector<Float64>> agaussmean_data);
    bool SetAGausssigmaData(UInt32 k, std::vector<std::vector<Float64>> agausssigma_data);

    bool GetTplPriorData(std::string tplname,
                         std::vector<Float64> redshifts,
                         TPriorZEList &zePriorData,
                         Int32 outsideZRangeExtensionMode=0);

    bool SetBeta(Float64 beta);

    bool mInitFailed = false;

private:

    TPriorTZEList m_data;
    std::vector<std::string> m_tplnames;

    UInt32 m_nZ = 24;
    Float64 m_dz = 0.25;
    Float64 m_z0 = 0.0;

    UInt32 m_nEbv = 10;
    Float64 m_ebv0 = 0.0;

    Float64 m_beta = -1;

    Float64 m_priorminval = DBL_MIN;
};

}

#endif

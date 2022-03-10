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
#ifndef _REDSHIFT_COMMON_FLAG_
#define _REDSHIFT_COMMON_FLAG_

#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/common/mutex.h"

#include <vector>
#include <cstdarg>

#define Flag (CFlagWarning::GetInstance())

namespace NSEpic
{

class CFlagWarning : public CSingleton<CFlagWarning>
{
public:
    typedef enum WarningCode
        {
        WARNING_NONE=0,
        AIR_VACCUM_CONVERSION_IGNORED, //1
        CRANGE_VALUE_OUTSIDERANGE, //2
        CRANGE_VECTBORDERS_OUTSIDERANGE, //3
        CRANGE_NO_INTERSECTION, //4
        FINDER_NO_PEAKS, //5
        STDESTIMATION_NO_MATCHING, //6
        STDESTIMATION_FAILED, //7
        MULTIROLL_STRTAG_NOTFOUND, //8
        LINEMATCHING_REACHED_ENDLOOP, //9
        FORCE_LOGSAMPLING_FFT,//10
        IGNORELINESSUPPORT_DISABLED_FFT,//11
        FORCE_FROMSPECTRUM_NEG_CONTINUUMAMP,//12
        INVALID_MERIT_VALUES,//13
        AIR_VACCUM_REACHED_MAX_ITERATIONS, //14
        ASYMFIT_NAN_PARAMS,//15
        DELTAZ_COMPUTATION_FAILED, //16
        INVALID_FOLDER_PATH,//17
        TPL_NAME_EMPTY,//18
	RELIABILITY_NEEDS_TENSORFLOW//19
        } WarningCode;

    void warning(WarningCode c, std::string message);
    void warning(WarningCode c, const char* format, ... );
    
    const std::vector<std::pair<WarningCode, std::string>>& getListMessages();
    Int32 getBitMask();
    void resetFlag();

private:
    friend class CSingleton<CFlagWarning>;

    CFlagWarning();
    ~CFlagWarning();

    UInt32          m_flags;
    Char*           m_WorkingBuffer;
    std::vector<std::pair<WarningCode, std::string>> m_messageList;
};

typedef std::vector<std::pair<CFlagWarning::WarningCode, std::string>> TWarningMsgList;

}
#endif

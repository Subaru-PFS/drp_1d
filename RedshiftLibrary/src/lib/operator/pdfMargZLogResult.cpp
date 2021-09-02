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
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>
#include "RedshiftLibrary/log/log.h"
#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace NSEpic;


CPdfMargZLogResult::CPdfMargZLogResult()
{
  this->m_type = "CPdfMargZLogResult";

}


CPdfMargZLogResult::CPdfMargZLogResult(const TFloat64List & redshifts):
    Redshifts(redshifts),
    countTPL(redshifts.size()), // assumed 1 model per z
    valProbaLog(redshifts.size(), -DBL_MAX),
    valEvidenceLog(-1.0)
{
  this->m_type = "CPdfMargZLogResult";
}

Int32 CPdfMargZLogResult::getIndex( Float64 z ) const
{
    Int32 solutionIdx=-1;
    for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
    {
        if( Redshifts[i2]==z )
        {
            solutionIdx = i2;
            break;
        }
    }
    return solutionIdx;
}



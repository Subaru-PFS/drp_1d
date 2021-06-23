#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <RedshiftLibrary/log/log.h>
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



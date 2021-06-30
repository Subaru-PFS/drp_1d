#include "RedshiftLibrary/operator/pdfLogresult.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

using namespace NSEpic;


CPdfLogResult::CPdfLogResult()
{

}

CPdfLogResult::~CPdfLogResult()
{

}

void CPdfLogResult::SetSize( UInt32 n )
{
    Redshifts.resize(n);
    valProbaLog.resize(n);
    Overlap.resize(n);
    Status.resize(n);
}


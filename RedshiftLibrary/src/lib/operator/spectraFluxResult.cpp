#include <RedshiftLibrary/operator/spectraFluxResult.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace NSEpic;


CSpectraFluxResult::CSpectraFluxResult()
{
  this->m_type = "CSpectraFluxResult";
}

CSpectraFluxResult::~CSpectraFluxResult()
{

}


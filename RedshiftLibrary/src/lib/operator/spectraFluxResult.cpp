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

}

CSpectraFluxResult::~CSpectraFluxResult()
{
	//m_optio = 1;
}

void CSpectraFluxResult::getData(const std::string& name, double **data, int *size) const
{
  *size = wavel.size();
  if(name.compare("BestContinuumLambda") == 0)
    {
      *data = const_cast<double *>(wavel.data());
    }
 else if(name.compare("BestContinuumFlux") == 0)
    {
      *data = const_cast<double *>(fluxes.data());
    }
 else Log.LogError("Unkwnown data %s",name.c_str());

}

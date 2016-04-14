#include <epic/redshift/operator/spectraFluxResult.h>

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

Void CSpectraFluxResult::Save( const CDataStore& store, std::ostream& stream ) const
{
	UInt32 p = 13;
	std::string tx;
	if        (m_optio ==0) { tx="# Wavelength \tFlux"; }
	else if (m_optio ==1) { tx="# Wavelength \tFilteredFlux";   }
	else if (m_optio ==2) { tx="# Wavelength \tErrorFiltered";  }

	stream << tx<< std::endl;
	for ( Int32 i=0; i<fluxes.size(); i++)
	{
		stream.precision(10);
		stream  <<  wavel[i] << "\t" ;

		stream.precision(p);
		stream<< fluxes[i] << std::endl;
	}
}

Void CSpectraFluxResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
	stream.precision(10);
    stream << "CSpectraFluxResult" << "\t" << fluxes.size() << std::endl;

}

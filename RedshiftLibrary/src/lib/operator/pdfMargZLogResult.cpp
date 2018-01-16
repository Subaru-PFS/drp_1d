#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace NSEpic;


CPdfMargZLogResult::CPdfMargZLogResult()
{
    valEvidenceLog=-1;
}

CPdfMargZLogResult::~CPdfMargZLogResult()
{

}

Void CPdfMargZLogResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream << "#z_tested \tLog P( z | priors )" <<std::endl;
    stream << "#EvidenceLog=" << valEvidenceLog << std::endl;
	for ( Int32 i=0; i<Redshifts.size(); i++)
    {
    		stream.precision(10);
        stream <<  Redshifts[i] << "\t" << valProbaLog[i] << std::endl;
    }
}

Void CPdfMargZLogResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
	stream.precision(20);
    stream << "CPdfMargZLogResult" << "\t" << Redshifts.size() << std::endl;
}

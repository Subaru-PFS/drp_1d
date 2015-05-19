#include <epic/redshift/operator/correlationresult.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CCorrelationResult )

CCorrelationResult::CCorrelationResult()
{

}

CCorrelationResult::~CCorrelationResult()
{

}

Void CCorrelationResult::Save( std::ostream& stream ) const
{
    stream <<  "#Redshifts\tCorrelation\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << "\t" << std::scientific << Correlation[i] << std::fixed << "\t" << Overlap[i] << std::endl;
    }
}

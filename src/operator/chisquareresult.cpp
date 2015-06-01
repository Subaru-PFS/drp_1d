#include <epic/redshift/operator/chisquareresult.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CChisquareResult )

CChisquareResult::CChisquareResult()
{

}

CChisquareResult::~CChisquareResult()
{

}

Void CChisquareResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << "\t" << std::scientific << ChiSquare[i] << std::fixed << "\t" << Overlap[i] << std::endl;
    }
}

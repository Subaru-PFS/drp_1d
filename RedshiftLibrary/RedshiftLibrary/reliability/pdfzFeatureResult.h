#ifndef _REDSHIFT_RELIABILITY_PDFZFEATURERESULT_
#define _REDSHIFT_RELIABILITY_PDFZFEATURERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <boost/unordered_map.hpp>

using namespace std;
namespace NSEpic
{

class CPdfzFeatureResult : public COperatorResult
{

public:

	CPdfzFeatureResult();
    virtual ~CPdfzFeatureResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

	typedef boost::unordered_map<const std::string, Float64> Mapz;
	Mapz mapzfeatures;

};


}

#endif

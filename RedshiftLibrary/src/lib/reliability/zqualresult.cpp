#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/datastore.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>

#include "RedshiftLibrary/reliability/zqualresult.h"
#include "RedshiftLibrary/reliability/pdfzPredictResult.h"

using namespace std;
using namespace NSEpic;

CQualzResult::CQualzResult()
{

}

CQualzResult::~CQualzResult()
{

}

void CQualzResult::Save (std::ostream& stream ) const
{
  /*
	std::string scope = store.GetScope( *this ) + "zReliability/result.zpredict";
	auto zqualResults = store.GetGlobalResult(scope.c_str());

	auto zQual_vect =  std::dynamic_pointer_cast<const CPdfzPredictResult>( zqualResults.lock() );

    if(!zQual_vect)
    {
        Log.LogError( "CQualzResult: no results retrieved from scope: %s", scope.c_str());
        return;
    }


	stream << "#Predicted_Label(zReliability)" << "\t" << zQual_vect->m_predLabel << std::endl;
	stream << "#Class \t PosteriorProbability_Prediction" << "\t"<< std::endl;
	for ( Int32 i = 0; i< zQual_vect->m_Labels.size(); i++)
	{
		stream
		<< zQual_vect->m_Labels[i] << "\t"
		<< std::setprecision(10)
		<<  zQual_vect->m_posterior->data[i]<< "\t"
		<< std::endl;
	}
	stream << "#Learner \t Score_Prediction " << "\t"<< std::endl;
	for ( Int32 i = 0; i< zQual_vect->m_score->size; i++)
	{
		stream
		<< "\n#Learner_"<< (i+1) << "\t"
		<< std::setprecision(10)
		<<  zQual_vect->m_score->data[i] << "\t"
		<< std::endl;
	}
  */
}

void CQualzResult::SaveLine( std::ostream& stream ) const
{
  /*
	std::string scope = store.GetScope( *this ) + "zReliability/result.zpredict";
	auto zqualResults = store.GetGlobalResult(scope.c_str());

	auto zQual_vect =  std::dynamic_pointer_cast<const CPdfzPredictResult>( zqualResults.lock() );

    if(!zQual_vect)
    {
        Log.LogError( "CQualzResult: no results retrieved from scope: %s", scope.c_str());
        return;
    }

	stream
	<< "#zReliability_Pred \t" <<zQual_vect->m_predLabel<< "\t"
	<< "PosteriorProba_Pred \t" << std::setprecision(10)<<zQual_vect->m_predProba << "\t"
	<<store.GetSpectrumName()
	<<std::endl;
  */

}

Bool CQualzResult::GetPredictedLabel( const CDataStore& store, std::string& predLabel ) const
{
    predLabel = "";

    std::string scope = store.GetScope( *this ) + "zReliability/result.zpredict";
    auto zqualResults = store.GetGlobalResult(scope.c_str());

    auto zQual_vect =  std::dynamic_pointer_cast<const CPdfzPredictResult>( zqualResults.lock() );

    if(!zQual_vect)
    {
        Log.LogError( "CQualzResult: no results retrieved from scope: %s", scope.c_str());
        return false;
    }

    predLabel = zQual_vect->m_predLabel;

    return true;
}

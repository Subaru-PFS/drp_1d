// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
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

#include <RedshiftLibrary/reliability/pdfzPredictResult.h>

#include <RedshiftLibrary/processflow/datastore.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace NSEpic;


CPdfzPredictResult::CPdfzPredictResult()
{

}

CPdfzPredictResult::~CPdfzPredictResult()
{

}


void CPdfzPredictResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream << "#Predicted_Label(zReliability) "<<" \t" << m_predLabel.c_str() << std::endl;
	stream << "#Class \t PosteriorProbability_Prediction" << "\t"<< std::endl;
    for ( Int32 i = 0; i< m_Labels.size(); i++)
    {
        stream
        << m_Labels[i] << "\t"
        << std::setprecision(10)
        << m_posterior->data[i]<< "\t"
        << std::endl;
    }
	stream << "\n#Learner \t Score_Prediction " << "\t"<< std::endl;
	for ( Int32 i = 0; i< m_score->size; i++)
	{
		stream
		<< "Learner_"<< (i+1) << "\t"
		<< std::setprecision(10)
        << m_score->data[i] << "\t"
		<< std::endl;
	}
}


void CPdfzPredictResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
	stream
    << "#zReliability_Pred \t" <<m_predLabel.c_str() << "\t"
    << "PosteriorProba_Pred \t" << std::setprecision(10)<< m_predProba << "\t"
    << store.GetSpectrumName()
    << std::endl;

}

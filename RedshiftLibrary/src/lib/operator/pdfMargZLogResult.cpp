#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <RedshiftLibrary/log/log.h>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace NSEpic;


CPdfMargZLogResult::CPdfMargZLogResult()
{
    valEvidenceLog=-1;
}

CPdfMargZLogResult::~CPdfMargZLogResult()
{

}

void CPdfMargZLogResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream << "#z_tested \tLog P( z | priors )" <<std::endl;
    stream << "#EvidenceLog=" << valEvidenceLog << std::endl;
    for ( Int32 i=0; i<Redshifts.size(); i++)
    {
            stream.precision(10);
        stream <<  Redshifts[i] << "\t" << valProbaLog[i] << std::endl;
    }
}


Int32 CPdfMargZLogResult::Load( std::string filePath )
{
    // Clear current list
    Redshifts.clear();
    valProbaLog.clear();

    std::ifstream file;
    file.open( filePath, std::ifstream::in );
    if( file.rdstate() & ios_base::failbit )
    {
        return -3;
    }

    std::string line;

    // Read file line by line
    while( std::getline( file, line ) )
    {
        if( boost::starts_with( line, "#" ) )
        {
            continue;
        }
        boost::char_separator<char> sep(" \t");

        // Tokenize each line
        typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() )
        {
            // Parse
            Float64 r = 0.0;
            try
            {
                r = boost::lexical_cast<Float64>(*it);
            }
            catch (boost::bad_lexical_cast&)
            {
                return -1;
            }

            // Parse
            ++it;
            Float64 c = 0.0;
            if( it != tok.end() )
            {
                try
                {
                    c = boost::lexical_cast<Float64>(*it);
                }
                catch (boost::bad_lexical_cast&)
                {
                    return -1;
                }
            }
            else
            {
                return -2;
            }

            Redshifts.push_back( r );
            valProbaLog.push_back( c );
        }
    }
    return 0;
}

void CPdfMargZLogResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
	stream.precision(20);
    stream << "CPdfMargZLogResult" << "\t" << Redshifts.size() << std::endl;
}

Int32 CPdfMargZLogResult::getIndex( Float64 z ) const
{
    Int32 solutionIdx=-1;
    for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
    {
        if( Redshifts[i2]==z )
        {
            solutionIdx = i2;
            break;
        }
    }
    return solutionIdx;
}


  void CPdfMargZLogResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
  {
  }
  void CPdfMargZLogResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
  {
  }
  void CPdfMargZLogResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const{}
  void CPdfMargZLogResult::getCandidateData(const int& rank,const std::string& name, double **data, int *size) const{}

  void CPdfMargZLogResult::getData(const std::string& name, Int32& v) const{}
  void CPdfMargZLogResult::getData(const std::string& name, Float64& v) const{}
  void CPdfMargZLogResult::getData(const std::string& name, std::string& v) const{}
  void CPdfMargZLogResult::getData(const std::string& name, double **data, int *size) const
  {

    if(name.compare("pdf_zgrid") == 0)
      {
        *size = Redshifts.size();
        *data = const_cast<double *>(Redshifts.data());
      }
    if(name.compare("pdf_probaLog") == 0)
      {
        *size = valProbaLog.size();
        *data = const_cast<double *>(valProbaLog.data());
      }
  }

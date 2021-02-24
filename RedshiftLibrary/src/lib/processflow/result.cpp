#include <RedshiftLibrary/processflow/result.h>

using namespace NSEpic;

void COperatorResult::SaveJSON(std::ostream& stream) const
{
  // does nothing, -> no need to cast COperatorResult to LineModelResult in COperatorResultStore::SaveAllResults
}


void COperatorResult::SaveFloat64(std::ostream& stream,Float64 data) const
{
    if (data != data ) stream << "null";
    else stream << data;
}

void COperatorResult::SaveTFloat64List(std::ostream& stream,std::string name,TFloat64List data) const 
{
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
        SaveFloat64(stream,data[i]);
        if( i< data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveTFloat64ListOfList(std::ostream& stream,std::string name,std::vector<TFloat64List> data) const
{
  if(data.size()>0){
    stream <<  "\""<<name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      stream << "[";
      for ( int ip=0; ip<data[i].size(); ip++)
      {
        SaveFloat64(stream,data[i][ip]);
        if ( ip< data[i].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<data.size() - 1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveInt32Vector(std::ostream& stream,std::string name,std::vector<Int32> data) const
{
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      stream <<  data[i];
      if( i<data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveStringVector(std::ostream& stream,std::string name,std::vector<std::string> data) const 
{
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      stream << "\"" << data[i] << "\""; 
      if( i< data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveStringVectorOfVector(std::ostream& stream,std::string name,std::vector<std::vector<std::string>> data) const
{
if(data.size()>0){
    stream <<  "\""<< name <<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      std::vector<std::string> line_list;
      line_list = data[i];
      stream << "[";
      for ( int ki=0; ki<line_list.size(); ki++)
      {
        stream <<  "\"" << line_list[ki] << "\"";
        if (ki < line_list.size()-1) stream <<",";
      }
      stream << "]";
      if (i<data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveTContinuumIndexListVector(std::ostream& stream,std::string name,std::vector<CContinuumIndexes::TContinuumIndexList> data) const
{
  if(data.size()>0){
    stream <<  "\""<<name << "Color\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<data[i].size(); kci++)
      {
	      SaveFloat64(stream,data[i][kci].Color);
        if ( kci<data[i].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<data.size() - 1) stream << ",";
    }

    stream << "],"<<std::endl ;
    stream <<  "\"" << name << "Break\" :[";
    for ( int i=0; i<data.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<data[i].size(); kci++)
      {
	      SaveFloat64(stream,data[i][kci].Break);	
        if ( kci<data[i].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<data.size() - 1) stream << ",";
    }
    stream << "]" ;
  }
}


  void COperatorResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
  {
  }
  void COperatorResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
  {
  }
  void COperatorResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const{}
  void COperatorResult::getCandidateData(const int& rank,const std::string& name, double **data, int *size) const{}
void COperatorResult::getCandidateData(const int& rank,const std::string& name, std::string *data, int *size) const{}
void COperatorResult::getCandidateData(const int& rank,const std::string& name, int **data, int *size) const{}

  void COperatorResult::getData(const std::string& name, Int32& v) const{}
  void COperatorResult::getData(const std::string& name, Float64& v) const{}
  void COperatorResult::getData(const std::string& name, std::string& v) const{}
  void COperatorResult::getData(const std::string& name, double **data, int *size) const{}
  void COperatorResult::getData(const std::string& name, int **data, int *size) const{}
void COperatorResult::getData(const std::string& name, std::string *data, int *size) const{}

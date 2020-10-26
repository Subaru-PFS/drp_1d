#include <RedshiftLibrary/processflow/result.h>

using namespace NSEpic;

COperatorResult::COperatorResult()
{

}

COperatorResult::~COperatorResult()
{

}

void COperatorResult::SaveJSON(const CDataStore& store, std::ostream& stream) const
{
  // does nothing, -> no need to cast COperatorResult to LineModelResult in COperatorResultStore::SaveAllResults
}


void COperatorResult::SaveFloat64(std::ostream& stream,Float64 data) const
{
    if (data != data ) stream << "null";
    else stream << data;
}

void COperatorResult::SaveTFloat64List(std::ostream& stream,std::string name,TFloat64List data, TInt32List order) const 
{
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<order.size(); i++)
    {
        SaveFloat64(stream,data[order[i]]);
        if( i< order.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveTFloat64ListOfList(std::ostream& stream,std::string name,std::vector<TFloat64List> data, TInt32List order) const
{
  if(data.size()>0){
    stream <<  "\""<<name<<"\" : [";
    for ( int i=0; i<order.size(); i++)
    {
      stream << "[";
      for ( int ip=0; ip<data[order[i]].size(); ip++)
      {
        SaveFloat64(stream,data[order[i]][ip]);
        if ( ip< data[order[i]].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<order.size() - 1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveInt32Vector(std::ostream& stream,std::string name,std::vector<Int32> data, TInt32List order) const
{
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<order.size(); i++)
    {
      stream <<  data[order[i]];
      if( i<order.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveStringVector(std::ostream& stream,std::string name,std::vector<std::string> data, TInt32List order) const 
{
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<order.size(); i++)
    {
      stream << "\"" << data[order[i]] << "\""; 
      if( i< order.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveStringVectorOfVector(std::ostream& stream,std::string name,std::vector<std::vector<std::string>> data, TInt32List order) const
{
if(data.size()>0){
    stream <<  "\""<< name <<"\" : [";
    for ( int i=0; i<order.size(); i++)
    {
      std::vector<std::string> line_list;
      line_list = data[order[i]];
      stream << "[";
      for ( int ki=0; ki<line_list.size(); ki++)
      {
        stream <<  "\"" << line_list[ki] << "\"";
        if (ki < line_list.size()-1) stream <<",";
      }
      stream << "]";
      if (i<order.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveTContinuumIndexListVector(std::ostream& stream,std::string name,std::vector<CContinuumIndexes::TContinuumIndexList> data, TInt32List order) const
{
  if(data.size()>0){
    stream <<  "\""<<name << "Color\" : [";
    for ( int i=0; i<order.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<data[order[i]].size(); kci++)
      {
	      SaveFloat64(stream,data[order[i]][kci].Color);
        if ( kci<data[order[i]].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<order.size() - 1) stream << ",";
    }

    stream << "],"<<std::endl ;
    stream <<  "\"" << name << "Break\" :[";
    for ( int i=0; i<order.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<data[order[i]].size(); kci++)
      {
	      SaveFloat64(stream,data[order[i]][kci].Break);	
        if ( kci<data[order[i]].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<order.size() - 1) stream << ",";
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

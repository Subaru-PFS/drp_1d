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

void COperatorResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}

void COperatorResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}

void COperatorResult::SaveFloat64(std::ostream& stream,Float64 data) const
{
    if (data != data ) stream << "null";
    else stream << data;
}

void COperatorResult::SaveTFloat64List(std::ostream& stream,std::string name,TFloat64List data, TFloat64List order) const 
{
  Bool b = order.size()==0;
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      if(b)
        SaveFloat64(stream,data[i]);
      else //case when data is well sorted
        SaveFloat64(stream,data[order[i]]);
      if( i< data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveTFloat64ListOfList(std::ostream& stream,std::string name,std::vector<TFloat64List> data, TFloat64List order) const
{
    if(data.size()>0){
    stream <<  "\""<<name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      stream << "[";
      for ( int ip=0; ip<data[order[i]].size(); ip++)
      {
        SaveFloat64(stream,data[order[i]][ip]);
        if ( ip< data[order[i]].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<data.size() - 1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveInt32Vector(std::ostream& stream,std::string name,std::vector<Int32> data, TFloat64List order) const
{
  Bool b = order.size()==0;
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      if(b)
        stream <<  data[i];
      else //case when data is well sorted
        stream <<  data[order[i]];
      if( i< data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveStringVector(std::ostream& stream,std::string name,std::vector<std::string> data, TFloat64List order) const 
{
  Bool b = order.size()==0;
  if(data.size()>0){
    stream <<  "\""<< name<<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      if(b)
        stream << "\"" << data[i] << "\""; 
      else
        stream << "\"" << data[order[i]] << "\""; 
      if( i< data.size()-1) stream << ",";
    }
    stream << "]" ;
  }
}

void COperatorResult::SaveStringVectorOfVector(std::ostream& stream,std::string name,std::vector<std::vector<std::string>> data, TFloat64List order) const
{
if(data.size()>0){
    Bool b = order.size()==0;
    stream <<  "\""<< name <<"\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      std::vector<std::string> line_list;
      if(b)
        line_list = data[i];
      else
        line_list = data[order[i]];
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

void COperatorResult::SaveTContinuumIndexListVector(std::ostream& stream,std::string name,std::vector<CContinuumIndexes::TContinuumIndexList> data, TFloat64List order) const
{
  if(data.size()>0){
    stream <<  "\""<<name << "Color\" : [";
    for ( int i=0; i<data.size(); i++)
    {
      stream << "[";
      for(Int32 kci=0; kci<data[order[i]].size(); kci++)
      {
	      SaveFloat64(stream,data[order[i]][kci].Color);
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
      for(Int32 kci=0; kci<data[order[i]].size(); kci++)
      {
	      SaveFloat64(stream,data[order[i]][kci].Break);	
        if ( kci<data[order[i]].size() - 1) stream << ",";
      }
      stream << "]";
      if ( i<data.size() - 1) stream << ",";
    }
    stream << "]" ;
  }
}

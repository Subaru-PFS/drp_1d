#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/rule.h>

#include <cstdarg>

using namespace NSEpic;
using namespace std;

/**
 * \brief Sets itself as disabled per default.
 **/
CRule::CRule ( )
{
  Enabled = false;
}

/**
 * Checks the elements and corrects if necessary.
 **/
void CRule::Apply( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  m_LineModelElements = LinemodelElements;
  if ( ! Enabled )
    {
      return;
    }
  if ( Check( LinemodelElements ) )
    {
      return;
    }
  Correct( LinemodelElements );
}

/**
 * \brief Returns the first index of m_Elements where calling the element's FindElementIndex method with LineCatalogIndex argument does not return -1.
 **/
Int32 CRule::FindElementIndex( Int32 LineCatalogIndex )
{
  Int32 idx = -1;
  for( UInt32 iElts=0; iElts<m_LineModelElements.size(); iElts++ )
    {
      if( m_LineModelElements[iElts]->FindElementIndex( LineCatalogIndex )!=-1 )
	{
	  idx = iElts;
	  break;
        }
    }
  return idx;
}

/**
 * \brief Returns the first index of m_Elements where calling the element's FindElementIndex method with LineTagStr argument does not return -1.
 **/
Int32 CRule::FindElementIndex( std::string LineTagStr, Int32 linetype )
{
  Int32 idx = -1;
  for( UInt32 iElts=0; iElts<m_LineModelElements.size(); iElts++ )
    {
      //if( m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetype){
      //      continue;
      //  }
      if( m_LineModelElements[iElts]->FindElementIndex( LineTagStr )!=-1 )
	{
	  idx = iElts;
	  break;
        }
    }
  return idx;
}

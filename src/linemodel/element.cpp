#include <epic/redshift/linemodel/element.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>


using namespace NSEpic;

CLineModelElement::CLineModelElement()
{
    m_OutsideLambdaRange = true;
    m_OutsideLambdaRangeOverlapThreshold = 0.1;
    //example: 0.1 means 10% of the line is allowed to be outside the spectrum with the line still considered inside the lambda range
}

CLineModelElement::~CLineModelElement()
{
}

Int32 CLineModelElement::FindElementIndex(Int32 LineCatalogIndex)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_LineCatalogIndexes.size(); iElts++ )
    {
        if(m_LineCatalogIndexes[iElts] == LineCatalogIndex){
            idx = iElts;
            break;
        }
    }

    return idx;
}

Int32 CLineModelElement::GetSize()
{
    return (Int32)m_LineCatalogIndexes.size();
}

bool CLineModelElement::IsOutsideLambdaRange()
{
    return m_OutsideLambdaRange;
}

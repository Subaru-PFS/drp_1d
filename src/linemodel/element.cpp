#include <epic/redshift/linemodel/element.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>


using namespace NSEpic;

CLineModelElement::CLineModelElement()
{
    m_OutsideLambdaRange = true;
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

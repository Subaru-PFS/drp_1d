#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_
#define _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <string>

namespace __NS__
{

class CTemplate : public CSpectrum
{

public:

    enum ECategory
    {
        nCategory_Emission = 0,
        nCategory_Galaxy,
        nCategory_Star,
        nCategory_Qso,
        nCategory_Count,
        nCategory_None = -1
    };

    CTemplate();
    CTemplate( const char* name, ECategory category );
    ~CTemplate();

    ECategory           GetCategory() const;
    const std::string&  GetName() const;

    static const char*  GetCategoryName( ECategory cat );

private:

    ECategory   m_Category;
    std::string m_Name;
};


}

#endif

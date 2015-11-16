#ifndef _REDSHIFT_RAY_CATALOG_
#define _REDSHIFT_RAY_CATALOG_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/ray/ray.h>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * /ingroup Redshift
 * Line catalog allow to store multiple lines description in a single text file.
 *
 * - Each line of the file represent a single Line
 * - Each line begenning with a # is a comment, and is skipped by the parser
 * - Format for each line is as follow:
 *        [Position in agstrum]   [Name of the line]                   [A/E]       [W/S]
 *        ex: 10320   [SII]                   E       W
 */
class CRayCatalog : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CRayCatalog )

public:

    typedef std::vector<CRay> TRayVector;

    CRayCatalog();
    ~CRayCatalog();


    Bool Add( const CRay& r );
    Bool Load( const char* filePath );
    const TRayVector& GetList() const;
    const TRayVector GetFilteredList(Int32 typeFilter = -1, Int32 forceFilter=-1) const;
    void Sort();
    void ConvertVacuumToAir();
private:

    TRayVector m_List;

};


}

#endif

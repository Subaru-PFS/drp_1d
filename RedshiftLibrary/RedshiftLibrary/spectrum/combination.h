#ifndef _REDSHIFT_SPECTRUM_COMBINATION_
#define _REDSHIFT_SPECTRUM_COMBINATION_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/spectrum/spectrum.h>


namespace NSEpic
{

class CSpectrumAxis;
class CMask;

/**
 * \ingroup Redshift
 */
class CSpectrumCombination
{

public:

    CSpectrumCombination( );
    ~CSpectrumCombination( );

    Int32 Combine(std::vector<std::shared_ptr<CSpectrum>>, CSpectrum &);

private:

};


}

#endif

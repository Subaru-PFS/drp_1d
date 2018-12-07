#ifndef _REDSHIFT_CONTINUUM_WAVELETSDF_
#define  _REDSHIFT_CONTINUUM_WAVELETSDF_

#include <RedshiftLibrary/continuum/continuum.h>

#include <RedshiftLibrary/common/datatypes.h>

namespace NSEpic
{

class CSpectrumFluxAxis;

class CContinuumDF : public CContinuum
{

public:
	explicit CContinuumDF( std::string binPath );
	CContinuumDF();
	~CContinuumDF();

	Bool RemoveContinuum ( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis);

	Bool m_optclus = false;

private:

   float* estimatedBaseline;
   float* estimatedBaseline_extd;
   float* extendedData_;


   void mirror_( const CSpectrum& s, UInt32 nb, float*  array);
   std::string temporaryPath ( const CSpectrum& s, UInt32 nall);
   std::string m_dfBinPath;

};

}

#endif

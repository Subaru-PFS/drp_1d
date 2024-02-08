// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_SPECTRUM_LSFFACTORY_
#define _REDSHIFT_SPECTRUM_LSFFACTORY_

#define LSFFactory (CLSFFactory::GetInstance())
//#define LSFArgsFactory (CLSFArgsFactory::GetInstance())

#include <map>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/spectrum/LSFVariableWidth.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"

namespace NSEpic {
class CLSF;
// the idea  is to define one single factory instance accross app instances
// now CLSFFactory, as defined here, "multiplexes" between the different LSFs
class CLSFFactory : public CSingleton<CLSFFactory> {

public:
  using CreateLSFFn = std::shared_ptr<CLSF> (*)(
      const std::shared_ptr<const TLSFArguments> &args);

  void Register(const std::string &name, CreateLSFFn fn_makeLSF);

  std::shared_ptr<CLSF>
  Create(const std::string &name,
         const std::shared_ptr<const TLSFArguments> &args = nullptr) {
    return m_FactoryMap.at(name)(args);
  };

private:
  friend class CSingleton<CLSFFactory>;

  CLSFFactory();
  ~CLSFFactory() = default;

  typedef std::map<std::string, CreateLSFFn> FactoryMap;
  FactoryMap m_FactoryMap;
};
inline CLSFFactory::CLSFFactory() {
  Register("gaussianConstantWidth", &CLSFGaussianConstantWidth::make_LSF);
  Register("gaussianConstantResolution",
           &CLSFGaussianConstantResolution::make_LSF);
  Register("gaussianNISPSIM2016", &CLSFGaussianNISPSIM2016::make_LSF);
  Register("gaussianNISPVSSPSF201707", &CLSFGaussianNISPVSSPSF201707::make_LSF);
  Register("gaussianVariableWidth", &CLSFGaussianVariableWidth::make_LSF);
}

inline void CLSFFactory::Register(const std::string &name,
                                  CreateLSFFn fn_makeLSF) {
  m_FactoryMap[name] = fn_makeLSF;
}
} // namespace NSEpic
#endif

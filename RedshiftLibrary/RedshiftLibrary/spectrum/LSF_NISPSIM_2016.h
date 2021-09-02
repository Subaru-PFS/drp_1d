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
#ifndef _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPSIM2016_COMPONENTS_
#define _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPSIM2016_COMPONENTS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/ray/lineprofile.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CLSFGaussianNISPSIM2016 : public CLSF
  {
    public:
        CLSFGaussianNISPSIM2016();

        Float64             GetWidth(Float64 lambda) const override;
        bool                IsValid() const override;
        void                SetSourcesizeDispersion(Float64 sigma);

        static std::shared_ptr<CLSF> make_LSF(const std::shared_ptr<const TLSFArguments>& args);
  };
inline
std::shared_ptr<CLSF>   CLSFGaussianNISPSIM2016::make_LSF(const std::shared_ptr<const TLSFArguments>& args)
{
     return std::make_shared<CLSFGaussianNISPSIM2016>();
}

}
#endif
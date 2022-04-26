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
#ifndef _REDSHIFT_SPECTRUM_LSF_
#define _REDSHIFT_SPECTRUM_LSF_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
namespace NSEpic {
class CLineProfile;

struct TLSFArguments {
  // std::string type;
  virtual ~TLSFArguments(){};
  TLSFArguments() = default;
  TLSFArguments(const TLSFArguments &other) = default;
  TLSFArguments(TLSFArguments &&other) = default;
  TLSFArguments &operator=(const TLSFArguments &other) = default;
  TLSFArguments &operator=(TLSFArguments &&other) = default;
};

struct TLSFGaussianVarWidthArgs : virtual TLSFArguments {
  // std::string type = "GaussianVariableWidth";
  TFloat64List lambdas;
  TFloat64List width;
  TLSFGaussianVarWidthArgs(TFloat64List _lambdas, TFloat64List _width)
      : lambdas(_lambdas), width(_width) {}
};

struct TLSFGaussianConstantWidthArgs : virtual TLSFArguments {
  Float64 width;
  TLSFGaussianConstantWidthArgs(
      const std::shared_ptr<const CParameterStore> &parameterStore)
      : width(parameterStore->GetScoped<Float64>("LSF.width")) {}

  TLSFGaussianConstantWidthArgs(Float64 _width) : width(_width) {}
};

struct TLSFGaussianConstantResolutionArgs : virtual TLSFArguments {
  Float64 resolution;
  TLSFGaussianConstantResolutionArgs(
      const std::shared_ptr<const CParameterStore> &parameterStore)
      : resolution(parameterStore->GetScoped<Float64>("LSF.resolution")) {}

  TLSFGaussianConstantResolutionArgs(Float64 _resolution)
      : resolution(_resolution) {}
};

struct TLSFGaussianNISPVSSPSF201707Args : virtual TLSFArguments {
  Float64 sourcesize;
  TLSFGaussianNISPVSSPSF201707Args(
      const std::shared_ptr<const CParameterStore> &parameterStore)
      : sourcesize(parameterStore->GetScoped<Float64>("LSF.sourcesize")) {}
};
/**
 * \ingroup Redshift
 */
class CLSF {

public:
  enum TLSFType {
    GaussianConstantWidth,
    GaussianConstantResolution,
    GaussianNISPSIM2016,
    GaussianNISPVSSPSF201707,
    GaussianVariableWidth
  };
  CLSF(TLSFType name);
  CLSF(TLSFType name, CLineProfile_ptr &&profile);
  virtual ~CLSF() = default;

  CLSF(const CLSF &other) = default;
  CLSF(CLSF &&other) = default;
  CLSF &operator=(const CLSF &other) = default;
  CLSF &operator=(CLSF &&other) = default;
  // GetWidth requires observed wavelength, not restframe
  virtual Float64 GetWidth(Float64 lambda) const = 0;
  virtual bool checkAvailability(Float64 lambda) const {
    return true;
  }; // default to true
  virtual bool IsValid() const = 0;
  Float64 GetLineProfile(Float64 lambda, Float64 lambda0 = 0.) const;
  Float64 GetLineProfile(Float64 lambda, Float64 lambda0, Float64 sigma0) const;
  TFloat64List getRestFrameProfileVector(Float64 lambda0_rest, Float64 z) const;

  const CLineProfile &GetProfile() const;
  const TLSFType m_name;

protected:
  CLineProfile_ptr m_profile;
};

} // namespace NSEpic
#endif

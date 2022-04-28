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
#ifndef _REDSHIFT_RAY_AIRVACUUM_
#define _REDSHIFT_RAY_AIRVACUUM_

#include "RedshiftLibrary/common/datatypes.h"

#include <vector>

namespace NSEpic {
typedef struct conversionParams {
  const Float64 a;
  const Float64 b1;
  const Float64 b2;
  const Float64 c1;
  const Float64 c2;
  const std::string name;
  conversionParams(const Float64 _a, const Float64 _b1, const Float64 _b2,
                   const Float64 _c1, const Float64 _c2,
                   const std::string &_name = "Undefined")
      : a(_a), b1(_b1), b2(_b2), c1(_c1), c2(_c2), name(_name) {}
} conversionParams;
static const conversionParams Morton2000Params{
    8.34254e-5, 2.406147e-2, 1.5998e-4, 130.0, 38.9, "Morton2000"};
static const conversionParams VacCiddor1996Params{
    0.0, 5.792105E-2, 1.67917E-3, 238.0185, 57.362, "VacCiddor1996"};
static const conversionParams PeckReeder1972Params{
    0.0, 5.791817E-2, 1.67909E-3, 238.0185, 57.362, "PeckReeder1972"};
static const conversionParams Edlen1966Params{
    8.34213E-5, 2.406030E-2, 1.5997E-4, 130.0, 38.9, "Edlen1966"};
static const conversionParams Edlen1953Params{
    6.4328E-5, 2.94981E-2, 2.5540E-4, 146.0, 41.0, "Edlen1953"};
class CAirVacuum {
public:
  CAirVacuum(const conversionParams &params)
      : m_a(params.a), m_b1(params.b1), m_b2(params.b2), m_c1(params.c1),
        m_c2(params.c2), m_name(params.name){};

  // rule of 5 defaults
  CAirVacuum(const CAirVacuum &other) = default;
  CAirVacuum(CAirVacuum &&other) = default;
  CAirVacuum &operator=(const CAirVacuum &other) = default;
  CAirVacuum &operator=(CAirVacuum &&other) = default;
  virtual ~CAirVacuum() = default;

  virtual TFloat64List AirToVac(const TFloat64List &waveAir) const;
  virtual TFloat64List VacToAir(const TFloat64List &waveVac) const;

  TFloat64List AirRefractiveIndex(const TFloat64List &waveVac) const;

protected:
  virtual void CheckWaveRange(const TFloat64List &wave) const {};

  const Float64 m_a;
  const Float64 m_b1;
  const Float64 m_b2;
  const Float64 m_c1;
  const Float64 m_c2;
  const std::string m_name;
};
class CAirVacuumConverter {
public:
  static std::shared_ptr<CAirVacuum> Get(const std::string &ConverterName);
};

class CAirVacEdlen1953 : public CAirVacuum {
public:
  CAirVacEdlen1953() : CAirVacuum(Edlen1953Params){};
  void CheckWaveRange(const TFloat64List &wave) const override;
};

class CAirVacEdlen1966 : public CAirVacuum {
public:
  CAirVacEdlen1966() : CAirVacuum(Edlen1966Params){};
  void CheckWaveRange(const TFloat64List &wave) const override;
};

class CAirVacPeckReeder1972 : public CAirVacuum {
public:
  CAirVacPeckReeder1972() : CAirVacuum(PeckReeder1972Params){};
  void CheckWaveRange(const TFloat64List &wave) const override;
};

class CAirVacCiddor1996 : public CAirVacuum {
public:
  CAirVacCiddor1996() : CAirVacuum(VacCiddor1996Params){};
  void CheckWaveRange(const TFloat64List &wave) const override;
};

/*
 *   Morton (2000) actually comes from Birch and Downs (1994, Metrologia, 31,
 * 315)
 */
class CAirVacMorton2000 : public CAirVacuum {
public:
  CAirVacMorton2000() : CAirVacuum(Morton2000Params){};
  void CheckWaveRange(const TFloat64List &wave) const override;
  TFloat64List AirToVac(const TFloat64List &waveAir) const override;
};

} // namespace NSEpic

#endif
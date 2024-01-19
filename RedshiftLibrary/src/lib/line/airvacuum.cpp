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
#include "RedshiftLibrary/line/airvacuum.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

std::shared_ptr<CAirVacuum>
CAirVacuumConverter::Get(const std::string &converterName) {
  if (converterName == "edlen1953")
    return (std::make_shared<CAirVacEdlen1953>());
  else if (converterName == "edlen1966")
    return (std::make_shared<CAirVacEdlen1966>());
  else if (converterName == "PeckReeder1972")
    return (std::make_shared<CAirVacPeckReeder1972>());
  else if (converterName == "ciddor1996")
    return (std::make_shared<CAirVacCiddor1996>());
  else if (converterName == "morton2000")
    return (std::make_shared<CAirVacMorton2000>());
  else {
    THROWG(INTERNAL_ERROR,
           Formatter() << "unknown air->vacuum conversion: " << converterName);
  }
}

TFloat64List CAirVacuum::VacToAir(const TFloat64List &waveVac) const {
  auto sz = waveVac.size();

  CheckWaveRange(waveVac);
  Log.LogDetail(Formatter() << "CAirVacuum::VacToAir converting wavelengths "
                               "from vacuum to air using translation from "
                            << m_name);

  const TFloat64List n = AirRefractiveIndex(waveVac);

  TFloat64List waveAir(sz);
  // performs waveVacList/nList
  std::transform(waveVac.begin(), waveVac.end(), n.begin(), waveAir.begin(),
                 std::divides<Float64>());

  return (waveAir);
}

TFloat64List CAirVacuum::AirToVac(const TFloat64List &waveAir) const {
  const Float64 minprecision = 1e-9;
  const Int32 maxiter = 30;

  TFloat64List waveVac = waveAir;
  TFloat64List delta, deltaAir;
  Int32 iter = 0;
  Float64 precision;

  CheckWaveRange(waveAir);
  Log.LogDetail(Formatter()
                << "CAirVacuum::AirToVac converting wavelengths from air to "
                   "vacuum by iterating the inverse translation from "
                << m_name);

  do {
    // delta = VacToAir(waveVac) - waveAir;
    delta = VacToAir(waveVac);
    std::transform(delta.begin(), delta.end(), waveAir.begin(), delta.begin(),
                   std::minus<Float64>());

    // waveVac = waveVac - delta
    std::transform(waveVac.begin(), waveVac.end(), delta.begin(),
                   waveVac.begin(), std::minus<Float64>());

    precision =
        *std::max_element(delta.begin(), delta.end(), [](Float64 x, Float64 y) {
          return (abs(x) < abs(y));
        });
    if (precision <= minprecision)
      break;

  } while (++iter < maxiter);

  if (iter == maxiter) {
    Flag.warning(WarningCode::AIR_VACCUM_REACHED_MAX_ITERATIONS,
                 Formatter()
                     << "CAirVacuum::" << __func__
                     << " reach max iteration, with precision: " << precision);
  }

  return (waveVac);
}

/*
 *  direct inversion following N.E Piskunov
 * (https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)
 *
 *   n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) +
 * 0.0001599740894897 / (38.92568793293 - s2), where s = 1E4 / λair and the
 * conversion is: λvac = λair * n.
 */
TFloat64List CAirVacMorton2000::AirToVac(const TFloat64List &waveAir) const {
  CheckWaveRange(waveAir);
  Log.LogDetail(
      Formatter()
      << "CAirVacuum::AirToVac converting wavelengths from air to vacuum using "
         "an approximation of the inverse translation from "
      << m_name);

  auto refractiveindex = [](Float64 w) {
    Float64 s = 1e8 / (w * w);
    return (1. + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s) +
            0.0001599740894897 / (38.92568793293 - s));
  };

  TFloat64List n(waveAir.size());

  std::transform(waveAir.begin(), waveAir.end(), n.begin(), refractiveindex);

  TFloat64List waveVac(waveAir.size());
  std::transform(waveAir.begin(), waveAir.end(), n.begin(), waveVac.begin(),
                 std::multiplies<Float64>());

  return (waveVac);
}

/*
 *  computes air refraction index from wavelength expressed in Angstroem
 */
TFloat64List CAirVacuum::AirRefractiveIndex(const TFloat64List &waveVac) const {
  auto refractiveindex = [this](Float64 w) {
    Float64 s = 1e8 / (w * w);
    return (1. + m_a + m_b1 / (m_c1 - s) + m_b2 / (m_c2 - s));
  };

  TFloat64List n(waveVac.size());

  std::transform(waveVac.begin(), waveVac.end(), n.begin(), refractiveindex);

  return (n);
}

void CAirVacEdlen1953::CheckWaveRange(const TFloat64List &wave) const {
  auto it = std::min_element(wave.begin(), wave.end());
  if (it == wave.end())
    THROWG(INTERNAL_ERROR, " Empty vector");

  if (*it < 2000.0) {
    THROWG(INTERNAL_ERROR, " some wavelengths "
                           "are below 2000 Angstroem");
  }
}

void CAirVacEdlen1966::CheckWaveRange(const TFloat64List &wave) const {
  auto it = std::min_element(wave.begin(), wave.end());
  if (it == wave.end())
    THROWG(INTERNAL_ERROR, "Empty vector");

  if (*it < 2000.0) {
    THROWG(INTERNAL_ERROR, "some wavelengths "
                           "are below 2000 Angstroem");
  }
}

void CAirVacPeckReeder1972::CheckWaveRange(const TFloat64List &wave) const {
  auto it = std::min_element(wave.begin(), wave.end());
  if (it == wave.end())
    THROWG(INTERNAL_ERROR, "Empty vector");

  if (*it < 3000.0) {
    THROWG(INTERNAL_ERROR, "some "
                           "wavelengths are below 2300 Angstroem");
  }

  if (*std::max_element(wave.begin(), wave.end()) > 19000.) {
    THROWG(INTERNAL_ERROR, "some "
                           "wavelengths are above 16900 Angstroem");
  }
}

void CAirVacCiddor1996::CheckWaveRange(const TFloat64List &wave) const {
  auto it = std::min_element(wave.begin(), wave.end());
  if (it == wave.end())
    THROWG(INTERNAL_ERROR, "Empty vector");

  if (*it < 3000.0) {
    THROWG(INTERNAL_ERROR, "some "
                           "wavelengths are below 3000 Angstroem");
  }

  if (*std::max_element(wave.begin(), wave.end()) > 19000.) {
    THROWG(INTERNAL_ERROR, "some "
                           "wavelengths are above 16900 Angstroem");
  }
}

void CAirVacMorton2000::CheckWaveRange(const TFloat64List &wave) const {
  auto it = std::min_element(wave.begin(), wave.end());
  if (it == wave.end())
    THROWG(INTERNAL_ERROR, "Empty vector");

  if (*it < 3000.0) {
    THROWG(INTERNAL_ERROR, "some wavelengths "
                           "are below 3000 Angstroem");
  }

  if (*std::max_element(wave.begin(), wave.end()) > 19000.) {
    THROWG(INTERNAL_ERROR, "some wavelengths "
                           "are above 16900 Angstroem");
  }
}

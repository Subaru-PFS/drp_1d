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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/airvacuum.h"


#include <boost/test/unit_test.hpp>


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(CAirVacuum)

Float64 precision = 1e-2;
TFloat64List lambda = {6330., 16000.};

BOOST_AUTO_TEST_CASE(Edlen1953)
{
  
  BOOST_CHECK_NO_THROW(CAirVacEdlen1953 converter);

  CAirVacEdlen1953 converter; 
  TFloat64List n = converter.AirRefractiveIndex(lambda);
  BOOST_CHECK_CLOSE((n[0]-1.)*1e8, 27651.65, precision);
  BOOST_CHECK_CLOSE((n[1]-1.)*1e8, 27320.098, precision);

  TFloat64List lambdaAir = converter.VacToAir(lambda);
  BOOST_CHECK_CLOSE(lambdaAir[0], 6328.25, precision);
  BOOST_CHECK_CLOSE(lambdaAir[1], 15995.63, precision);

  TFloat64List lambdaVac = converter.AirToVac(lambdaAir);
  BOOST_CHECK_CLOSE(lambdaVac[0], lambda[0], precision);
  BOOST_CHECK_CLOSE(lambdaVac[1], lambda[1], precision);

}

BOOST_AUTO_TEST_CASE(Edlen1966)
{
  
  BOOST_CHECK_NO_THROW(CAirVacEdlen1966 converter);

  CAirVacEdlen1966 converter; 
  TFloat64List n = converter.AirRefractiveIndex(lambda);
  BOOST_CHECK_CLOSE((n[0]-1.)*1e8, 27651.74, precision);
  BOOST_CHECK_CLOSE((n[1]-1.)*1e8, 27321.24, precision);

  TFloat64List lambdaAir = converter.VacToAir(lambda);
  BOOST_CHECK_CLOSE(lambdaAir[0], 6328.25, precision);
  BOOST_CHECK_CLOSE(lambdaAir[1], 15995.63, precision);

  TFloat64List lambdaVac = converter.AirToVac(lambdaAir);
  BOOST_CHECK_CLOSE(lambdaVac[0], lambda[0], precision);
  BOOST_CHECK_CLOSE(lambdaVac[1], lambda[1], precision);

}

BOOST_AUTO_TEST_CASE(PeckReeder1972)
{
  
  BOOST_CHECK_NO_THROW(CAirVacPeckReeder1972 converter);

  CAirVacPeckReeder1972 converter; 
  TFloat64List n = converter.AirRefractiveIndex(lambda);
  BOOST_CHECK_CLOSE((n[0]-1.)*1e8, 27651.65, precision);
  BOOST_CHECK_CLOSE((n[1]-1.)*1e8, 27320.7267, precision);

  TFloat64List lambdaAir = converter.VacToAir(lambda);
  BOOST_CHECK_CLOSE(lambdaAir[0], 6328.25, precision);
  BOOST_CHECK_CLOSE(lambdaAir[1], 15995.63, precision);

  TFloat64List lambdaVac = converter.AirToVac(lambdaAir);
  BOOST_CHECK_CLOSE(lambdaVac[0], lambda[0], precision);
  BOOST_CHECK_CLOSE(lambdaVac[1], lambda[1], precision);

}

BOOST_AUTO_TEST_CASE(Ciddor1996)
{
  
  BOOST_CHECK_NO_THROW(CAirVacCiddor1996 converter);

  CAirVacCiddor1996 converter; 
  TFloat64List n = converter.AirRefractiveIndex(lambda);
  BOOST_CHECK_CLOSE((n[0]-1.)*1e8, 27653.02, precision);
  BOOST_CHECK_CLOSE((n[1]-1.)*1e8, 27322.08, precision);

  TFloat64List lambdaAir = converter.VacToAir(lambda);
  BOOST_CHECK_CLOSE(lambdaAir[0], 6328.25, precision);
  BOOST_CHECK_CLOSE(lambdaAir[1], 15995.63, precision);

  TFloat64List lambdaVac = converter.AirToVac(lambdaAir);
  BOOST_CHECK_CLOSE(lambdaVac[0], lambda[0], precision);
  BOOST_CHECK_CLOSE(lambdaVac[1], lambda[1], precision);

}

BOOST_AUTO_TEST_CASE(Morton2000)
{
  
  BOOST_CHECK_NO_THROW(CAirVacMorton2000 converter);

  CAirVacMorton2000 converter; 
  TFloat64List n = converter.AirRefractiveIndex(lambda);
  BOOST_CHECK_CLOSE((n[0]-1.)*1e8, 27653.099, precision);
  BOOST_CHECK_CLOSE((n[1]-1.)*1e8, 27322.58, precision);

  TFloat64List lambdaAir = converter.VacToAir(lambda);
  BOOST_CHECK_CLOSE(lambdaAir[0], 6328.25, precision);
  BOOST_CHECK_CLOSE(lambdaAir[1], 15995.63, precision);

  TFloat64List lambdaVac = converter.AirToVac(lambdaAir);
  BOOST_CHECK_CLOSE(lambdaVac[0], lambda[0], precision);
  BOOST_CHECK_CLOSE(lambdaVac[1], lambda[1], precision);

}
//}
//BOOST_AUTO_TEST_SUITE(CAirVacuumConverter)

BOOST_AUTO_TEST_CASE(Converter)
{
  BOOST_CHECK_THROW(CAirVacuumConverter::Get(""), GlobalException);
  BOOST_CHECK_NO_THROW(CAirVacuumConverter::Get("Morton2000"));

  auto converter1=CAirVacuumConverter::Get("Morton2000"); 
  CAirVacMorton2000 converter2; 
  
  TFloat64List n1 = converter1->AirRefractiveIndex(lambda);
  TFloat64List n2 = converter2.AirRefractiveIndex(lambda);
  BOOST_CHECK_EQUAL(n1[0], n2[0]);
  BOOST_CHECK_EQUAL(n1[1], n2[1]);

  TFloat64List lambdaAir1 = converter1->VacToAir(lambda);
  TFloat64List lambdaAir2 = converter2.VacToAir(lambda);
  BOOST_CHECK_EQUAL(lambdaAir1[0], lambdaAir2[0]);
  BOOST_CHECK_EQUAL(lambdaAir1[1], lambdaAir2[1]);

  TFloat64List lambdaVac1 = converter1->AirToVac(lambdaAir1);
  TFloat64List lambdaVac2 = converter2.AirToVac(lambdaAir1);
  BOOST_CHECK_EQUAL(lambdaVac1[0], lambdaVac2[0]);
  BOOST_CHECK_EQUAL(lambdaVac1[1], lambdaVac2[1]);

}

}

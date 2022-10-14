
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
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

const std::string jsonString1 =
    "{\"LSF\" : {\"LSFType\" : \"GaussianConstantWidth\" , \"width\" : 9.5}}";

const std::string jsonString2 =
    "{\"LSF\" : {\"LSFType\" : \"GaussianConstantResolution\" , \"resolution\" "
    ": 0.1}}";

const std::string jsonString3 =
    "{\"LSF\" : {\"LSFType\" : \"GaussianNISPSIM2016\" }}";

const std::string jsonString4 =
    "{\"LSF\" : {\"LSFType\" : \"GaussianNISPVSSPSF201707\" , \"sourcesize\" "
    ": 0.3}}";

const std::string jsonString5 =
    "{\"LSF\" : {\"LSFType\" : \"GaussianVariableWidth\"}}";

BOOST_AUTO_TEST_SUITE(myInputContext_test)

BOOST_AUTO_TEST_CASE(LSF_test) {
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> paramStore =
      std::make_shared<CParameterStore>(scopeStack);
  paramStore->FromString(jsonString1);

  MyInputContext ctx(paramStore);
  std::shared_ptr<CLSF> lsf = ctx.createLSF();
  BOOST_CHECK(lsf->m_name == CLSF::TLSFType::GaussianConstantWidth);

  paramStore->FromString(jsonString2);
  lsf = ctx.createLSF();
  BOOST_CHECK(lsf->m_name == CLSF::TLSFType::GaussianConstantResolution);

  paramStore->FromString(jsonString3);
  lsf = ctx.createLSF();
  BOOST_CHECK(lsf->m_name == CLSF::TLSFType::GaussianNISPSIM2016);

  paramStore->FromString(jsonString4);
  lsf = ctx.createLSF();
  BOOST_CHECK(lsf->m_name == CLSF::TLSFType::GaussianNISPVSSPSF201707);

  paramStore->FromString(jsonString5);
  lsf = ctx.createLSF();
  BOOST_CHECK(lsf->m_name == CLSF::TLSFType::GaussianVariableWidth);
}

BOOST_AUTO_TEST_SUITE_END()
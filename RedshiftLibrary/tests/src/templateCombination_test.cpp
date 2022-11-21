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
#include "RedshiftLibrary/operator/tplcombination.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "test-config.h"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(TemplateCombination)

BOOST_AUTO_TEST_CASE(inverseMatrix) {
  /*
      Int32 dim = 3;
      cov = gsl_matrix_alloc (dim, dim);
      TFloatList _cov(1.0, 0.6, 0.0,
                      0.0, 1.5, 1.0,
                      0.0, 1.0, 1.0 );
      Int32 count = 0;
      for(Int32 i = 0; i < dim; i++)
      {
          for(Int32 j = 0; j < dim; j++)
          {
              gsl_matrix_set (cov, i, j, _cov[count]);
              count++;
          }
      }

      //COperatorTplcombination::InverseMatrix(cov);
      TFloatList inv_cov(
              1,  -1.2,   1.2,
              0,   2.0,  -2.0,
              0,  -2.0,   3.0);
      count = 0;
      for (Int32 i = 0; i < dim; ++i)
      {
          for (Int32 j = 0; j < dim; ++j)
          {
              Float64 v = gsl_matrix_get(&inv.matrix,i,j));
              //BOOST_CHECK(v == inv_cov[count]);
              count++;
          }
      }
     */
}

BOOST_AUTO_TEST_SUITE_END()
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
#ifndef PYCONV_H
#define PYCONV_H
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"

namespace NSEpic {
class PC {
public:
  static void get(const TFloat64List &vec, double **data, int *size) {
    *data = const_cast<double *>(vec.data());
    *size = vec.size();
  }
  static void getasl(const TAxisSampleList &vec, double **data, int *size) {
    *data = const_cast<double *>(vec.data());
    *size = vec.size();
  }

  static void get(const TInt32List &vec, int **data, int *size) {
    *data = const_cast<int *>(vec.data());
    *size = vec.size();
  }
  static void get(const TFloat32List &vec, float **data, int *size) {
    *data = const_cast<float *>(vec.data());
    *size = vec.size();
  }

  static void get(const TBoolList &vec, short **data, int *size) {
    *data = (short *)malloc(sizeof(short) * vec.size());

    for (std::size_t i = 0; i < vec.size(); i++)
      (*data)[i] = vec[i];
    *size = vec.size();
  }

  static void get(const TMaskList &vec, unsigned char **data, int *size) {
    *data = const_cast<unsigned char *>(vec.data());
    *size = vec.size();
  }

  static void get(const TMaskList &vec, short **data, int *size) {
    *data = (short *)malloc(sizeof(short) * vec.size());
    *size = (int)vec.size();
    for (std::size_t i = 0; i < vec.size(); i++)
      (*data)[i] = (short)vec[i];
  }
};
} // namespace NSEpic
#endif

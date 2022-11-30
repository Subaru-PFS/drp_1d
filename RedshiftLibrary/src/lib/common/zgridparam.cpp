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
// ============================================================================using
// namespace std;
#include "RedshiftLibrary/common/zgridparam.h"
#include "RedshiftLibrary/common/vectorOperations.h"

using namespace NSEpic;

TFloat64List CZGridParam::getZGrid(bool logsampling) const {

  TFloat64Range range(zmin, zmax);
  if (!isnan(zcenter))
    return range.spanCenteredWindow(zcenter, logsampling, zstep);

  auto grid = logsampling ? range.SpreadOverLogZplusOne(zstep)
                          : range.SpreadOver(zstep);

  grid.back() = std::min(grid.back(), zmax); // cope with rounding errors

  return grid;
}

TFloat64List CZGridListParams::getZGrid(bool logsampling) const {

  if (m_zparams.empty())
    return TFloat64List();

  TFloat64List zgrid = m_zparams.front().getZGrid(logsampling);

  // if needed create and insert 2nd pass subgrids into main grid
  for (Int32 i = 1; i < m_zparams.size(); i++) {
    // create a centered sub grid around Zcand
    TFloat64List subgrid = m_zparams[i].getZGrid(logsampling);
    // insert it into main gridgetZGrid
    insertSubgrid(subgrid, zgrid);
  }
  return zgrid;
}

std::tuple<Int32, Int32> CZGridListParams::insertSubgrid(TFloat64List &subgrid,
                                                         TFloat64List &zgrid) {
  const Float64 epsilon = 1E-8;
  TFloat64Range range_epsilon = {subgrid.front() - epsilon,
                                 subgrid.back() + epsilon};
  range_epsilon.IntersectWith(zgrid);
  Int32 imin = -1;
  Int32 imax = -1;
  bool b = range_epsilon.getClosedIntervalIndices(zgrid, imin, imax, false);
  if (!b) // range not included in the main range
    THROWG(INTERNAL_ERROR, "range not inside base grid ");

  // deal with subgrid front or end samples when equal to coarse grid samples
  auto subgrid_start = subgrid.cbegin();
  auto subgrid_end = subgrid.cend();
  if (range_epsilon.GetBegin() == zgrid.front()) {
    ++imin;
    ++subgrid_start;
    subgrid.front() = zgrid.front();
  }
  if (range_epsilon.GetEnd() == zgrid.back()) {
    --imax;
    --subgrid_end;
    subgrid.back() = zgrid.back();
  }

  Int32 ndup = imax - imin + 1;
  insertWithDuplicates(zgrid, imin, subgrid_start, subgrid_end, ndup);

  return std::make_tuple(imin, ndup);
}

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

class CZGridParam {
public:
  CZGridParam() = default;
  CZGridParam(const TFloat64Range &range, Float64 step, Float64 center = NAN)
      : zmin(range.GetBegin()), zmax(range.GetEnd()), zstep(step),
        zcenter(center){};

  TFloat64List getZGrid(bool logsampling) const;

  Float64 zmin = NAN;
  Float64 zmax = NAN;
  Float64 zstep = NAN;
  Float64 zcenter = NAN;
};

typedef std::vector<CZGridParam> TZGridListParams;

class CZGridListParams {
public:
  //  CZGridListParams() = default;
  CZGridListParams(const TZGridListParams &params)
      : m_zparams(params){};
  CZGridListParams(TZGridListParams &&params)
      : m_zparams(std::move(params)){};

  CZGridParam &operator[](const Int32 i) { return m_zparams[i]; };
  const CZGridParam &operator[](const Int32 i) const { return m_zparams[i]; };
  Int32 size() const { return m_zparams.size(); }

  TFloat64List getZGrid(bool logsampling) const;

  // returns insertion index & mumber of overlapped indices
  static std::tuple<Int32, Int32> insertSubgrid(const TFloat64List &subgrid,
                                                TFloat64List &zgrid);

private:
  TZGridListParams m_zparams;
};

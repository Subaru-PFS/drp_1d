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
typedef struct ZGridParameters {
  ZGridParameters(const TFloat64Range &range, Float64 step,
                  Float64 center = NAN)
      : zmin(range.GetBegin()), zmax(range.GetEnd()), zstep(step),
        zcenter(center){};
  ZGridParameters(Float64 min, Float64 max, Float64 step, Float64 center = NAN)
      : zmin(min), zmax(max), zstep(step), zcenter(center){};

  ZGridParameters() = default;
  Float64 zmin = NAN;
  Float64 zmax = NAN;
  Float64 zstep = NAN;
  Float64 zcenter = NAN;
} ZGridParameters;

typedef struct TPdf {
  TPdf(const TFloat64List &zgrid_, const TFloat64List &probaLog_)
      : zgrid(zgrid_), probaLog(probaLog_){};
  TPdf() = default;
  TFloat64List zgrid;
  TFloat64List probaLog;
} TPdf;

typedef std::vector<ZGridParameters> TZGridListParams;

class CLogZPdfResult : public COperatorResult {
public:
  CLogZPdfResult();

  CLogZPdfResult(const TFloat64List &redshifts,
                 const TZGridListParams &zparams);

  CLogZPdfResult(const TZGridListParams &zparamList,
                 const TFloat64List &valproba_mixed); // for ML learning

  bool isZGridCoherent() const;
  const TFloat64List buildLogZPdfGrid(bool logsampling, bool fine = false,
                                      Int32 idx = 0) const;
  const TFloat64List buildLogMixedZPdfGrid(bool logsampling) const;

  const CLogZPdfResult getLogZPdf_fine(bool logsampling) const;
  static void interpolateLargeGridOnFineGrid(
      const TFloat64List &coarseGrid, const TFloat64List &fineGrid,
      const TFloat64List &entityLargeGrid, TFloat64List &entityFineGrid);

  static TPdf getFineGridPdf(
      const TZGridListParams &pdf_params,
      //			     const Float64* pdf_probaLog, Int32 n);
      const TFloat64List &pdf_probaLog, bool logsampling = true);
  static TPdf getMixedGridPdf(const TZGridListParams &pdf_params,
                              const TFloat64List &pdf_probaLog,
                              bool logsampling = true);

  TFloat64List Redshifts;
  TFloat64List valProbaLog;
  Float64 valEvidenceLog = NAN;
  Float64 valMargEvidenceLog = NAN;

  TFloat64List zcenter;
  TFloat64List zmin;
  TFloat64List zmax;
  TFloat64List zstep;

  TZGridListParams zparams;

private:
  void setZGridParams(const TZGridListParams &paramList);
};

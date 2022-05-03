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

class CLineModelSolution : public COperatorResult {
public:
  CLineModelSolution();

  TFloat64List ElementId; // id of the linemodel element it is part of
  TFloat64List Amplitudes;
  std::vector<CLine> Lines;
  TFloat64List AmplitudesUncertainties;       // noise sigma
  TFloat64List FittingError; // ModelLeastSquare error under each line
  TFloat64List
      CenterContinuumFlux; // Continuum flux value at the center of each line
  TFloat64List ContinuumError;        // Continuum error value for each line
  TFloat64List Sigmas;                // width for each line
  TFloat64List Fluxs;                 // Flux for each line
  TFloat64List FluxErrors;            // Flux error for each line
  TFloat64List FluxDirectIntegration; // Flux obtained by direct integration for
                                      // each line
  TFloat64List FluxDirectIntegrationError; // Flux obtained by direct
                                           // integration for each line
  TInt32List lineId;

  Float64 snrHa;
  Float64 lfHa;
  Float64 snrOII;
  Float64 lfOII;
  Int32 NLinesAboveSnrCut;

  TFloat64List LambdaObs; // observed position in Angstrom
  TFloat64List Velocity;  // dispersion velocity in km/s
  TFloat64List Offset;    // line offset in km/s
  TBoolList OutsideLambdaRange;
  std::vector<TInt32Range> fittingIndexRange;
  TStringList fittingGroupInfo;

  Float64 LyaWidthCoeff;
  Float64 LyaAlpha;
  Float64 LyaDelta;

  Float64 AbsorptionVelocity;
  Float64 EmissionVelocity;
  Float64 Redshift;

  TFloat64List continuum_pCoeff0;
  TFloat64List continuum_pCoeff1;
  TFloat64List continuum_pCoeff2;

  Int32 nDDL;
  void fillLineIds();
};

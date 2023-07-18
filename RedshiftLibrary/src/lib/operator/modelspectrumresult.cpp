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
#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iomanip> // std::setprecision
#include <string>

#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;

/**
 * \brief Sets the model to CSpectrum ( spc ).
 **/
CModelSpectrumResult::CModelSpectrumResult(const CSpectrum &spc,
                                           const std::string &obsId)
    : COperatorResult("CModelSpectrumResult") {
  ModelLambda[obsId] = spc.GetSpectralAxis().GetSamplesVector();
  ModelFlux[obsId] = spc.GetFluxAxis().GetSamplesVector();
  // probably can add model params as class variable here..
}

CModelSpectrumResult::CModelSpectrumResult(CSpectrum &&spc,
                                           const std::string &obsId)
    : COperatorResult("CModelSpectrumResult") {
  ModelLambda[obsId] = std::move(spc.GetSpectralAxis().GetSamplesVector());
  ModelFlux[obsId] = std::move(spc.GetFluxAxis().GetSamplesVector());
}

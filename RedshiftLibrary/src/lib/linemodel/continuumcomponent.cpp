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

#include "RedshiftLibrary/linemodel/continuumcomponent.h"

using namespace NSEpic;

const std::unordered_map<std::string, TContinuumComponent::EContinuumComponent>
    TContinuumComponent::stringToEnum = {
        {"tplFit", TContinuumComponent::EContinuumComponent::tplFit},
        {"tplFitAuto", TContinuumComponent::EContinuumComponent::tplFitAuto},
        {"noContinuum", TContinuumComponent::EContinuumComponent::noContinuum},
        {"fromSpectrum",
         TContinuumComponent::EContinuumComponent::fromSpectrum},
        {"powerLaw", TContinuumComponent::EContinuumComponent::powerLaw},
        {"powerLawAuto",
         TContinuumComponent::EContinuumComponent::powerLawAuto}};

const std::unordered_map<TContinuumComponent::EContinuumComponent, std::string>
    TContinuumComponent::enumToString = {
        {TContinuumComponent::EContinuumComponent::tplFit, "tplFit"},
        {TContinuumComponent::EContinuumComponent::tplFitAuto, "tplFitAuto"},
        {TContinuumComponent::EContinuumComponent::noContinuum, "noContinuum"},
        {TContinuumComponent::EContinuumComponent::fromSpectrum,
         "fromSpectrum"},
        {TContinuumComponent::EContinuumComponent::powerLaw, "powerLaw"},
        {TContinuumComponent::EContinuumComponent::powerLawAuto,
         "powerLawAuto"}};

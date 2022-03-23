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
#ifndef _REDSHIFT_RAY_LINETAGS_
#define _REDSHIFT_RAY_LINETAGS_
#include <string>
namespace NSEpic
{

class linetags{

public:
    const char* halpha_abs = "HalphaA";
    const char* hbeta_abs = "HbetaA";
    const char* hgamma_abs = "HgammaA";
    const char* hdelta_abs = "HdeltaA";
    const char* h8_abs = "H8A";
    const char* h9_abs = "H9A";
    const char* h10_abs = "H10A";
    const char* h11_abs = "H11A";


    const char* halpha_em = "Halpha";
    const char* hbeta_em = "Hbeta";
    const char* hgamma_em = "Hgamma";
    const char* hdelta_em = "Hdelta";
    const char* h8_em = "H8";
    const char* h9_em = "H9";
    const char* h10_em = "H10";
    const char* h11_em = "H11";


    const char* oII3726_em = "[OII]3726";
    const char* oII3729_em = "[OII]3729";

    const char* oIIIa_em = "[OIII](doublet-1)";
    const char* oIIIb_em = "[OIII](doublet-1/3)";

    const char* cIII1907_em = "[CIII]1907";
    const char* cIII1909_em = "[CIII]1909";

    const char* lya_em = "LyAE";

    const char* niia_em = "[NII](doublet-1)";
    const char* niib_em = "[NII](doublet-1/2.95)";


};

}

#endif

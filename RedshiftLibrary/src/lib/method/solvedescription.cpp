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
#include "RedshiftLibrary/method/solvedescription.h"
using namespace NSEpic;

/**
 * \brief Returns a string describing the names and allowed values for the parameters of the Linemodel method.
 **/
const std::string CSolveDescription::GetDescription(const std::string& method)
{
    std::string desc;

    if (method == "linemodelsolve")
      {
        desc = "Method linemodel:\n";
        
        desc.append("\tparam: linemodel.linetypefilter = {""no"", ""E"", ""A""}\n");
        desc.append("\tparam: linemodel.lineforcefilter = {""no"", ""S""}\n");
        desc.append("\tparam: linemodel.fittingmethod = {""hybrid"", ""individual""}\n");
        desc.append("\tparam: linemodel.continuumcomponent = {""fromspectrum"", ""tplfit"", ""nocontinuum"", ""zero""}\n");
        desc.append("\tparam: linemodel.continuumfit.method = {""templatefitting"", ""templatefittinglog""}\n");
        desc.append("\tparam: linemodel.continuumfit.ismfit = {""false"", ""true""}\n");
        desc.append("\tparam: linemodel.continuumfit.igmfit = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.continuumfit.count = <float value>\n");
        desc.append("\tparam: linemodel.continuumfit.ignorelinesupport = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.continuumfit.priors.betaA = <float value>\n");
        desc.append("\tparam: linemodel.continuumfit.priors.betaTE = <float value>\n");
        desc.append("\tparam: linemodel.continuumfit.priors.betaZ = <float value>\n");
        desc.append("\tparam: linemodel.continuumfit.priors.catalog_dirpath = <path>\n");
        desc.append("\tparam: linemodel.secondpasslcfittingmethod = {""-1"", ""svdlc"",""svdlcp2""}\n");
        desc.append("\tparam: linemodel.rigidity = {""rules"", ""tplcorr"", ""tplshape""}\n");
        desc.append("\tparam: linemodel.tplratio_catalog = <relative path>\n");
        desc.append("\tparam: linemodel.tplratio_ismfit = {""false"",""true""}\n");

        desc.append("\tparam: linemodel.tplratio.priors.betaA = <float value>\n");
        desc.append("\tparam: linemodel.tplratio.priors.betaTE = <float value>\n");
        desc.append("\tparam: linemodel.tplratio.priors.betaZ = <float value>\n");
        desc.append("\tparam: linemodel.tplratio.priors.catalog_dirpath = <path>\n");


        desc.append("\tparam: linemodel.enableLSF = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""velocitydriven"",  ""combined"",  ""nispvsspsf201707"", ""fixed""}\n");
        desc.append("\tparam: linemodel.nsigmasupport = <float value>\n");
        desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
        desc.append("\tparam: linemodel.velocityemission = <float value>\n");
        desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
        desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
        desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""superstrong"", ""ratiorange"", ""ciiiratio"", ""no""}\n");
        desc.append("\tparam: linemodel.improveBalmerFit = {""false"", ""true""}\n");
        desc.append("\tparam: linemodel.extremacount = <float value>\n");
        desc.append("\tparam: linemodel.extremacountB = <float value>\n");
        desc.append("\tparam: linemodel.extremacutprobathreshold = <float value> (-1:disabled)\n");
        desc.append("\tparam: linemodel.velocityfit = {""false"", ""true""}\n");
        desc.append("\tparam: linemodel.emvelocityfitmin = <float value>\n");
        desc.append("\tparam: linemodel.emvelocityfitmax = <float value>\n");
        desc.append("\tparam: linemodel.emvelocityfitstep = <float value>\n");
        desc.append("\tparam: linemodel.absvelocityfitmin = <float value>\n");
        desc.append("\tparam: linemodel.absvelocityfitmax = <float value>\n");
        desc.append("\tparam: linemodel.absvelocityfitstep = <float value>\n");

        desc.append("\tparam: linemodel.lyaforcefit = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.lyaforcedisablefit = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.lyafit.asymfitmin = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.asymfitmax = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.asymfitstep = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.widthcoefffitmin = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.widthcoefffitmax = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.widthcoefffitstep = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.deltafitmin = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.deltafitmax = <float value>\n");
        desc.append("\tparam: linemodel.lyafit.deltafitstep = <float value>\n");

        //first pass
        desc.append("\tparam: linemodel.firstpass.largegridstep = <float value>, deactivated if negative or zero\n");
        desc.append("\tparam: linemodel.firstpass.tplratio_ismfit = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.firstpass.multiplecontinuumfit_disable = {""false"",""true""}\n");

        //second pass
        desc.append("\tparam: linemodel.skipsecondpass = {""false"",""true""}\n");
        desc.append("\tparam: linemodel.secondpass.continuumfit = {""fromfirstpass"", ""retryall"", ""refitfirstpass""}\n");
        desc.append("\tparam: linemodel.secondpass.halfwindowsize = <float value>\n");

        desc.append("\tparam: linemodel.pdfcombination = {""marg"", ""bestchi2""}\n");
        desc.append("\tparam: linemodel.pdf.margampcorr = {""false"", ""true""}\n");

        desc.append("\tparam: linemodel.stronglinesprior = <float value>, penalization factor = positive value or -1 to deactivate\n");
        desc.append("\tparam: linemodel.haprior = <float value>, penalization factor = positive value (typical 1e-1 to 1e-5) or -1 to deactivate\n");
        desc.append("\tparam: linemodel.euclidnhaemittersStrength = <float value>, prior strength factor = positive value (typically 1 to 5) or -1 to deactivate\n");
        desc.append("\tparam: linemodel.modelpriorzStrength = <float value>, prior strength factor = positive value (typically 1 to 5) or -1 to deactivate\n");

        desc.append("\tparam: linemodel.saveintermediateresults = {""false"", ""true""}\n");
      }
    else if(method == "templatefittingsolve")
      {
        desc = "Method templatefittingsolve:\n";

        desc.append("\tparam: templatefittingsolve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
        desc.append("\tparam: templatefittingsolve.overlapThreshold = <float value>\n");
        desc.append("\tparam: templatefittingsolve.interpolation = {""precomputedfinegrid"", ""lin""}\n");
        desc.append("\tparam: templatefittingsolve.extinction = {""false"", ""true""}\n");
        desc.append("\tparam: templatefittingsolve.dustfit = {""false"", ""true""}\n");
        desc.append("\tparam: templatefittingsolve.pdfcombination = {""marg"", ""bestchi2""}\n");
        desc.append("\tparam: templatefittingsolve.saveintermediateresults = {""false"", ""true""}\n");
      }
    else if(method=="tplcombinationsolve")
      {
        desc = "Method tplcombinationsolve:\n";

        desc.append("\tparam: tplcombinationsolve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
        desc.append("\tparam: tplcombinationsolve.overlapThreshold = <float value>\n");
        desc.append("\tparam: tplcombinationsolve.interpolation = {""precomputedfinegrid"", ""lin""}\n");
        //desc.append("\tparam: tplcombinationsolve.extinction = {""false"", ""true""}\n");
        //desc.append("\tparam: tplcombinationsolve.dustfit = {""false"", ""true""}\n");
        //desc.append("\tparam: tplcombinationsolve.pdfcombination = {""marg"", ""bestchi2""}\n");
        desc.append("\tparam: tplcombinationsolve.saveintermediateresults = {""false"", ""true""}\n");
      }
    return desc;

}

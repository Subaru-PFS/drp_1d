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
#include <unsupported/Eigen/NumericalDiff>

#include "RedshiftLibrary/linemodel/gaussianfit/lbfgsbfitter.h"

using namespace std;
using namespace Eigen;
using namespace LBFGSpp;
using namespace NSEpic;

CLbfgsbFitter::CLeastSquare::CLeastSquare(
    CLbfgsbFitter &fitter, const TInt32List &EltsIdx, CLine::EType lineType,
    Float64 redshift, const TInt32List &xInds, Int32 velA_idx, Int32 velE_idx,
    Int32 lbdaOffset_idx, Int32 pCoeff_idx, Float64 normFactor, Float64 normVel,
    Float64 normLbdaOffset)
    : Functor<Float64, Dynamic, 1>(), m_fitter(&fitter), m_EltsIdx(&EltsIdx),
      m_lineType(lineType), m_redshift(redshift), m_xInds(&xInds),
      m_velA_idx(velA_idx), m_velE_idx(velE_idx),
      m_lbdaOffset_idx(lbdaOffset_idx), m_pCoeff_idx(pCoeff_idx),
      m_normFactor(normFactor), m_normVel(normVel),
      m_normLbdaOffset(normLbdaOffset),
      m_spectralAxis(&fitter.getSpectrum().GetSpectralAxis()),
      m_noContinuumFluxAxis(&fitter.getModel().getSpcFluxAxisNoContinuum()),
      m_continuumFluxAxis(&fitter.getModel().getContinuumFluxAxis()),
      m_ErrorNoContinuum(&fitter.getSpectrum().GetErrorAxis()) {

  // init normalized polynomial
  if (m_fitter->m_enableAmplitudeOffsets) {
    Float64 x_center = 0.5 * ((*m_spectralAxis)[m_xInds->back()] +
                              (*m_spectralAxis)[m_xInds->front()]);
    Float64 scale = 0.5 * ((*m_spectralAxis)[m_xInds->back()] -
                           (*m_spectralAxis)[m_xInds->front()]);
    m_pCoeffs = CPolynomCoeffsNormalized(x_center, scale);
  }

  // precompute data sum square
  Float64 sumSquare = 0.;
  for (Int32 i = 0; i < m_xInds->size(); i++) {
    Float64 xi, yi, ei, ei2;
    Int32 idx = (*m_xInds)[i];
    yi = (*m_noContinuumFluxAxis)[idx] * m_normFactor;
    ei = (*m_ErrorNoContinuum)[idx] * m_normFactor;
    ei2 = ei * ei;
    sumSquare += yi * yi / ei2;

    m_sumSquareData = sumSquare;
  }
}

// compute Sum Squared diff
void CLbfgsbFitter::CLeastSquare::operator()(const VectorXd &x,
                                             ValueType &val) const {

  CPolynomCoeffsNormalized pCoeffs = unpack(x);

  m_fitter->m_spectraIndex
      .reset(); // temporary multiobs implementation (might be
  //  useless -> investigate)
  val(0, 0) = ComputeLeastSquare(pCoeffs);

  Log.LogDebug(Formatter() << "LeastSquare = " << val(0, 0));
}

/*
// compute Sum Squared diff and gradient
CLbfgsbFitter::CLeastSquare::ValueType
CLbfgsbFitter::CLeastSquare::operator()(const VectorXd &x, VectorXd &grad) {

  TFloat64List amps;
  Float64 redshift;
  TFloat64List pCoeffs;
  std::tie(amps, redshift, pCoeffs) = unpack(x);

  Float64 lsq = ComputeLeastSquareAndGrad(amps, redshift, pCoeffs, grad);

  Log.LogDebug(Formatter() << "LeastSquare = " << lsq);
  for (size_t i = 0; i < m_EltsIdx->size(); ++i)
    Log.LogDebug(Formatter() << "LeastSquare grad amp[" << i <<"] = " <<
grad[i]); if (m_lineType == CLine::EType::nType_All) { Log.LogDebug(Formatter()
<< "LeastSquare grad velA= " << grad[m_velA_idx]); Log.LogDebug(Formatter() <<
"LeastSquare grad velE= " << grad[m_velE_idx]); } else if (m_lineType ==
CLine::EType::nType_Absorption) Log.LogDebug(Formatter() << "LeastSquare grad
velA= "
<< grad[m_velA_idx]); else Log.LogDebug(Formatter() << "LeastSquare grad velE= "
<< grad[m_velE_idx]); Log.LogDebug(Formatter() << "LeastSquare grad zoffset= "
<< grad[m_lbdaOffset_idx]); Log.LogDebug(Formatter() << "LeastSquare grad p
coeffs= "
                           << grad[m_lbdaOffset_idx + 1] << " "
                           << grad[m_lbdaOffset_idx + 2] << " "
                           << grad[m_lbdaOffset_idx + 3]);
  ValueType ret{{lsq}};
  return ret;
}
*/

CPolynomCoeffsNormalized
CLbfgsbFitter::CLeastSquare::unpack(const VectorXd &x) const {

  // unpack param vector
  TFloat64List amps(x.begin(), x.begin() + m_EltsIdx->size());
  for (Int32 i = 0; i < m_EltsIdx->size(); ++i) {
    Log.LogDebug(Formatter() << "amplitude[" << i << "]= " << amps[i]);
    m_fitter->m_ElementsVector->SetElementAmplitude((*m_EltsIdx)[i], amps[i],
                                                    0.0);
  }

  if (m_fitter->m_enableVelocityFitting) {

    if (m_velA_idx != undefIdx)
      Log.LogDebug(Formatter()
                   << "velocity Abs = " << x[m_velA_idx] * m_normVel);
    if (m_velE_idx != undefIdx)
      Log.LogDebug(Formatter()
                   << "velocity Em  = " << x[m_velE_idx] * m_normVel);

    for (Int32 eltIndex : *m_EltsIdx) {
      auto &elt_param = m_fitter->getElementParam()[eltIndex];

      if (elt_param->IsEmission())
        elt_param->SetVelocityEmission(x[m_velE_idx] * m_normVel);
      else
        elt_param->SetVelocityAbsorption(x[m_velA_idx] * m_normVel);
    }
  }

  if (m_fitter->m_enableLambdaOffsetsFit) {
    const Float64 delta_offset = x[m_lbdaOffset_idx] * m_normLbdaOffset;
    Log.LogDebug(Formatter() << "lambda offset = " << delta_offset);
    /* Float64 const redshift = m_redshift + (1.0 + m_redshift) * delta_offset;
    Log.LogDebug(Formatter() << "redshift=" << redshift); */
    for (Int32 eltIndex : *m_EltsIdx) {
      auto &elt_param = m_fitter->getElementParam()[eltIndex];

      for (Int32 line_idx = 0; line_idx != elt_param->size(); ++line_idx) {
        Float64 offset = elt_param->m_Lines[line_idx].GetOffset();
        offset /= SPEED_OF_LIGHT_IN_VACCUM;
        offset += (1 + offset) * delta_offset;
        offset *= SPEED_OF_LIGHT_IN_VACCUM;
        elt_param->setLambdaOffset(line_idx, offset);

        /* // set redshift in SYMIGM profiles
        const TInt32List &igm_indices = elt_param->m_asymLineIndices;
        for (Int32 idx : igm_indices) {
          auto igmp = elt_param->GetSymIgmParams();
          igmp.m_redshift = redshift;
          elt_param->SetSymIgmParams(igmp, idx);
        } */
      }
    }
  }

  // reset the support (lines outside range)
  m_fitter->m_spectraIndex.reset();
  for (Int32 eltIndex : *m_EltsIdx) {
    auto &elt_ptr = m_fitter->getElementList()[eltIndex];
    elt_ptr->prepareSupport(*m_spectralAxis, m_redshift,
                            m_fitter->getLambdaRange());
  }
  m_fitter->m_ElementsVector->computeGlobalOutsideLambdaRange();

  CPolynomCoeffsNormalized pCoeffs = m_pCoeffs;
  if (m_fitter->m_enableAmplitudeOffsets) {
    pCoeffs.a0 = x[m_pCoeff_idx];
    pCoeffs.a1 = x[m_pCoeff_idx + 1];
    pCoeffs.a2 = x[m_pCoeff_idx + 2];
    Log.LogDebug(Formatter() << "p coeffs  = " << pCoeffs.a0 << " "
                             << pCoeffs.a1 << " " << pCoeffs.a2);
  }

  return pCoeffs;
}

Float64 CLbfgsbFitter::CLeastSquare::ComputeLeastSquare(
    const CPolynomCoeffsNormalized &pCoeffs) const {
  // compute least square term
  Float64 sumSquare = m_sumSquareData;
  for (Int32 i = 0; i < m_xInds->size(); i++) {
    Float64 xi, yi, ei, ei2;
    Int32 idx = (*m_xInds)[i];
    xi = (*m_spectralAxis)[idx];
    yi = (*m_noContinuumFluxAxis)[idx] * m_normFactor;
    ei = (*m_ErrorNoContinuum)[idx] * m_normFactor;
    ei2 = ei * ei;

    // compute model value
    Float64 fval = 0.;
    for (auto &eltIndex : *m_EltsIdx) {
      auto &elt = m_fitter->getElementList()[eltIndex];
      // linemodel value
      Float64 mval =
          elt->getModelAtLambda(xi, m_redshift, (*m_continuumFluxAxis)[idx]);
      fval += mval;
    }

    if (m_fitter->m_enableAmplitudeOffsets)
      fval += pCoeffs.getValue(xi);

    // add squared diff
    sumSquare += (fval * fval - 2.0 * yi * fval) / ei2;
  }

  return sumSquare;
}

Float64 CLbfgsbFitter::CLeastSquare::ComputeLeastSquareAndGrad(
    const TFloat64List &amps, Float64 redshift,
    const CPolynomCoeffsNormalized &pCoeffs, VectorXd &grad) const {

  // compute least square term
  Float64 sumSquare = m_sumSquareData;
  for (Int32 i = 0; i < m_xInds->size(); i++) {
    Float64 xi, yi, ei, ei2;
    Int32 idx = (*m_xInds)[i];
    xi = (*m_spectralAxis)[idx];
    yi = (*m_noContinuumFluxAxis)[idx] * m_normFactor;
    ei = (*m_ErrorNoContinuum)[idx] * m_normFactor;
    ei2 = ei * ei;

    // compute model value and gradient
    Float64 fval = 0.;
    TFloat64List ampsGrad;
    Float64 velocityAGrad = 0.0;
    Float64 velocityEGrad = 0.0;
    Float64 offsetGrad = 0.;
    TFloat64List pCoeffGrad;
    auto amps_it = amps.cbegin();
    for (auto &eltIndex : *m_EltsIdx) {
      auto &elt = m_fitter->getElementList()[eltIndex];
      const auto &elt_param = elt->getElementParam();
      // linemodel value &  amplitude derivative
      Float64 mval =
          elt->getModelAtLambda(xi, redshift, (*m_continuumFluxAxis)[idx]);
      ampsGrad.push_back(0.0);
      if (*amps_it != 0.0)
        ampsGrad.back() = mval / *amps_it;
      ++amps_it;
      fval += mval;

      // velocity derivative
      if (m_fitter->m_enableVelocityFitting) {
        if (elt_param->GetElementType() == CLine::EType::nType_Absorption) {
          velocityAGrad += elt->GetModelDerivVelAtLambda(
              xi, redshift, (*m_continuumFluxAxis)[idx]);
        } else {
          velocityEGrad += elt->GetModelDerivVelAtLambda(
              xi, redshift, (*m_continuumFluxAxis)[idx]);
        }
      }

      // offset derivative
      if (m_fitter->m_enableLambdaOffsetsFit)
        offsetGrad += elt->GetModelDerivZAtLambda(xi, redshift,
                                                  (*m_continuumFluxAxis)[idx]);
    }

    if (m_fitter->m_enableAmplitudeOffsets)
      fval += pCoeffs.getValueAndGrad(xi, pCoeffGrad);

    // add squared diff
    sumSquare += (fval * fval - 2.0 * yi * fval) / ei2;
    Float64 residual = -2.0 * (yi - fval) / ei2;
    for (Int32 eltIndex = 0; eltIndex < m_EltsIdx->size(); ++eltIndex) {
      // squared diff derivative wrt amplitudes
      grad[eltIndex] += residual * ampsGrad[eltIndex];

      // squared diff derivative wrt velocity
      if (m_fitter->m_enableVelocityFitting) {
        if (m_fitter->getElementParam()[eltIndex]->GetElementType() ==
            CLine::EType::nType_Absorption) {
          grad[m_velA_idx] += residual * velocityAGrad;
        } else {
          grad[m_velE_idx] += residual * velocityEGrad;
        }
      }

      // squared diff derivative wrt offset
      grad[m_lbdaOffset_idx] += residual * offsetGrad;
    }

    // squared diff derivative wrt polynome coeffs
    if (m_fitter->m_enableAmplitudeOffsets) {
      for (size_t i = 0; i <= TPolynomCoeffs::degree; ++i)
        grad[m_pCoeff_idx + i] += residual * pCoeffGrad[i];
    }
  }

  return sumSquare;
}

void CLbfgsbFitter::resetSupport(Float64 redshift) {

  // set velocity at max value (to set largest line overlapping)
  if (Context.GetParameterStore()->GetScoped<bool>("lineModel.velocityFit")) {
    for (Int32 j = 0; j < getElementList().size(); j++) {
      getElementParam()[j]->m_VelocityEmission = m_velfitMaxE;
      getElementParam()[j]->m_VelocityAbsorption = m_velfitMaxA;
      m_enlarge_line_supports = MAX_LAMBDA_OFFSET;
    }
  }

  CAbstractFitter::resetSupport(redshift);
}

void CLbfgsbFitter::doFit(Float64 redshift) {

  // fit the amplitudes of each element independently, unless there is overlap
  CHybridFitter::fitAmplitudesHybrid(redshift);
}

// since fitAmplitudesLinSolve is already fitting lambda offsets,
// this is a simple wrapper
void CLbfgsbFitter::fitAmplitudesLinSolveAndLambdaOffset(
    TInt32List EltsIdx, bool enableOffsetFitting, Float64 redshift) {

  fitAmplitudesLinSolvePositive(EltsIdx, redshift);
}

// overriding the SVD linear fitting, but here it is not linear inversion
void CLbfgsbFitter::fitAmplitudesLinSolvePositive(const TInt32List &EltsIdx,
                                                  Float64 redshift) {

  if (EltsIdx.size() < 1)
    THROWG(ErrorCode::INTERNAL_ERROR, "empty Line element list to fit");

  Int32 nddl = EltsIdx.size();
  Int32 param_idx = nddl;
  Int32 velA_idx = undefIdx;
  Int32 velE_idx = undefIdx;
  CLine::EType lineType = CLine::EType::nType_All;
  if (m_enableVelocityFitting) {
    ++nddl;
    // position of velocity in the param vector
    Int32 velocity_param_idx = param_idx++;

    // check if mixed types (Abs & Em)
    lineType = getElementParam()[EltsIdx.front()]->GetElementType();
    if (lineType == CLine::EType::nType_Absorption)
      velA_idx = velocity_param_idx;
    else
      velE_idx = velocity_param_idx;
    for (Int32 eltIndex : EltsIdx)
      if (lineType != getElementParam()[eltIndex]->GetElementType()) {
        lineType = CLine::EType::nType_All;
        ++nddl; // 2 velocity parameter
        velA_idx = velocity_param_idx;
        velE_idx = velocity_param_idx + 1;
        ++param_idx;
      }
  }

  Int32 lbdaOffset_param_idx =
      undefIdx; // position of lambda offset in the param vector
  if (m_enableLambdaOffsetsFit) {
    ++nddl;
    lbdaOffset_param_idx = param_idx++;
  }

  TInt32List xInds = getElementList().getSupportIndexes(EltsIdx);
  if (xInds.size() < 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "no observed samples for the line Element to fit");

  Int32 pCoeff_param_idx =
      undefIdx; // position of polynomial coeffs in the param vector
  if (m_enableAmplitudeOffsets) {
    nddl += TPolynomCoeffs::degree + 1;
    pCoeff_param_idx = param_idx;
    param_idx += TPolynomCoeffs::degree + 1;
  }

  // check N DOF
  Int32 n = xInds.size();
  if (n < nddl) {
    Flag.warning(WarningCode::LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT,
                 Formatter() << __func__ << " LBFGSB ill ranked:"
                             << " number of samples = " << n
                             << ", number of parameters to fit = " << nddl);
    for (Int32 eltIndex : EltsIdx)
      m_ElementsVector->SetElementAmplitude(eltIndex, 0., INFINITY);
    if (m_enableAmplitudeOffsets) {
      for (Int32 eltIndex : EltsIdx)
        m_ElementsVector->getElementParam()[eltIndex]->SetPolynomCoeffs(
            {0., 0., 0.});
    }
    return;
  }

  // Normalize
  // Float64 maxabsval = DBL_MIN;
  const auto &noContinuumFluxAxis = getModel().getSpcFluxAxisNoContinuum();
  Int32 maxabsval_idx =
      *std::max_element(xInds.cbegin(), xInds.cend(), [&](Int32 l, Int32 r) {
        return std::abs(noContinuumFluxAxis[l]) <
               std::abs(noContinuumFluxAxis[r]);
      });

  Float64 normFactor = 1.0 / std::abs(noContinuumFluxAxis[maxabsval_idx]);
  Float64 normVel = m_velfitMaxE;
  Float64 normLbdaOffset = m_LambdaOffsetMax / SPEED_OF_LIGHT_IN_VACCUM;

  CLeastSquare func(*this, EltsIdx, lineType, redshift, xInds, velA_idx,
                    velE_idx, lbdaOffset_param_idx, pCoeff_param_idx,
                    normFactor, normVel, normLbdaOffset);

  // bounds for all params
  /////////////////////////
  VectorXd lb(nddl);
  VectorXd ub(nddl);

  // amplitude bounds
  const auto amp_indices = seq(0, EltsIdx.size() - 1);
  lb(amp_indices) = VectorXd::Zero(EltsIdx.size());
  for (size_t i = 0; i < EltsIdx.size(); ++i) {
    Float64 ampMax = INFINITY;
    auto &elt_param = getElementParam()[EltsIdx[i]];
    if (elt_param->GetElementType() == CLine::EType::nType_Absorption &&
        elt_param->GetAbsLinesLimit() > 0.0)
      ampMax =
          elt_param->GetAbsLinesLimit() / elt_param->GetMaxNominalAmplitude();
    ub[i] = ampMax;
  }

  // velocity bounds
  if (m_enableVelocityFitting) {
    if (lineType == CLine::EType::nType_All) {
      const auto vel_indices = seq(velA_idx, velE_idx);
      lb(vel_indices) << m_velfitMinA / normVel, m_velfitMinE / normVel;
      ub(vel_indices) << m_velfitMaxA / normVel, m_velfitMaxE / normVel;
    } else if (lineType == CLine::EType::nType_Absorption) {
      lb[velA_idx] = m_velfitMinA / normVel;
      ub[velA_idx] = m_velfitMaxA / normVel;
    } else {
      lb[velE_idx] = m_velfitMinE / normVel;
      ub[velE_idx] = m_velfitMaxE / normVel;
    }
  }

  // TODO adapt bounds to lambdaRange (if line goes out of border #8669)
  if (m_enableLambdaOffsetsFit) {
    // offset bounds
    lb[lbdaOffset_param_idx] =
        m_LambdaOffsetMin / SPEED_OF_LIGHT_IN_VACCUM / normLbdaOffset;
    ub[lbdaOffset_param_idx] =
        m_LambdaOffsetMax / SPEED_OF_LIGHT_IN_VACCUM / normLbdaOffset;
  }

  // polynomial bounds
  if (m_enableAmplitudeOffsets) {
    constexpr size_t ncoeff = TPolynomCoeffs::degree + 1;
    auto pCoeff_indices = seq(pCoeff_param_idx, pCoeff_param_idx + ncoeff - 1);
    auto pCoeffMin = lb(pCoeff_indices);
    auto pCoeffMax = ub(pCoeff_indices);
    for (size_t i = 0; i < ncoeff; ++i) {
      pCoeffMin[i] = -INFINITY;
      pCoeffMax[i] = INFINITY;
    }
  }

  // initial guess for all params
  ///////////////////////////////
  VectorXd v_xGuess(nddl);

  // compute amplitudes initial guess using CSvdFitter
  m_spectraIndex.reset();
  for (auto eltIndex : EltsIdx) {
    auto &elt_param = getElementParam()[eltIndex];
    // set velocity guess
    if (elt_param->GetElementType() == CLine::EType::nType_Absorption) {
      elt_param->SetVelocityAbsorption(m_velIniGuessA);
    } else {
      elt_param->SetVelocityEmission(m_velIniGuessE);
    }

    // reset support with initial velocities
    auto &elt_ptr = getElementList()[eltIndex];
    elt_ptr->prepareSupport(getSpectrum().GetSpectralAxis(), redshift,
                            getLambdaRange());
  }
  m_ElementsVector->computeGlobalOutsideLambdaRange();
  m_spectraIndex.reset(); // dummy reset, TO BE CORRECTED for full multiobs
  CSvdFitter::fitAmplitudesLinSolvePositive(EltsIdx, redshift);
  Float64 max_snr = -INFINITY;
  for (size_t i = 0; i != EltsIdx.size(); ++i) {
    auto &elt_param = getElementParam()[EltsIdx[i]];
    if (elt_param->isNotFittable()) {
      // the initial velocity renders the line outside range, set amplitude at
      // zero and do not update the max snr
      v_xGuess[i] = 0.0;
      continue;
    }
    v_xGuess[i] = elt_param->GetElementAmplitude() * normFactor;
    if (std::isnan(v_xGuess[i]))
      THROWG(ErrorCode::INTERNAL_ERROR,
             "NAN amplitude for LBFGSB fitter initial guess");
    // retrive max SNR amplitude:
    auto sigma = (elt_param->GetElementAmplitudeError() * normFactor);
    auto snr = v_xGuess[i] / sigma;
    max_snr = std::max(max_snr, snr);
  }

  // if all amplitudes SNR are too small, keep the guess values
  // since the precise fit may fail because of too noisy least-square
  Float64 min_snr_threshold = 1.0;
  if (max_snr < min_snr_threshold)
    return;

  // polynomial coeffs initial guess
  if (m_enableAmplitudeOffsets) {
    const auto &pCoeffs =
        getElementParam()[EltsIdx.front()]->GetPolynomCoeffs();
    auto pCoeffsNormalized = func.getPcoeffs();
    pCoeffsNormalized.setCoeffs(pCoeffs.a0 * normFactor,
                                pCoeffs.a1 * normFactor,
                                pCoeffs.a2 * normFactor);
    v_xGuess[pCoeff_param_idx] = pCoeffsNormalized.a0;
    v_xGuess[pCoeff_param_idx + 1] = pCoeffsNormalized.a1;
    v_xGuess[pCoeff_param_idx + 2] = pCoeffsNormalized.a2;
  }

  // velocity initial guess
  if (m_enableVelocityFitting) {
    if (lineType == CLine::EType::nType_All) {
      const auto vel_indices = seq(velA_idx, velE_idx);
      v_xGuess(vel_indices) << m_velIniGuessA / normVel,
          m_velIniGuessE / normVel;
    } else if (lineType == CLine::EType::nType_Absorption)
      v_xGuess[velA_idx] = m_velIniGuessA / normVel;
    else
      v_xGuess[velE_idx] = m_velIniGuessE / normVel;
  }

  // lambda offset inital guess
  if (m_enableLambdaOffsetsFit)
    v_xGuess[lbdaOffset_param_idx] = 0.;

  // Set up solver parameters
  LBFGSBParam<Float64> param; // New parameter class
  /*  param.epsilon = 1; //large to use only rel tolerance
    //param.epsilon_rel = 1e-5;
    param.ftol = 1e-2;
    param.max_iterations = 100;
    // param.max_linesearch = 30;
  */
  // mimic scipy params
  // param.epsilon = 1e-5;
  // param.epsilon_rel = 0.0;
  param.epsilon = 0.0;
  param.epsilon_rel = 1e-5;
  param.delta = 2.22045e-9;
  param.past = 1;

  // param.m = 10;
  param.max_iterations = 15000;

  param.max_linesearch = 20;
  // Less sure about these ones:
  param.ftol = 1e-3;
  param.wolfe = 0.9;
  param.min_step = 0.0;
  param.max_step = std::numeric_limits<double>::infinity();
  // Edit: update this to zero following comments below.
  param.max_submin = 0;

  // Create solver and function object
  LBFGSBSolver<Float64, LineSearchMoreThuente> solver(param);

  NumericalDiff<CLeastSquare> func_with_num_diff(func, 1e-12);

  auto myfunc = [&func_with_num_diff](const VectorXd &x, VectorXd &grad) {
    CLeastSquare::JacobianType J(1, x.size());
    func_with_num_diff.df(x, J);
    grad = J;
    CLeastSquare::ValueType val;
    func_with_num_diff(x, val);
    return val(0, 0);
  };

  Float64 fx;
  VectorXd v_xResult = v_xGuess;
  try {
    // v_xResult will be overwritten to be the best point found
    int niter = solver.minimize(myfunc, v_xResult, fx, lb, ub);
  } catch (const AmzException &e) {
    // throw again, since the exception was thrown inside amazed
    // CLbfgsbFitter::CLeastSquare
    throw;
  } catch (const std::exception &e) {
    // the error was raised from lbfgsbpp
    Flag.warning(WarningCode::LBFGSPP_ERROR, Formatter()
                                                 << "in LBFGSPP, " << e.what());
    // reset the result to the initial guess
    v_xResult = v_xGuess;
  }

  // store fitted amplitude, with temporary unknown SNR
  Float64 amp_std = NAN;
  for (Int32 i = 0; i < EltsIdx.size(); ++i)
    m_ElementsVector->SetElementAmplitude(EltsIdx[i], v_xResult[i] / normFactor,
                                          amp_std);

  // store fitted velocity dispersion (line width)
  if (m_enableVelocityFitting) {
    if (lineType == CLine::EType::nType_All) {
      Float64 velocityA = v_xResult[velA_idx] * normVel;
      Float64 velocityE = v_xResult[velE_idx] * normVel;
      for (Int32 eltIndex : EltsIdx) {
        auto &elt_param = getElementParam()[eltIndex];
        if (elt_param->GetElementType() == CLine::EType::nType_Absorption)
          elt_param->SetVelocityAbsorption(velocityA);
        else
          elt_param->SetVelocityEmission(velocityE);
      }
    } else if (lineType == CLine::EType::nType_Absorption) {
      Float64 velocity = v_xResult[velA_idx] * normVel;
      for (Int32 eltIndex : EltsIdx)
        getElementParam()[eltIndex]->SetVelocityAbsorption(velocity);
    } else {
      Float64 velocity = v_xResult[velE_idx] * normVel;
      for (Int32 eltIndex : EltsIdx)
        getElementParam()[eltIndex]->SetVelocityEmission(velocity);
    }
  }

  // store velocity offset (line offset)
  if (m_enableLambdaOffsetsFit) {
    for (Int32 eltIndex : EltsIdx) {
      auto &elt_param = m_ElementsVector->getElementParam()[eltIndex];
      for (Int32 line_idx = 0; line_idx != elt_param->size(); ++line_idx) {
        Float64 offset = elt_param->m_Lines[line_idx].GetOffset();
        offset /= SPEED_OF_LIGHT_IN_VACCUM;
        offset +=
            (1 + offset) * v_xResult[lbdaOffset_param_idx] * normLbdaOffset;
        offset *= SPEED_OF_LIGHT_IN_VACCUM;
        elt_param->setLambdaOffset(line_idx, offset);
      }
    }
  }

  // store polynomial coeffs (continuum under line)
  if (m_enableAmplitudeOffsets) {
    auto pCoeffsNormalized = func.getPcoeffs();
    pCoeffsNormalized.a0 = v_xResult[pCoeff_param_idx];
    pCoeffsNormalized.a1 = 0.0;
    pCoeffsNormalized.a2 = 0.0;
    if (TPolynomCoeffs::degree > 0)
      pCoeffsNormalized.a1 = v_xResult[pCoeff_param_idx + 1];
    if (TPolynomCoeffs::degree > 1)
      pCoeffsNormalized.a2 = v_xResult[pCoeff_param_idx + 2];
    Float64 a0 = 0.0;
    Float64 a1 = 0.0;
    Float64 a2 = 0.0;
    pCoeffsNormalized.getCoeffs(a0, a1, a2);
    for (Int32 eltIndex : EltsIdx)
      m_ElementsVector->getElementParam()[eltIndex]->SetPolynomCoeffs(
          {a0 / normFactor, a1 / normFactor, a2 / normFactor});
  }

  // compute snr (analytical formula from doi:10.1364/AO.46.005374)
  // --------------------------------------------------------------
  //  Var(amp) = 3/2 * var(spectra) /
  //       (sqrt(pi) * sigma_in_wavelength_unit * pixsize_in_wavelength_unit)
  //
  // note: this is only for one line (one Gaussian profile)
  // Thus it does not take into account the correlation between overlapping
  // lines, and possibly underestimate the uncertainty in such a case
  // reset the support (lines outside range)
  m_spectraIndex.reset();
  for (Int32 eltIndex : EltsIdx) {
    auto &elt_ptr = getElementList()[eltIndex];
    elt_ptr->prepareSupport(getSpectrum().GetSpectralAxis(), redshift,
                            getLambdaRange());
  }
  m_ElementsVector->computeGlobalOutsideLambdaRange();
  m_spectraIndex.reset();
  getModel().refreshModel(); // recompute model with estimated params (needed
                             // for noise estimation from residual)

  for (Int32 i = 0; i < EltsIdx.size(); ++i) {
    m_spectraIndex.reset();
    auto const &elt = getElementList()[EltsIdx[i]];
    Float64 const amp = v_xResult[i] / normFactor;

    if (elt->IsOutsideLambdaRange()) {
      m_ElementsVector->SetElementAmplitude(EltsIdx[i], NAN, NAN);
      break;
    }

    auto const &[_, sigma] = elt->getObservedPositionAndLineWidth(redshift);
    Float64 const noise_stdev =
        getModelResidualRmsUnderElements({EltsIdx[i]}, true, false);
    Float64 pixsize = 0.0;
    Int32 nlines = 0;
    for (Int32 line_idx = 0; line_idx < elt->GetSize(); ++line_idx) {
      if (elt->IsOutsideLambdaRangeLine(line_idx))
        continue;
      auto const range = elt->getSupportSubElt(line_idx);
      m_spectraIndex.reset(); // dummy reset, TO BE CORRECTED for full multiobs
      pixsize += getSpectrum().GetSpectralAxis().GetMeanResolution(range);
      ++nlines;
    }
    pixsize /= nlines;
    Float64 const amp_std =
        noise_stdev * sqrt(3.0 / (2.0 * sigma * pixsize * sqrt(M_PI)));

    m_ElementsVector->SetElementAmplitude(EltsIdx[i], amp, amp_std);

    // we should divide the amplitude uncertainty by the nominal amp,  since it
    // has been multplied by it in SetElementAmplitude....
    auto const &param = m_ElementsVector->getElementParam()[EltsIdx[i]];
    for (Int32 line_idx = 0; line_idx < elt->GetSize(); ++line_idx)
      param->m_FittedAmplitudesStd[line_idx] /=
          param->m_NominalAmplitudes[line_idx];
  }
}

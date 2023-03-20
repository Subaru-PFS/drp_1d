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
#include "RedshiftLibrary/linemodel/gaussianfit/lbfgsfitter.h"

using namespace std;
using Eigen::VectorXd;
using namespace LBFGSpp;
using namespace NSEpic;

CLbfgsFitter::CLeastSquare::CLeastSquare(
    CLbfgsFitter &fitter, const TInt32List &EltsIdx, Int32 lineType,
    Float64 redshift, const TInt32List &xInds, Int32 velA_idx, Int32 velE_idx,
    Int32 lbdaOffset_idx, Float64 normFactor)
    : m_fitter(fitter), m_EltsIdx(EltsIdx), m_lineType(lineType),
      m_redshift(redshift), m_xInds(xInds), m_velA_idx(velA_idx),
      m_velE_idx(velE_idx), m_lbdaOffset_idx(lbdaOffset_idx),
      m_normFactor(normFactor),
      m_spectralAxis(fitter.m_inputSpc.GetSpectralAxis()),
      m_noContinuumFluxAxis(fitter.m_model->getSpcFluxAxisNoContinuum()),
      m_continuumFluxAxis(fitter.m_model->getContinuumFluxAxis()),
      m_ErrorNoContinuum(fitter.m_inputSpc.GetErrorAxis()) {}

// compute Sum Squared diff and gradient
Float64 CLbfgsFitter::CLeastSquare::operator()(const VectorXd &x,
                                               VectorXd &grad) {

  // unpack param vector
  TFloat64List amps(x.begin(), x.begin() + m_EltsIdx.size());
  for (Int32 i = 0; i < m_EltsIdx.size(); ++i) {
    auto &elt = m_fitter.m_Elements[m_EltsIdx[i]];
    elt->SetFittedAmplitude(amps[i], 0.0);
  }

  Int32 lambdaOffset_idx = m_EltsIdx.size() + 1;

  Float64 velocityA = NAN;
  Float64 velocityE = NAN;
  if (m_lineType == CLine::nType_All) {
    velocityA = x[m_velA_idx];
    velocityE = x[m_velE_idx];
    lambdaOffset_idx++;
    for (Int32 eltIndex : m_EltsIdx) {
      auto &elt = m_fitter.m_Elements[eltIndex];
      if (elt->GetElementType() == CLine::nType_Absorption) {
        elt->SetVelocityAbsorption(velocityA);
      } else {
        elt->SetVelocityEmission(velocityE);
      }
    }
  } else if (m_lineType == CLine::nType_Absorption) {
    velocityA = x[m_velA_idx];
    for (Int32 eltIndex : m_EltsIdx) {
      auto &elt = m_fitter.m_Elements[eltIndex];
      elt->SetVelocityAbsorption(velocityA);
    }
  } else {
    velocityE = x[m_velE_idx];
    for (Int32 eltIndex : m_EltsIdx) {
      auto &elt = m_fitter.m_Elements[eltIndex];
      elt->SetVelocityEmission(velocityE);
    }
  }

  Float64 redshift = m_redshift + (1.0 + m_redshift) * x[lambdaOffset_idx];

  // set redshift in SYMIGM profiles
  for (Int32 eltIndex : m_EltsIdx) {
    auto &elt = m_fitter.m_Elements[eltIndex];
    const TInt32List &igm_indices = elt->getIgmLinesIndices();
    for (Int32 idx : igm_indices) {
      auto igmp = elt->GetSymIgmParams();
      igmp.m_redshift = redshift;
      elt->SetSymIgmParams(igmp, idx);
    }
  }

  TFloat64List pCoeffs;
  size_t ncoeff = 0;
  if (m_fitter.m_enableAmplitudeOffsets) {
    ncoeff = m_fitter.m_AmplitudeOffsetsDegree + 1;
    auto ibegin = x.begin() + lambdaOffset_idx + 1;
    pCoeffs = TFloat64List(ibegin, ibegin + ncoeff);
  }

  return ComputeLeastSquare(amps, redshift, pCoeffs, grad);
}

Float64 CLbfgsFitter::CLeastSquare::ComputeLeastSquare(const TFloat64List &amps,
                                                       Float64 redshift,
                                                       TFloat64List pCoeffs,
                                                       VectorXd &grad) {

  // compute least square term
  Float64 sumSquare = 0.;
  for (Int32 i = 0; i < m_xInds.size(); i++) {
    Float64 xi, yi, ei, ei2;
    Int32 idx = m_xInds[i];
    xi = m_spectralAxis[idx];
    yi = m_noContinuumFluxAxis[idx] * m_normFactor;
    ei = m_ErrorNoContinuum[idx] * m_normFactor;
    ei2 = ei * ei;

    // compute model value and gradient
    Float64 fval = 0.;
    TFloat64List ampsGrad;
    Float64 velocityAGrad = 0.0;
    Float64 velocityEGrad = 0.0;
    Float64 offsetGrad = 0.;
    auto amps_it = amps.cbegin();
    for (auto &eltIndex : m_EltsIdx) {
      auto &elt = m_fitter.m_Elements[eltIndex];
      // linemodel value &  amplitude derivative
      Float64 mval =
          elt->getModelAtLambda(xi, redshift, m_continuumFluxAxis[idx]);
      ampsGrad.push_back(mval / *amps_it);
      ++amps_it;
      fval += mval;

      // velocity derivative
      if (elt->GetElementType() == CLine::nType_Absorption) {
        velocityAGrad += elt->GetModelDerivVelAtLambda(
            xi, redshift, m_continuumFluxAxis[idx]);
      } else {
        velocityEGrad += elt->GetModelDerivVelAtLambda(
            xi, redshift, m_continuumFluxAxis[idx]);
      }

      // offset derivative
      offsetGrad +=
          elt->GetModelDerivZAtLambda(xi, redshift, m_continuumFluxAxis[idx]);
    }

    TFloat64List pCoeffGrad;
    if (m_fitter.m_enableAmplitudeOffsets) {
      pCoeffGrad.reserve(pCoeffs.size());
      // polynome value
      Float64 polynome = 0.;
      Float64 xipower = 1.;
      for (auto coeff : pCoeffs) {
        polynome += coeff * xipower;
        pCoeffGrad.push_back(xipower); // polynome derivative
        xipower *= xi;
      }
      fval += polynome;
    }

    // add squared diff
    sumSquare += (fval * fval - 2.0 * yi * fval) / ei2;
    Float64 residual = -2.0 * (yi - fval) / ei2;
    for (Int32 eltIndex = 0; eltIndex < m_EltsIdx.size(); ++eltIndex) {
      // squared diff derivative wrt amplitudes
      grad[eltIndex] += residual * ampsGrad[eltIndex];

      // squared diff derivative wrt velocity
      if (m_fitter.m_Elements[eltIndex]->GetElementType() ==
          CLine::nType_Absorption) {
        grad[m_velA_idx] += residual * velocityAGrad;
      } else {
        grad[m_velE_idx] += residual * velocityEGrad;
      }

      // squared diff derivative wrt offset
      grad[m_lbdaOffset_idx] += residual * offsetGrad;
    }

    // squared diff derivative wrt polynome coeffs
    if (m_fitter.m_enableAmplitudeOffsets) {
      for (size_t i = 0; i != pCoeffs.size(); ++i)
        grad[m_lbdaOffset_idx + i + 1] += residual * pCoeffGrad[i];
    }
  }

  return sumSquare;
}

void CLbfgsFitter::fit(Float64 redshift) {

  // fit the amplitudes of each element independently, unless there is overlap
  CHybridFitter::fitAmplitudesHybrid(redshift);
}

// since fitAmplitudesLinSolve is already fitting lambda offsets,
// this is a simple wrapper
Int32 CLbfgsFitter::fitAmplitudesLinSolveAndLambdaOffset(
    TInt32List EltsIdx, std::vector<Float64> &ampsfitted,
    std::vector<Float64> &errorsfitted, bool enableOffsetFitting,
    Float64 redshift) {

  Int32 ret =
      fitAmplitudesLinSolve(EltsIdx, ampsfitted, errorsfitted, redshift);

  return ret;
}

// overriding the SVD linear fitting, but here it is not linear inversion
Int32 CLbfgsFitter::fitAmplitudesLinSolve(const TInt32List &EltsIdx,
                                          TFloat64List &ampsfitted,
                                          TFloat64List &errorsfitted,
                                          Float64 redshift) {

  Int32 nddl = EltsIdx.size() + 2; // + velocity and lambda offset
  Int32 lbdaOffset_idx =
      EltsIdx.size() + 1; // position of lambda offset in the param vector

  // check if mixed types (Abs & Em)
  Int32 velA_idx = undefIdx;
  Int32 velE_idx = undefIdx;
  auto lineType = m_Elements[EltsIdx.front()]->GetElementType();
  if (lineType = CLine::nType_Absorption)
    velA_idx = EltsIdx.size();
  else
    velE_idx = EltsIdx.size();
  for (Int32 eltIndex = 1; eltIndex < EltsIdx.size(); ++eltIndex)
    if (lineType != m_Elements[eltIndex]->GetElementType()) {
      lineType = CLine::nType_All;
      nddl++; // 2 velocity parameter
      velA_idx = EltsIdx.size();
      velE_idx = velA_idx + 1;
      lbdaOffset_idx++;
    }

  TInt32List xInds = m_Elements.getSupportIndexes(EltsIdx);
  if (xInds.size() < 1)
    return -1;

  Int32 n = xInds.size();
  if (n < nddl) {
    ampsfitted.resize(nddl);
    errorsfitted.resize(nddl);
    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      ampsfitted[iddl] = NAN;
      errorsfitted[iddl] = NAN;
    }
    return -1;
  }

  Int32 idxAmpOffset = undefIdx;
  if (m_enableAmplitudeOffsets) {
    nddl += m_AmplitudeOffsetsDegree + 1;
    // find the amplitudeOffset Support that corresponds to these elts
    idxAmpOffset = m_Elements.getIndexAmpOffset(xInds[0]);
  }

  // Normalize
  // Float64 maxabsval = DBL_MIN;
  const auto &noContinuumFluxAxis = m_model->getSpcFluxAxisNoContinuum();
  Int32 maxabsval_idx =
      *std::max_element(xInds.cbegin(), xInds.cend(), [&](Int32 l, Int32 r) {
        return std::abs(noContinuumFluxAxis[l]) <
               std::abs(noContinuumFluxAxis[r]);
      });
  Float64 normFactor = 1.0 / std::abs(noContinuumFluxAxis[maxabsval_idx]);

  // velocity bounds
  VectorXd v_velocityMin;
  VectorXd v_velocityMax;
  if (lineType == CLine::nType_All) {
    v_velocityMin << m_velfitMinA, m_velfitMinE;
    v_velocityMax << m_velfitMaxA, m_velfitMaxE;
  } else if (lineType == CLine::nType_Absorption) {
    v_velocityMin << m_velfitMinA;
    v_velocityMax << m_velfitMaxA;
  } else {
    v_velocityMin << m_velfitMinE;
    v_velocityMax << m_velfitMaxE;
  }

  // amplitude bounds
  VectorXd v_ampsMax(EltsIdx.size());
  VectorXd v_ampsMin = VectorXd::Zero(EltsIdx.size());
  for (Int32 eltIndex = 1; eltIndex < EltsIdx.size(); ++eltIndex) {
    Float64 ampMax = INFINITY;
    auto &elt = m_Elements[eltIndex];
    if (elt->GetElementType() == CLine::nType_Absorption &&
        elt->GetAbsLinesLimit() > 0.0)
      ampMax = elt->GetAbsLinesLimit() / elt->GetMaxNominalAmplitude();
    v_ampsMax << ampMax;
  }

  // offset bounds
  Float64 offsetMin = m_LambdaOffsetMin * SPEED_OF_LIGHT_IN_VACCUM;
  Float64 offsetMax = m_LambdaOffsetMax * SPEED_OF_LIGHT_IN_VACCUM;

  // all bounds
  VectorXd lb(nddl);
  VectorXd ub(nddl);
  lb << v_ampsMin, v_velocityMin, offsetMin;
  ub << v_ampsMax, v_velocityMax, offsetMax;

  // polynome bounds
  if (m_enableAmplitudeOffsets) {
    size_t ncoeff = m_AmplitudeOffsetsDegree + 1;
    VectorXd v_pCoeffMin(ncoeff);
    VectorXd v_pCoeffMax(ncoeff);
    for (size_t i = 0; i < ncoeff; ++i) {
      v_pCoeffMin << -INFINITY;
      v_pCoeffMax << INFINITY;
    }
    lb << v_pCoeffMin;
    ub << v_pCoeffMax;
  }

  // amplitudes initial guess
  TFloat64List ampsGuess;
  ampsGuess.reserve(EltsIdx.size());
  for (Int32 eltIndex : EltsIdx) {
    auto &elt = m_Elements[eltIndex];
    // set guess velocity
    if (elt->GetElementType() == CLine::nType_Absorption) {
      elt->SetVelocityAbsorption(m_velIniGuessA);
    } else {
      elt->SetVelocityEmission(m_velIniGuessE);
    }
    fitAmplitude(eltIndex, redshift);
    ampsGuess.push_back(elt->GetElementAmplitude() * normFactor);
    if (std::isnan(ampsGuess.back()))
      return -1;
  }

  // velocity initial guess
  VectorXd v_ampsGuess = VectorXd::Map(ampsGuess.data(), ampsGuess.size());
  VectorXd v_velocityGuess;
  if (lineType == CLine::nType_All)
    v_velocityGuess << m_velIniGuessA, m_velIniGuessE;
  else if (lineType == CLine::nType_Absorption)
    v_velocityGuess << m_velIniGuessA;
  else
    v_velocityGuess << m_velIniGuessE;

  VectorXd v_xGuess(nddl);
  v_xGuess << v_ampsGuess, v_velocityGuess, 0.;

  if (m_enableAmplitudeOffsets)
    v_xGuess << VectorXd::Zero(m_AmplitudeOffsetsDegree + 1);

  // Set up solver parameters
  LBFGSBParam<Float64> param; // New parameter class
  param.epsilon = 1e-6;
  param.max_iterations = 100;

  // Create solver and function object
  LBFGSBSolver<Float64> solver(param); // New solver class
  CLeastSquare fun(*this, EltsIdx, lineType, redshift, xInds, velA_idx,
                   velE_idx, lbdaOffset_idx, normFactor);

  // v_xGuess will be overwritten to be the best point found
  Float64 fx;
  int niter = solver.minimize(fun, v_xGuess, fx, lb, ub);

  // store fitted amplitude
  // TODO SNR computation
  Float64 snr = NAN;
  ampsfitted.assign(EltsIdx.size(), NAN);
  for (Int32 i = 0; i < EltsIdx.size(); ++i) {
    ampsfitted[i] = v_xGuess[i] / normFactor;
    m_Elements[EltsIdx[i]]->SetFittedAmplitude(ampsfitted[i], snr);
  }

  // store fitted velocity dispersion (line width)
  if (lineType == CLine::nType_All) {
    Float64 velocityA = v_xGuess[velA_idx];
    Float64 velocityE = v_xGuess[velE_idx];
    for (Int32 i = 0; i < EltsIdx.size(); ++i) {
      auto &elt = m_Elements[EltsIdx[i]];
      if (elt->GetElementType() == CLine::nType_Absorption)
        elt->SetVelocityAbsorption(velocityA);
      else
        elt->SetVelocityEmission(velocityE);
    }
  } else if (lineType == CLine::nType_Absorption) {
    Float64 velocity = v_xGuess[velA_idx];
    for (Int32 i = 0; i < EltsIdx.size(); ++i)
      m_Elements[EltsIdx[i]]->SetVelocityAbsorption(velocity);
  } else {
    Float64 velocity = v_xGuess[velE_idx];
    for (Int32 i = 0; i < EltsIdx.size(); ++i)
      m_Elements[EltsIdx[i]]->SetVelocityEmission(velocity);
  }

  // store velocity offset (line offset)
  for (Int32 eltIndex : EltsIdx) {
    auto &elt = m_Elements[eltIndex];
    for (Int32 i = 0; i < elt->GetSize(); ++i) {
      Float64 offset = elt->GetOffset(i);
      offset *= SPEED_OF_LIGHT_IN_VACCUM;
      offset += (1 + offset) * v_xGuess[lbdaOffset_idx];
      offset /= SPEED_OF_LIGHT_IN_VACCUM;
      elt->SetAllOffsets(offset);
    }
  }

  // TODO store polynome coeffs (continuum under line)
  if (m_enableAmplitudeOffsets) {
    Float64 x0 = v_xGuess[lbdaOffset_idx + 1] / normFactor;
    Float64 x1 = 0.0;
    Float64 x2 = 0.0;
    if (m_AmplitudeOffsetsDegree > 0)
      x1 = v_xGuess[lbdaOffset_idx + 2] / normFactor;
    if (m_AmplitudeOffsetsDegree > 1)
      x2 = v_xGuess[lbdaOffset_idx + 3] / normFactor;
    m_Elements.setAmplitudeOffsetsCoeffsAt(idxAmpOffset, {x0, x1, x2});
  }

  return EltsIdx.size(); // nb of amplitudes of same sign (all >0)
}

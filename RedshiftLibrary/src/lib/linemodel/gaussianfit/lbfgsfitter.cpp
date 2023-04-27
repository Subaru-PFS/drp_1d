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
#include "RedshiftLibrary/linemodel/gaussianfit/LineSearchMoreThuente.h"
#include <unsupported/Eigen/NumericalDiff>

using namespace std;
using namespace Eigen;
using namespace LBFGSpp;
using namespace NSEpic;

CLbfgsFitter::CLeastSquare::CLeastSquare(
    CLbfgsFitter &fitter, const TInt32List &EltsIdx, Int32 lineType,
    Float64 redshift, const TInt32List &xInds, Int32 velA_idx, Int32 velE_idx,
    Int32 lbdaOffset_idx, Int32 pCoeff_idx, Float64 normFactor)
    : Functor<Float64, Dynamic, 1>(), m_fitter(&fitter), m_EltsIdx(&EltsIdx),
      m_lineType(lineType), m_redshift(redshift), m_xInds(&xInds),
      m_velA_idx(velA_idx), m_velE_idx(velE_idx),
      m_lbdaOffset_idx(lbdaOffset_idx), m_pCoeff_idx(pCoeff_idx),
      m_normFactor(normFactor),
      m_spectralAxis(&fitter.m_inputSpc.GetSpectralAxis()),
      m_noContinuumFluxAxis(&fitter.m_model->getSpcFluxAxisNoContinuum()),
      m_continuumFluxAxis(&fitter.m_model->getContinuumFluxAxis()),
      m_ErrorNoContinuum(&fitter.m_inputSpc.GetErrorAxis()) {

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
void CLbfgsFitter::CLeastSquare::operator()(const VectorXd &x,
                                            ValueType &val) const {

  TFloat64List amps;
  Float64 redshift;
  TPolynomCoeffs pCoeffs;
  std::tie(amps, redshift, pCoeffs) = unpack(x);

  val(0, 0) = ComputeLeastSquare(amps, redshift, pCoeffs);

  Log.LogDebug(Formatter() << "LeastSquare = " << val(0, 0));
}

/*
// compute Sum Squared diff and gradient
CLbfgsFitter::CLeastSquare::ValueType
CLbfgsFitter::CLeastSquare::operator()(const VectorXd &x, VectorXd &grad) {

  TFloat64List amps;
  Float64 redshift;
  TFloat64List pCoeffs;
  std::tie(amps, redshift, pCoeffs) = unpack(x);

  Float64 lsq = ComputeLeastSquareAndGrad(amps, redshift, pCoeffs, grad);

  Log.LogDebug(Formatter() << "LeastSquare = " << lsq);
  for (size_t i = 0; i < m_EltsIdx->size(); ++i)
    Log.LogDebug(Formatter() << "LeastSquare grad amp[" << i <<"] = " <<
grad[i]); if (m_lineType == CLine::nType_All) { Log.LogDebug(Formatter() <<
"LeastSquare grad velA= " << grad[m_velA_idx]); Log.LogDebug(Formatter() <<
"LeastSquare grad velE= " << grad[m_velE_idx]); } else if (m_lineType ==
CLine::nType_Absorption) Log.LogDebug(Formatter() << "LeastSquare grad velA= "
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

std::tuple<TFloat64List, Float64, TPolynomCoeffs>
CLbfgsFitter::CLeastSquare::unpack(const VectorXd &x) const {

  // unpack param vector
  TFloat64List amps(x.begin(), x.begin() + m_EltsIdx->size());
  for (Int32 i = 0; i < m_EltsIdx->size(); ++i) {
    Log.LogDebug(Formatter() << "amplitude[" << i << "]= " << amps[i]);
    auto &elt = m_fitter->m_Elements[(*m_EltsIdx)[i]];
    elt->SetElementAmplitude(amps[i], 0.0);
  }

  if (m_fitter->m_enableVelocityFitting) {

    if (m_velA_idx != undefIdx)
      Log.LogDebug(Formatter() << "velocity Abs = " << x[m_velA_idx]);
    if (m_velE_idx != undefIdx)
      Log.LogDebug(Formatter() << "velocity Em  = " << x[m_velE_idx]);

    for (Int32 eltIndex : *m_EltsIdx) {
      auto &elt = m_fitter->m_Elements[eltIndex];
      if (elt->GetIsEmission())
        elt->SetVelocityEmission(x[m_velE_idx]);
      else
        elt->SetVelocityAbsorption(x[m_velA_idx]);
    }
  }

  Float64 redshift = m_redshift;
  if (m_fitter->m_enableLambdaOffsetsFit) {
    redshift += (1.0 + m_redshift) * x[m_lbdaOffset_idx];
    Log.LogDebug(Formatter() << "redshift=" << redshift);
  }

  // set redshift in SYMIGM profiles
  for (Int32 eltIndex : *m_EltsIdx) {
    auto &elt = m_fitter->m_Elements[eltIndex];
    const TInt32List &igm_indices = elt->getIgmLinesIndices();
    for (Int32 idx : igm_indices) {
      auto igmp = elt->GetSymIgmParams();
      igmp.m_redshift = redshift;
      elt->SetSymIgmParams(igmp, idx);
    }
  }

  TPolynomCoeffs pCoeffs;
  if (m_fitter->m_enableAmplitudeOffsets) {
    size_t ncoeff = 0;
    ncoeff = TPolynomCoeffs::degree + 1;
    auto ibegin = x.begin() + m_pCoeff_idx;
    pCoeffs = TPolynomCoeffs(TFloat64List(ibegin, ibegin + ncoeff));
    Log.LogDebug(Formatter() << "p coeffs  = " << pCoeffs.a0 << " "
                             << pCoeffs.a1 << " " << pCoeffs.a2);
  }

  return std::make_tuple(std::move(amps), redshift, std::move(pCoeffs));
}

Float64 CLbfgsFitter::CLeastSquare::ComputeLeastSquare(
    const TFloat64List &amps, Float64 redshift,
    const TPolynomCoeffs &pCoeffs) const {
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
      auto &elt = m_fitter->m_Elements[eltIndex];
      // linemodel value
      Float64 mval =
          elt->getModelAtLambda(xi, redshift, (*m_continuumFluxAxis)[idx]);
      fval += mval;
    }

    if (m_fitter->m_enableAmplitudeOffsets)
      fval += pCoeffs.getValue(xi);

    // add squared diff
    sumSquare += (fval * fval - 2.0 * yi * fval) / ei2;
  }

  return sumSquare;
}

Float64 CLbfgsFitter::CLeastSquare::ComputeLeastSquareAndGrad(
    const TFloat64List &amps, Float64 redshift, const TPolynomCoeffs &pCoeffs,
    VectorXd &grad) const {

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
      auto &elt = m_fitter->m_Elements[eltIndex];
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
        if (elt->GetElementType() == CLine::nType_Absorption) {
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
        if (m_fitter->m_Elements[eltIndex]->GetElementType() ==
            CLine::nType_Absorption) {
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
      for (size_t i = 0; i <= pCoeffs.degree; ++i)
        grad[m_pCoeff_idx + i] += residual * pCoeffGrad[i];
    }
  }

  return sumSquare;
}

void CLbfgsFitter::resetSupport(Float64 redshift) {

  // set velocity at max value (to set largest line overlapping)
  if (Context.GetParameterStore()->GetScoped<bool>("linemodel.velocityfit")) {
    const Float64 velfitMaxE = Context.GetParameterStore()->GetScoped<Float64>(
        "linemodel.emvelocityfitmax");
    const Float64 velfitMaxA = Context.GetParameterStore()->GetScoped<Float64>(
        "linemodel.absvelocityfitmax");

    for (Int32 j = 0; j < m_Elements.size(); j++) {
      m_ElementParam[j]->m_VelocityEmission = velfitMaxE;
      m_ElementParam[j]->m_VelocityAbsorption = velfitMaxA;
    }
  }

  CAbstractFitter::resetSupport(redshift);
}

void CLbfgsFitter::doFit(Float64 redshift) {

  // fit the amplitudes of each element independently, unless there is overlap
  CHybridFitter::fitAmplitudesHybrid(redshift);
}

// since fitAmplitudesLinSolve is already fitting lambda offsets,
// this is a simple wrapper
void CLbfgsFitter::fitAmplitudesLinSolveAndLambdaOffset(
    TInt32List EltsIdx, bool enableOffsetFitting, Float64 redshift) {

  fitAmplitudesLinSolvePositive(EltsIdx, redshift);
}

// overriding the SVD linear fitting, but here it is not linear inversion
void CLbfgsFitter::fitAmplitudesLinSolvePositive(const TInt32List &EltsIdx,
                                                 Float64 redshift) {

  if (EltsIdx.size() < 1)
    THROWG(INTERNAL_ERROR, "empty Line element list to fit");

  Int32 nddl = EltsIdx.size();
  Int32 param_idx = nddl;
  Int32 velA_idx = undefIdx;
  Int32 velE_idx = undefIdx;
  Int32 lineType = 0;
  if (m_enableVelocityFitting) {
    ++nddl;
    // position of velocity in the param vector
    Int32 velocity_param_idx = param_idx++;

    // check if mixed types (Abs & Em)
    lineType = m_Elements[EltsIdx.front()]->GetElementType();
    if (lineType == CLine::nType_Absorption)
      velA_idx = velocity_param_idx;
    else
      velE_idx = velocity_param_idx;
    for (Int32 eltIndex : EltsIdx)
      if (lineType != m_Elements[eltIndex]->GetElementType()) {
        lineType = CLine::nType_All;
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

  TInt32List xInds = m_Elements.getSupportIndexes(EltsIdx);
  if (xInds.size() < 1)
    THROWG(INTERNAL_ERROR, "no observed samples for the line Element to fit");

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
    Flag.warning(WarningCode::LINEARFIT_RANK_DEFICIENT,
                 Formatter() << __func__ << " LBFGS ill ranked:"
                             << " number of samples = " << n
                             << ", number of parameters to fit = " << nddl);
    for (Int32 eltIndex : EltsIdx)
      m_Elements.SetElementAmplitude(eltIndex, 0., INFINITY);
    if (m_enableAmplitudeOffsets) {
      for (Int32 eltIndex : EltsIdx)
        m_Elements[eltIndex]->SetPolynomCoeffs({0., 0., 0.});
    }
    return;
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

  // bounds for all params
  /////////////////////////
  VectorXd lb(nddl);
  VectorXd ub(nddl);

  // amplitude bounds
  const auto amp_indices = seq(0, EltsIdx.size() - 1);
  lb(amp_indices) = VectorXd::Zero(EltsIdx.size());
  for (size_t i = 0; i < EltsIdx.size(); ++i) {
    Float64 ampMax = INFINITY;
    auto &elt = m_Elements[EltsIdx[i]];
    if (elt->GetElementType() == CLine::nType_Absorption &&
        elt->GetAbsLinesLimit() > 0.0)
      ampMax = elt->GetAbsLinesLimit() / elt->GetMaxNominalAmplitude();
    ub[i] = ampMax;
  }

  // velocity bounds
  if (m_enableVelocityFitting) {
    if (lineType == CLine::nType_All) {
      const auto vel_indices = seq(velA_idx, velE_idx);
      lb(vel_indices) << m_velfitMinA, m_velfitMinE;
      ub(vel_indices) << m_velfitMaxA, m_velfitMaxE;
    } else if (lineType == CLine::nType_Absorption) {
      lb[velA_idx] = m_velfitMinA;
      ub[velA_idx] = m_velfitMaxA;
    } else {
      lb[velE_idx] = m_velfitMinE;
      ub[velE_idx] = m_velfitMaxE;
    }
  }

  if (m_enableLambdaOffsetsFit) {
    // offset bounds
    lb[lbdaOffset_param_idx] = m_LambdaOffsetMin / SPEED_OF_LIGHT_IN_VACCUM;
    ub[lbdaOffset_param_idx] = m_LambdaOffsetMax / SPEED_OF_LIGHT_IN_VACCUM;
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

  // amplitudes initial guess
  for (auto eltIndex : EltsIdx) {
    auto &elt = m_Elements[eltIndex];
    // set velocity guess
    if (elt->GetElementType() == CLine::nType_Absorption) {
      elt->SetVelocityAbsorption(m_velIniGuessA);
    } else {
      elt->SetVelocityEmission(m_velIniGuessE);
    }
  }
  CSvdFitter::fitAmplitudesLinSolvePositive(EltsIdx, redshift);
  for (size_t i = 0; i != EltsIdx.size(); ++i) {
    // fitAmplitude(EltsIdx[i], redshift);
    auto &elt = m_Elements[EltsIdx[i]];
    v_xGuess[i] = elt->GetElementAmplitude() * normFactor;
    if (std::isnan(v_xGuess[i]))
      THROWG(INTERNAL_ERROR, "NAN amplitude");
  }
  // polynomial coeffs initial guess
  if (m_enableAmplitudeOffsets) {
    const auto &pCoeffs = m_Elements[EltsIdx.front()]->GetPolynomCoeffs();
    v_xGuess[pCoeff_param_idx] = pCoeffs.a0 * normFactor;
    v_xGuess[pCoeff_param_idx + 1] = pCoeffs.a1 * normFactor;
    v_xGuess[pCoeff_param_idx + 2] = pCoeffs.a2 * normFactor;
  }

  // velocity initial guess
  if (m_enableVelocityFitting) {
    if (lineType == CLine::nType_All) {
      const auto vel_indices = seq(velA_idx, velE_idx);
      v_xGuess(vel_indices) << m_velIniGuessA, m_velIniGuessE;
    } else if (lineType == CLine::nType_Absorption)
      v_xGuess[velA_idx] = m_velIniGuessA;
    else
      v_xGuess[velE_idx] = m_velIniGuessE;
    // v_xGuess[velE_idx] =  50.0;// DEBUGGING
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
  LBFGSBSolver<Float64, LineSearchMoreThuenteFork> solver(param);

  CLeastSquare func(*this, EltsIdx, lineType, redshift, xInds, velA_idx,
                    velE_idx, lbdaOffset_param_idx, pCoeff_param_idx,
                    normFactor);

  NumericalDiff<CLeastSquare> func_with_num_diff(func, 1e-12);

  auto myfunc = [&func_with_num_diff](const VectorXd &x, VectorXd &grad) {
    CLeastSquare::JacobianType J(1, x.size());
    func_with_num_diff.df(x, J);
    grad = J;
    CLeastSquare::ValueType val;
    func_with_num_diff(x, val);
    return val(0, 0);
  };

  // v_xGuess will be overwritten to be the best point found
  Float64 fx;
  int niter = solver.minimize(myfunc, v_xGuess, fx, lb, ub);

  // store fitted amplitude
  // TODO SNR computation
  Float64 snr = NAN;
  for (Int32 i = 0; i < EltsIdx.size(); ++i)
    m_Elements[EltsIdx[i]]->SetElementAmplitude(v_xGuess[i] / normFactor, snr);

  // store fitted velocity dispersion (line width)
  if (m_enableVelocityFitting) {
    if (lineType == CLine::nType_All) {
      Float64 velocityA = v_xGuess[velA_idx];
      Float64 velocityE = v_xGuess[velE_idx];
      for (Int32 eltIndex : EltsIdx) {
        auto &elt = m_Elements[eltIndex];
        if (elt->GetElementType() == CLine::nType_Absorption)
          elt->SetVelocityAbsorption(velocityA);
        else
          elt->SetVelocityEmission(velocityE);
      }
    } else if (lineType == CLine::nType_Absorption) {
      Float64 velocity = v_xGuess[velA_idx];
      for (Int32 eltIndex : EltsIdx)
        m_Elements[eltIndex]->SetVelocityAbsorption(velocity);
    } else {
      Float64 velocity = v_xGuess[velE_idx];
      for (Int32 eltIndex : EltsIdx)
        m_Elements[eltIndex]->SetVelocityEmission(velocity);
    }
  }

  // store velocity offset (line offset)
  if (m_enableLambdaOffsetsFit) {
    for (Int32 eltIndex : EltsIdx) {
      auto &elt = m_Elements[eltIndex];
      for (Int32 i = 0; i < elt->GetSize(); ++i) {
        Float64 offset = elt->GetOffset(i);
        offset /= SPEED_OF_LIGHT_IN_VACCUM;
        offset += (1 + offset) * v_xGuess[lbdaOffset_param_idx];
        offset *= SPEED_OF_LIGHT_IN_VACCUM;
        elt->SetAllOffsets(offset);
      }
    }
  }

  // store polynome coeffs (continuum under line)
  if (m_enableAmplitudeOffsets) {
    Float64 x0 = v_xGuess[pCoeff_param_idx] / normFactor;
    Float64 x1 = 0.0;
    Float64 x2 = 0.0;
    if (TPolynomCoeffs::degree > 0)
      x1 = v_xGuess[pCoeff_param_idx + 1] / normFactor;
    if (TPolynomCoeffs::degree > 1)
      x2 = v_xGuess[pCoeff_param_idx + 2] / normFactor;
    for (Int32 eltIndex : EltsIdx)
      m_Elements[eltIndex]->SetPolynomCoeffs({x0, x1, x2});
  }
}

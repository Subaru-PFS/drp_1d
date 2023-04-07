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
#ifndef _REDSHIFT_LBFGS_FITTER_
#define _REDSHIFT_LBFGS_FITTER_

#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/processflow/context.h"

#include <Eigen/Core>
#include <LBFGSB.h>

using Eigen::VectorXd;

namespace NSEpic {

// Generic functor
// See http://eigen.tuxfamily.org/index.php?title=Functors
// C++ version of a function pointer that stores meta-data about the function
template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor {

  // Information that tells the caller the numeric type (eg. double) and size
  // (input / output dim)
  typedef _Scalar Scalar;
  enum { // Required by numerical differentiation module
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };

  // Tell the caller the matrix sizes associated with the input, output, and
  // jacobian
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime>
      JacobianType;

  // Local copy of the number of inputs
  int m_inputs, m_values;

  // Two constructors:
  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

class CLbfgsFitter : public CHybridFitter {
public:
  class CLeastSquare : public Functor<Float64, Eigen::Dynamic, 1> {

  public:
    CLeastSquare(CLbfgsFitter &fitter, const TInt32List &EltsIdx,
                 Int32 lineType, Float64 redshift, const TInt32List &xInds,
                 Int32 velA_idx, Int32 velE_idx, Int32 lbdaOffset_idx,
                 Int32 pCoeff_indx, Float64 normFactor);

    // Float64 operator()(const VectorXd &x, VectorXd &grad);
    void operator()(const VectorXd &x, ValueType &retvalue) const;

  private:
    std::tuple<TFloat64List, Float64, TFloat64List>
    unpack(const VectorXd &x) const;
    Float64 ComputeLeastSquare(const TFloat64List &amps, Float64 redshift,
                               TFloat64List pCoeffs) const;

    Float64 ComputeLeastSquareAndGrad(const TFloat64List &amps,
                                      Float64 redshift, TFloat64List pCoeffs,
                                      VectorXd &grad) const;

    CLbfgsFitter *m_fitter;
    const TInt32List *m_EltsIdx;
    const Int32 m_lineType;
    const Float64 m_redshift;
    const TInt32List *m_xInds;
    const Int32 m_velA_idx;
    const Int32 m_velE_idx;
    const Int32 m_lbdaOffset_idx;
    const Int32 m_pCoeff_idx;
    const Float64 m_normFactor;
    const CSpectrumSpectralAxis *m_spectralAxis;
    const CSpectrumFluxAxis *m_noContinuumFluxAxis;
    const CSpectrumFluxAxis *m_continuumFluxAxis;
    const CSpectrumNoiseAxis *m_ErrorNoContinuum;
    Float64 m_sumSquareData;
  };

  using CHybridFitter::CHybridFitter;

  void fit(Float64 redshift) override;

  void fitAmplitudesLinSolveAndLambdaOffset(TInt32List EltsIdx,
                                            bool enableOffsetFitting,
                                            Float64 redshift) override;

  void fitAmplitudesLinSolvePositive(const TInt32List &EltsIdx,
                                     Float64 redshift) override;

private:
  const bool m_enableVelocityFitting =
      Context.GetParameterStore()->GetScoped<bool>("linemodel.velocityfit");

  const Float64 m_velfitMinE = Context.GetParameterStore()->GetScoped<Float64>(
      "linemodel.emvelocityfitmin");
  const Float64 m_velfitMaxE = Context.GetParameterStore()->GetScoped<Float64>(
      "linemodel.emvelocityfitmax");
  const Float64 m_velfitMinA = Context.GetParameterStore()->GetScoped<Float64>(
      "linemodel.absvelocityfitmin");
  const Float64 m_velfitMaxA = Context.GetParameterStore()->GetScoped<Float64>(
      "linemodel.absvelocityfitmax");

  const Float64 m_velIniGuessE =
      Context.GetParameterStore()->GetScoped<Float64>(
          "linemodel.velocityemission");
  const Float64 m_velIniGuessA =
      Context.GetParameterStore()->GetScoped<Float64>(
          "linemodel.velocityabsorption");
};
} // namespace NSEpic
#endif

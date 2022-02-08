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

#include "RedshiftLibrary/method/bayesianSolve.h"

using namespace NSEpic;

void CObjectSolve::InitRanges(std::shared_ptr<const CInputContext> inputContext)
{
  m_lambdaRange=inputContext->m_lambdaRange;//non-clamped

  //m_redshiftSampling could be overwritten if fftprocessing is activated
  // Warning for LineMeas :  we consider linemeas use the same redshiftsampling as the objectMethod
  // if linemeas is called in standalone, redshiftsampling must be present in parameters
  m_redshiftSampling=inputContext->GetParameterStore()->GetScoped<std::string>("redshiftsampling");

  TFloat64Range redshiftRange;
  Float64 redshiftStep;
  GetRedshiftSampling(inputContext, redshiftRange, redshiftStep);
  Log.LogInfo(Formatter()<<"Init redshift range with " << redshiftRange << " and"<< redshiftStep);      
  if(m_redshiftSampling=="log")
    m_redshifts = redshiftRange.SpreadOverLogZplusOne( redshiftStep ); //experimental: spreadover a grid at delta/(1+z), unusable because PDF needs regular z-step
  else  
    m_redshifts = redshiftRange.SpreadOver( redshiftStep );

}


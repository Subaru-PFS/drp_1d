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
#include <cstdarg>
#include <iostream>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/ray/ruleBalmerLinearSolver.h"
#include "RedshiftLibrary/ray/linetags.h"

using namespace NSEpic;
using namespace std;

CRuleBalmerLinearSolver::CRuleBalmerLinearSolver()
{
}

CRuleBalmerLinearSolver::~CRuleBalmerLinearSolver()
{
}

void CRuleBalmerLinearSolver::SetUp( Bool EnabledArgument, ... )
{
  Name = "balmerlinesolve";
  Enabled = EnabledArgument;
  //va_list Arguments;
  //va_start ( Arguments, EnabledArgument );
}

Bool CRuleBalmerLinearSolver::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

/**
 * \brief Check for and correct if necessary the Balmer rule on the m_Elements.
 **/
void CRuleBalmerLinearSolver::Correct( CLineModelElementList& LineModelElementList )
{
    linetags ltags;
    const CSpectrumSpectralAxis& spectralAxis = LineModelElementList.m_SpectrumModel.GetSpectralAxis();
    std::vector<Float64> lambdax;
    std::vector<Float64> continuumx;
    std::vector<Float64> datax;
    std::vector<Float64> errdatax;
    std::vector<Int32> iEltE;
    std::vector<Int32> iEltA;
    std::vector<std::string> linetags;
    linetags.push_back( ltags.halpha_em );
    linetags.push_back( ltags.hbeta_em );
    linetags.push_back( ltags.hgamma_em );
    linetags.push_back( ltags.hdelta_em );
    linetags.push_back( ltags.h8_em );
    linetags.push_back( ltags.h9_em );
    linetags.push_back( ltags.h10_em );
    linetags.push_back( ltags.h11_em );
    Int32 ilineE;
    Int32 ilineA;
    std::string tagE;
    std::string tagA;
    TFloat64List ampsE;
    TFloat64List ersE;
    TFloat64List ersFitE;
    TFloat64List ampsA;
    TFloat64List ersA;
    TFloat64List ersFitA;
    TBoolList ampsEwasfitted;
    Int32 nLines = 0;
    for( Int32 itag=0; itag<linetags.size(); itag++)
    {
        tagE = linetags[itag];
        tagA = linetags[itag];
        tagA.append( "A" );
        ilineE = LineModelElementList.FindElementIndex( tagE, CRay::nType_Emission );
        ilineA = LineModelElementList.FindElementIndex( tagA, CRay::nType_Absorption );
        if( ilineE==-1 )
        {
            Log.LogDebug( "Rule %s: no element with name %s.", Name.c_str(), tagE.c_str() );
            continue;
        }
        if( ilineA==-1 )
        {
            Log.LogDebug( "Rule %s: no element with name %s.", Name.c_str(), tagA.c_str() );
            continue;
        }
        Float64 ampE = LineModelElementList.m_Elements[ilineE]->GetFittedAmplitude( 0 );
        Float64 erE = LineModelElementList.m_Elements[ilineE]->GetFittedAmplitudeErrorSigma( 0 );
        Float64 erFitE = LineModelElementList.getModelErrorUnderElement( ilineE );
        erE = sqrt(erE*erE+erFitE*erFitE);
        Float64 ampA = LineModelElementList.m_Elements[ilineA]->GetFittedAmplitude( 0 );
        Float64 erA = LineModelElementList.m_Elements[ilineA]->GetFittedAmplitudeErrorSigma( 0 );
        Float64 erFitA = LineModelElementList.getModelErrorUnderElement( ilineA );
        erA = sqrt(erA*erA+erFitA*erFitA);
        Float64 amp;
        Float64 er;
        if( ampE>0.0 && ampE>ampA )
        {
            amp = ampE;
            er = erE;
            ampsEwasfitted.push_back( true );
        }
        else
        {
            if( ampA>0.0 )
            {
                amp = -ampA;
                er = erA;
                ampsEwasfitted.push_back( false );
            }
            else
            {
                continue;
            }
        }
        nLines++;
        ampsE.push_back( ampE );
        ersE.push_back( erE );
        ersFitE.push_back( erFitE );
        ampsA.push_back( ampA );
        ersA.push_back( erA );
        ersFitA.push_back( erFitA );
        iEltE.push_back( ilineE );
        iEltA.push_back( ilineA );
        Float64 lambda = LineModelElementList.m_Elements[ilineE]->GetRays()[0].GetPosition()*(1.0+LineModelElementList.GetRedshift());
        lambdax.push_back( lambda );
        Int32 Idx = spectralAxis.GetIndexAtWaveLength( lambda );
        continuumx.push_back( LineModelElementList.m_inputSpc.GetContinuumFluxAxis()[Idx] );
        datax.push_back( LineModelElementList.m_inputSpc.GetContinuumFluxAxis()[Idx] + amp );
        errdatax.push_back(er);
    }

    /*
    TFloat64List coeffs = BalmerModelLinSolve( lambdax, continuumx, datax, errdatax );
    // export for debug
    FILE* fspc = fopen( "BalmerLinSolve_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1.0;
    for(Int32 i=0; i<nLines; i++){
        Float64 ampRegE = coeffs[0]*lambdax[i] + coeffs[3] ;
        Float64 coeffRegA = (coeffs[2]*lambdax[i] + coeffs[1]);
        Float64 coeffA = (continuumx[i] - ampsA[i])/continuumx[i];
        fprintf( fspc, "%d %f %f %f %f %f %f\n", i, lambdax[i], errdatax[i], ampRegE, ampsE[i], coeffRegA, coeffA);
    }
    fclose( fspc );
    //*/

    //apply the absorption rule
    for( Int32 j=0; j<nLines; j++)
    {
        Int32 i1=nLines-1-j; //go through the lines, small wavelengths first
        Int32 i2 = i1-1;
        if( i2<0 )
        {
            break;
        }
        Float64 coeffA1 = (continuumx[i1] - ampsA[i1])/continuumx[i1];
        Float64 coeffA2 = (continuumx[i2] - ampsA[i2])/continuumx[i2];
        Float64 erA1 = ersA[i1]; //todo, check if this error value can be calculated better !
        Float64 erA2 = ersA[i2]; //...

        if(coeffA1<coeffA2)
        { //A1 should absorb less than A2
            Float64 R = 1.0;
            Float64 wA1 = 0.0;
            if(erA1!=0.0)
            {
                wA1 = 1.0/(erA1*erA1);
            }
            Float64 wA2 = 0.0;
            if(erA2!=0.0)
            {
                wA2 = 1.0/(erA2*erA2*R*R);
            }
            Float64 correctedCoeffA1 = (coeffA1*wA1 + coeffA2*wA2*R)/(wA1+wA2);
            Float64 correctedCoeffA2 = correctedCoeffA1/R;
            {  //make AL and EL higher...
                Float64 correctedAmpA2 = (1.0 - correctedCoeffA2)*continuumx[i2];
                if(correctedAmpA2 > 0){
                    Float64 correction = correctedAmpA2-ampsA[i2];
                    Float64 widthRatioAE = 1.0/3.0;
                    ampsA[i2] = ampsA[i2]+correction;
                    if( ampsEwasfitted[i2] )
                    {
                        ampsE[i2] = ampsE[i2]+correction;
                        LineModelElementList.m_Elements[iEltE[i2]]->SetFittedAmplitude(ampsE[i2]+ampsA[i2]*widthRatioAE, ersE[i2]);
                    }
                    LineModelElementList.m_Elements[iEltA[i2]]->SetFittedAmplitude(ampsA[i2], ersA[i2]);
                    //than re-fit the lines together, because the correction applied so far is for R=1.0, it could be R>1.0
                    if(ampsEwasfitted[i2])
                    {
                        Float64 AStepRatio = 0.15;
                        Int32 maxIterations = 6;
                        Int32 iterations = 0;
                        Float64 ACorrected=ampsA[i2];

                        Float64 ValidAcorrected = ampsA[i2];
                        Float64 ValidEcorrected = ampsE[i2];

                        Float64 overlap = 0.33;
                        std::vector<UInt32> indexesFitted;
                        std::vector<UInt32> EOverlapIdx = LineModelElementList.getOverlappingElements( iEltE[i2], indexesFitted, overlap );
                        std::vector<UInt32> EOverlapIdxA = LineModelElementList.getOverlappingElements( iEltA[i2], indexesFitted, overlap );
                        for(Int32 io=0; io<EOverlapIdxA.size(); io++)
                        {
                            EOverlapIdx.push_back(EOverlapIdxA[io]);
                        }
                        LineModelElementList.refreshModelUnderElements( EOverlapIdx );
                        Float64 previousFitErr = LineModelElementList.getModelErrorUnderElement( iEltE[i2] );
                        while(iterations<maxIterations)
                        {
                            iterations++;
                            ACorrected=ACorrected*(1.0+AStepRatio);
                            Float64 correction = ACorrected-ampsA[i2];
                            LineModelElementList.m_Elements[iEltA[i2]]->SetFittedAmplitude(ACorrected, ersA[i2]);
                            LineModelElementList.m_Elements[iEltE[i2]]->SetFittedAmplitude(ampsE[i2]+correction+ACorrected*widthRatioAE, ersE[i2]);
                            LineModelElementList.refreshModelUnderElements(EOverlapIdx);
                            Float64 newFitErr = LineModelElementList.getModelErrorUnderElement( iEltE[i2] );
                            if(newFitErr >= previousFitErr)
                            { //todo put a ratio threshold ?
                                iterations = maxIterations;
                            }
                            else
                            {
                                ValidAcorrected = ACorrected;
                                ValidEcorrected = ampsE[i2]+correction;
                            }
                        }
                        ampsA[i2] = ValidAcorrected;
                        ampsE[i2] = ValidEcorrected+ValidAcorrected*widthRatioAE;
                        LineModelElementList.m_Elements[iEltA[i2]]->SetFittedAmplitude(ValidAcorrected, ersA[i2]);
                        LineModelElementList.m_Elements[iEltE[i2]]->SetFittedAmplitude(ValidEcorrected+ValidAcorrected*widthRatioAE, ersE[i2]);
                    }
                }
            }
        }
    }
    /*
    TFloat64List coeffs = BalmerModelLinSolve( lambdax, continuumx, datax, errdatax );
    // export for debug
    FILE* fspc = fopen( "BalmerLinSolve_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1.0;
    for(Int32 i=0; i<nLines; i++){
        Float64 ampRegE = coeffs[0]*lambdax[i] + coeffs[3] ;
        Float64 coeffRegA = (coeffs[2]*lambdax[i] + coeffs[1]);
        Float64 coeffA = (continuumx[i] - ampsA[i])/continuumx[i];
        fprintf( fspc, "%d %f %f %f %f %f %f\n", i, lambdax[i], errdatax[i], ampRegE, ampsE[i], coeffRegA, coeffA);
    }
    fclose( fspc );
    //*/
    return;
}

/**
 * \brief Make a linear fit of the lines selected for Balmer rule correction.
 **/
TFloat64List CRuleBalmerLinearSolver::BalmerModelLinSolve( std::vector<Float64> lambdax, std::vector<Float64> continuumx, std::vector<Float64> datax, std::vector<Float64> errdatax )
{
  //Linear fit
  int i, n;
  Float64 fval;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  n = lambdax.size();
  Int32 nddl = 4;
  if(n<nddl)
    {
      TFloat64List empty;
      return empty;
    }

  X = gsl_matrix_alloc (n, nddl);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);

  c = gsl_vector_alloc (nddl);
  cov = gsl_matrix_alloc (nddl, nddl);

  for (i = 0; i < n; i++)
    {
      double yi, ei;
      yi = datax[i];
      ei = errdatax[i];

      for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
	  if( iddl==0 )
	    {
	      fval = lambdax[i];
	    }
	  if( iddl==1 )
	    {
	      fval = continuumx[i];
	    }
	  if( iddl==2 )
	    {
	      fval = lambdax[i]*continuumx[i];
	    }
	  if( iddl==3 )
	    {
	      fval = 1.0;
	    }
	  gsl_matrix_set (X, i, iddl, fval);
        }
      gsl_vector_set (y, i, yi);
      gsl_vector_set (w, i, 1.0/(ei*ei));
    }

  {
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free (work);
  }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
  if(1)
    {
      Log.LogDebug( "# best fit: Y = %g X + %g CX + %g XCX + %g", C(0), C(1), C(2), C(3) );
      Log.LogDebug( "# covariance matrix:\n" );
      Log.LogDebug( "[ %+.5e, %+.5e \n", COV(0,0), COV(0,1) );
      Log.LogDebug( "  %+.5e, %+.5e \n", COV(1,0), COV(1,1) );
      Log.LogDebug( "[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2) );
      Log.LogDebug( "  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2) );
      Log.LogDebug( "  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2) );
      Log.LogDebug( "# chisq/n = %g", chisq/n );
    }

  TFloat64List coeffs;
  coeffs.push_back(gsl_vector_get(c,(0)));
  coeffs.push_back(gsl_vector_get(c,(1)));
  coeffs.push_back(gsl_vector_get(c,(2)));
  coeffs.push_back(gsl_vector_get(c,(3)));

  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (w);
  gsl_vector_free (c);
  gsl_matrix_free (cov);

  return coeffs;
}


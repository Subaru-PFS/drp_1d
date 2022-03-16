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
#include "RedshiftLibrary/spectrum/LSFVariableWidth.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/exception.h"
using namespace NSEpic;
using namespace std;

//used to create mapping between width and lambda values
//Instead of this 
CLSFGaussianVariableWidth::CLSFGaussianVariableWidth(const std::shared_ptr<const TLSFGaussianVarWidthArgs>& args):
    CLSF(GaussianVariableWidth, std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())),
    m_width(args->width),
    m_spcAxis(args->lambdas)
{
    IsValid();
}

Float64 CLSFGaussianVariableWidth::GetWidth(Float64 lambda) const
{
    Int32 idx = -1;
    if(!checkAvailability(lambda))
    {
        throw GlobalException(INTERNAL_ERROR,"CLSFGaussianVariableWidth::GetWidth: lambda outside spectralAxis range");
    }

    TFloat64Index::getClosestLowerIndex(m_spcAxis.GetSamplesVector(), lambda, idx);

    if(m_spcAxis[idx]!=lambda)
    {//interpolation
        Float64 t = (lambda - m_spcAxis[idx])/(m_spcAxis[idx+1] - m_spcAxis[idx]);
        Float64 width_interp = m_width[idx] + (m_width[idx+1] - m_width[idx])*t;
        return width_interp;
    }
    return m_width[idx];
}

bool CLSFGaussianVariableWidth::IsValid() const
{
    if(!m_width.size()){
        throw GlobalException(BAD_COUNTMATCH, "CLSFGaussianVariableWidth::Width array cannot be null ");
    }
    if(m_spcAxis.GetSamplesCount() != m_width.size()) {
        throw GlobalException(BAD_COUNTMATCH, "CLSFGaussianVariableWidth::isValid Spectral axis size and width axis size do not match ");
    }
    for(Float64 w : m_width)
        if(w <= 0. ) return false;
    return true;
}

bool CLSFGaussianVariableWidth::checkAvailability(Float64 lambda)const
{   
    bool available = true;
    if(lambda<m_spcAxis[0] || lambda>m_spcAxis[m_width.size()-1])
        available = false;
    return available;
}
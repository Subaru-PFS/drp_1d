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
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"

#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
using namespace NSEpic;


std::shared_ptr<const COperatorResult> TplCombinationExtremaResult::getCandidate(const int& rank,const std::string& dataset) const
{ 
  if (dataset == "model_parameters" || dataset == "modeltplCombination_parameters")  return std::make_shared<const TTplCombinationResult>(this->m_ranked_candidates[rank].second);
  else if (dataset == "model")  return this->m_savedModelSpectrumResults[rank];
  else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
}
    
const std::string& TplCombinationExtremaResult::getCandidateDatasetType(const std::string& dataset) const 
{
      if (dataset == "model_parameters" || dataset == "modeltplCombination_parameters")      return this->m_ranked_candidates[0].second.getType();
      else if (dataset == "model")  return this->m_savedModelSpectrumResults[0]->getType();
      else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
}

bool TplCombinationExtremaResult::HasCandidateDataset(const std::string& dataset) const
{
  return (dataset == "model_parameters" || dataset == "model" || dataset == "modeltplCombination_parameters");
}

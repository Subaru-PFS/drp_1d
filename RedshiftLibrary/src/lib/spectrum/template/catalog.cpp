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
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/io/genericreader.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/continuum/median.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/continuum/waveletsdf.h"

#include "RedshiftLibrary/log/log.h"

#include <boost/filesystem.hpp>

#include <string>

using namespace NSEpic;
using namespace std;
using namespace boost::filesystem;


/**
 * Variable instantiator constructor.
 */
CTemplateCatalog::CTemplateCatalog( std::string cremovalmethod, Float64 mediankernelsize, Float64 waveletsScales, std::string waveletsDFBinPath, Bool sampling )
{
    m_continuumRemovalMethod = cremovalmethod;
    m_continuumRemovalMedianKernelWidth = mediankernelsize;
    m_continuumRemovalWaveletsNScales = waveletsScales;
    m_continuumRemovalWaveletsBinPath = waveletsDFBinPath;
    m_logsampling = sampling;
}


TTemplateConstRefList CTemplateCatalog::const_TTemplateRefList_cast(const TTemplateRefList & list)
{
    TTemplateConstRefList const_list;
    for (auto tpl : list) const_list.push_back(tpl);

    return const_list;
}

/**
 * Returns a list containing all templates as enumerated in the categoryList input.
 */
TTemplateRefList CTemplateCatalog::GetTemplateList_( const TStringList& categoryList ) const
{
    TTemplateRefList list;

    for( Int32 i=0; i<categoryList.size(); i++ )
    {
        for ( Int32 j=0; j<GetTemplateCount( categoryList[i] ); j++ )
        {
            list.push_back( GetList().at( categoryList[i] )[j] );
        }
    }

    return list;
}

std::shared_ptr<const CTemplate>  CTemplateCatalog::GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const
{
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        for( UInt32 j=0; j<GetTemplateCount( tplCategoryList[i] ); j++ )
        {
            std::shared_ptr<const CTemplate> tpl = GetTemplate( tplCategoryList[i], j );
            if(tpl->GetName() == tplName){
                return tpl;
            }
        }
    }
    throw std::runtime_error("Could not find template with name");
}


/**
 * Get a list of strings with the contents of m_List.
 */
TStringList CTemplateCatalog::GetCategoryList() const
{
    TStringList l;
    for (auto it : GetList()){
        l.push_back( it.first );
    }
    return l;
}

/**
 * Returns the size of the category entry in m_List.
 */
UInt32 CTemplateCatalog::GetTemplateCount( const std::string& category ) const
{   
    UInt32 l; 
    auto ret = GetList().find( category );

    if( GetList().find( category ) == GetList().end() )
        return 0;
    l = GetList().at( category ).size();

    return l;
}

/**
 * Adds the input to the list of templates, under its category. If the input doesn't have a category, function returns false. Also computes the template without continuum and adds it to the list of templates without continuum. Returns true.
 * @sampling here is relevant here especially for when rebinning
 * otherwise we will have to SetSampling ("lin") when we want to read the non-log template, and then call SetSampling("log")
 * when we want to add the log tpl
 */
void CTemplateCatalog::Add( const std::shared_ptr<CTemplate> & r)
{
    if( r->GetCategory().empty() )
      throw runtime_error("Template has no category");
        
    GetList()[r->GetCategory()].push_back( r );
}

void CTemplateCatalog::SetTemplate( const std::shared_ptr<CTemplate> & tpl, UInt32 i)
{
    GetList().at(tpl->GetCategory())[i] = tpl;
}

//adapt it to apply to all m_list
void CTemplateCatalog::InitIsmIgm(const std::string & calibrationPath, 
                                  std::shared_ptr<const CParameterStore> parameterStore,
                                  const std::shared_ptr<const CLSF>& lsf)
{
    Float64 ebmv_start=0.0;
    Float64 ebmv_step=0.1;
    UInt32 ebmv_n=10;
    parameterStore->Get( "ebmv.start", ebmv_start, 0. );
    parameterStore->Get( "ebmv.step", ebmv_step, 0.1 );
    parameterStore->Get( "ebmv.count", ebmv_n, 10 );
    //ISM
    auto ismCorrectionCalzetti = std::make_shared<CSpectrumFluxCorrectionCalzetti>();
    ismCorrectionCalzetti->Init(calibrationPath, ebmv_start, ebmv_step, ebmv_n);
    //IGM
    auto igmCorrectionMeiksin = std::make_shared<CSpectrumFluxCorrectionMeiksin>();
    igmCorrectionMeiksin->Init(calibrationPath, lsf);

    //push in all templates
    //backup current sampling
    Bool currentsampling = m_logsampling;
    for(Bool sampling : {0, 1}){
        m_logsampling = sampling;
        for(auto it : GetList())
        {             
            const TTemplateRefList  & TplList = it.second;
            for (auto tpl : TplList)
                tpl->m_ismCorrectionCalzetti = ismCorrectionCalzetti;
            if(it.first != "star")//no igm for stars
                for (auto tpl : TplList)
                    tpl->m_igmCorrectionMeiksin = igmCorrectionMeiksin;
        }
    }
    //put back the initial sampling:
    m_logsampling = currentsampling;
}
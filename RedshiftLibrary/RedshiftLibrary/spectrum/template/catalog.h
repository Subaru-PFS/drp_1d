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
#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_
#define _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#include <boost/filesystem.hpp>
#include <vector>

namespace NSEpic
{

class CTemplateCatalog
{

public:

    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8.0, std::string waveletsDFBinPath="", Bool sampling = 0 );
    
    void                    Add( const std::shared_ptr<CTemplate> & tpl);

    std::shared_ptr<const CTemplate>        GetTemplate( const std::string& category, UInt32 i ) const;
    std::shared_ptr<const CTemplate>        GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const;
    void    SetTemplate( const std::shared_ptr<CTemplate> & tpl,  UInt32 i);
    void    ClearTemplates(const std::string & category,  Bool opt_ortho, Bool opt_logsampling, UInt32 i, Bool alltemplates=false);
    void    ClearTemplateList(const std::string & category,  Bool opt_ortho, Bool opt_logsampling);

    TTemplateConstRefList GetTemplateList( const TStringList& categoryList ) const;
    TTemplateRefList GetTemplateList( const TStringList& categoryList );
    
    static TTemplateConstRefList const_TTemplateRefList_cast(const TTemplateRefList & list);
    TStringList             GetCategoryList() const;
    UInt32                  GetTemplateCount( const std::string& category ) const;
    UInt32                  GetNonNullTemplateCount( const std::string& category ) const;
    void                    InitIsmIgm(std::shared_ptr<const CParameterStore> parameterStore,
                                       const std::shared_ptr<const CLSF>& lsf);
    mutable Bool            m_logsampling = 0;//non-log by default
    mutable Bool            m_orthogonal = 0;//non-log by default
    Float64                 m_ortho_LSFWidth = NAN;
private:
    // this const version must stay private, since it returns non const templates.
    TTemplateRefList GetTemplateList_( const TStringList& categoryList ) const; 

    UInt32 GetNonNullTemplateCount( const std::string& category, Bool opt_ortho, Bool opt_logsampling ) const;

          TTemplatesRefDict &    GetList();
    const TTemplatesRefDict &    GetList() const;

          TTemplatesRefDict &    GetList( Bool opt_ortho, Bool opt_logsampling);
    const TTemplatesRefDict &    GetList( Bool opt_ortho, Bool opt_logsampling) const;

    typedef std::vector<std::vector<TTemplatesRefDict>> TTemplatesRefDictAA;
    TTemplatesRefDictAA m_ListMatrix{2, std::vector<TTemplatesRefDict>(2)};//row corresponds to original vs ortho; col corresponds to orig vs rebinned

    std::string m_continuumRemovalMethod;
    Float64 m_continuumRemovalMedianKernelWidth;
    Float64 m_continuumRemovalWaveletsNScales;
    std::string m_continuumRemovalWaveletsBinPath;

};

/**
 * Returns the contents of the i-th entry in the category item of m_List.
 */
inline 
std::shared_ptr<const CTemplate> CTemplateCatalog::GetTemplate( const std::string& category, UInt32 i ) const
{
    return GetList().at( category )[i];   
}

// non const getter returning mutable templates
inline 
TTemplateRefList CTemplateCatalog::GetTemplateList( const TStringList& categoryList )
{
    return GetTemplateList_(categoryList);
}

//  const getter returning const templates
inline
TTemplateConstRefList CTemplateCatalog::GetTemplateList( const TStringList& categoryList ) const
{
    return const_TTemplateRefList_cast( GetTemplateList_(categoryList)); 
}


 //below functions aim at avoid using if..else to access the right categoryList
inline 
const TTemplatesRefDict & CTemplateCatalog::GetList( Bool opt_ortho, Bool opt_logsampling) const
{
    return const_cast<CTemplateCatalog*>(this)->GetList(opt_ortho, opt_logsampling);
}

inline 
const TTemplatesRefDict & CTemplateCatalog::GetList() const
{
    return const_cast<CTemplateCatalog*>(this)->GetList();
}

inline 
TTemplatesRefDict & CTemplateCatalog::GetList( Bool opt_ortho, Bool opt_logsampling)
{
    return m_ListMatrix[opt_ortho][opt_logsampling];
}

inline 
TTemplatesRefDict & CTemplateCatalog::GetList()
{
    return m_ListMatrix[m_orthogonal][m_logsampling];
}

inline
UInt32 CTemplateCatalog::GetNonNullTemplateCount( const std::string& category) const
{
    return GetNonNullTemplateCount(category, m_orthogonal, m_logsampling);
}

}

#endif

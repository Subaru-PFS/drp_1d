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
#include "RedshiftLibrary/processflow/parameterstore.h"

namespace bpt = boost::property_tree;

namespace NSEpic
{
CParameterStore::CParameterStore(const TScopeStack& stack):
  CScopeStore(stack)
{

}

void CParameterStore::Set( const std::string& name, const TFloat64List& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.put_child( name, array );
}

void CParameterStore::Set( const std::string& name, const TStringList& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.put_child( name, array );
}

void CParameterStore::Set( const std::string& name, const TFloat64Range& v )
{
    TFloat64List list( 2 );

    list[0] = v.GetBegin();
    list[1] = v.GetEnd();

    Set( name, list );
}

void CParameterStore::Set( const std::string& name, const TInt64List& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.put_child( name, array );
}

void CParameterStore::Set( const std::string& name, const TBoolList& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.put_child( name, array );
}

void CParameterStore::Set( const std::string& name, Float64 v )
{
    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( name );

    m_PropertyTree.put( name, v );
}

void CParameterStore::Set( const std::string& name, Int64 v )
{
    boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>( name );

    m_PropertyTree.put( name, v );
}

void CParameterStore::Set( const std::string& name, bool v )
{
    boost::optional<bool> property = m_PropertyTree.get_optional<bool>( name );

    m_PropertyTree.put( name, v );
}

void CParameterStore::Set( const std::string& name, const std::string& v )
{
    boost::optional<std::string> property = m_PropertyTree.get_optional<std::string>( name );

    m_PropertyTree.put( name, v );
}


void CParameterStore::Save( const std::string& path ) const
{
    bpt::json_parser::write_json( path, m_PropertyTree );
}

void CParameterStore::FromString(const std::string& json)
{
    std::istringstream jsonstream(json);
    bpt::json_parser::read_json(jsonstream, m_PropertyTree);
}

bool CParameterStore::HasFFTProcessing(const std::string &objectType) const
{
    bool fft_processing = false;

    if(Has<bool>(objectType + ".templatefittingsolve.fftprocessing"))
        fft_processing |= Get<bool>(objectType + ".templatefittingsolve.fftprocessing");
    if(Has<bool>(objectType + ".linemodelsolve.linemodel.continuumfit.fftprocessing"))
        fft_processing |= Get<bool>(objectType + ".linemodelsolve.linemodel.continuumfit.fftprocessing");

    return fft_processing;
}

bool CParameterStore::HasToOrthogonalizeTemplates(const std::string &objectType) const
{  
    bool orthogonalize = Get<std::string>(objectType + ".method") == "linemodelsolve";
    if(orthogonalize){
        std::string continuumComponent = Get<std::string>(objectType + ".linemodelsolve.linemodel.continuumcomponent");
        orthogonalize &= (continuumComponent == "tplfit" || continuumComponent == "tplfitauto" );
    }
    return orthogonalize;
}
bool CParameterStore::EnableTemplateOrthogonalization(const std::string &objectType) const
{  
    bool enableOrtho = HasToOrthogonalizeTemplates(objectType);
    if(enableOrtho)
    {
        enableOrtho &= !Get<bool>(objectType + ".linemodelsolve.linemodel.continuumfit.ignorelinesupport");
    }
    return enableOrtho;
}

}

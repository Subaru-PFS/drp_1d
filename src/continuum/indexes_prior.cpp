
#include <epic/redshift/continuum/indexes_prior.h>

#include <iostream>

#include <epic/core/common/datatypes.h>
#include <epic/core/log/log.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>


namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

/**
 * Get the prior value for a gievn color, break pair
 */
Float64 CContinuumIndexesPrior::GetHeatmapVal( Int32 _index, Float64 _color, Float64 _break)
{
    //deal with nans
    if(_color!=_color)
    {
        return 1.0;
    }
    if(_break!=_break)
    {
        return 1.0;
    }

    Float64 color_max = m_tbl_color_min+m_tbl_color_n*m_tbl_color_step;

    if(_color>=color_max || _color<m_tbl_color_min)
    {
        return 0.0;
    }
    Float64 break_max = m_tbl_break_min+m_tbl_break_n*m_tbl_break_step;

    if(_break>=break_max || _break<m_tbl_break_min)
    {
        return 0.0;
    }
    Int32 icolor = (Int32)((_color-m_tbl_color_min)/m_tbl_color_step);
    Int32 ibreak = (Int32)((_break-m_tbl_break_min)/m_tbl_break_step);
    Float64 val = m_ciprior_table[_index][icolor][ibreak];
    return val/m_ciprior_max[_index];
}

/**
 * Init
 */
bool CContinuumIndexesPrior::Init( std::string calibrationPath )
{
    bfs::path calibrationFolder( calibrationPath.c_str() );
    //std::string ciprior_path = (calibrationFolder.append( "continuu_indexes_prior_map_20161130" )).string();
    //std::string ciprior_path = (calibrationFolder.append( "continuum_indexes_prior_map_201612_eucsimuHaOii" )).string();
    std::string ciprior_path = (calibrationFolder/"continuum_indexes_prior_map_201612_eucsimuHaOii").string();


    //std::string cfg_basename = "heatmap_ycolor_xbreak__blurred_fixedgrid_"; unused, harcoded cfg for now
    m_tbl_color_step = 0.1;
    m_tbl_color_min = -3.0;
    m_tbl_color_n = 50;
    m_tbl_break_step = 0.1;
    m_tbl_break_min = -3.0;
    m_tbl_break_n = 50;

    std::string dat_basename = "heatmap_ycolor_xbreak__blurred_fixedgrid_";
    m_nIndexes = 6;

    for(Int32 iIndex=0; iIndex<m_nIndexes; iIndex++)
    {
        TContinuumIndexData _contIndexData;
        Float64 maxVal = 0.0;
        std::string filename = boost::str(boost::format("%s%d.dat") %dat_basename %iIndex);
        //std::string filePath = bfs::path( ciprior_path ).append( filename).string();
        std::string filePath = (bfs::path( ciprior_path )/filename).string();
        std::ifstream file;
        file.open( filePath, std::ifstream::in );
        bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
        if(fileOpenFailed)
        {
            Log.LogError("Continuum indexes prior table, unable to load the calib. file... aborting!");
            return false;
        }

        std::string line;
        // Read file line by line
        Int32 kLine = 0;
        while( getline( file, line ) )
        {
            if( !boost::starts_with( line, "#" ) )
            {
                std::istringstream csvStream( line );
                std::vector<Float64> csvColumn;
                string csvElement;
                while( getline(csvStream, csvElement, '\t') )
                {
                    Float64 val = (Float64)std::atof(csvElement.c_str());
                    if(maxVal<val)
                    {
                        maxVal = val;
                    }
                    csvColumn.push_back(val);
                }
                if(csvColumn.size()!=m_tbl_break_n)
                {
                    Log.LogError("Continuum indexes prior table, not enough colums in the calib. file... aborting!");
                    return false;
                }

                _contIndexData.push_back(csvColumn);
                kLine++;
                if(kLine>=m_tbl_color_n)
                {
                    break;
                }
            }
        }
        file.close();

        m_ciprior_max.push_back(maxVal);
        m_ciprior_table.push_back(_contIndexData);
    }

    return true;
}


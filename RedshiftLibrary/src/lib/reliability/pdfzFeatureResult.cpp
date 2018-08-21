#include <RedshiftLibrary/reliability/pdfzFeatureResult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace NSEpic;


CPdfzFeatureResult::CPdfzFeatureResult()
{

}

CPdfzFeatureResult::~CPdfzFeatureResult()
{

}


void CPdfzFeatureResult::Save( const CDataStore& store, std::ostream& stream ) const
{
	std::setprecision(20);
    stream << "#zPDF_descriptors \tValue" << std::endl;

	boost::unordered_map<const std::string , Float64>::const_iterator it = mapzfeatures.begin();
	for ( it=mapzfeatures.begin(); it!=mapzfeatures.end(); ++it) {
		std::string map_key = (it->first);
		Float64 map_content = (it->second);
        //std::cout
        //         << " \t map_key: " << map_key << "\n"
        //        << " \t map_content: " << map_content << "\n";

		stream << map_key << "\t" << map_content << std::endl;
	}

}


void CPdfzFeatureResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
	//stream.precision(20);
	//stream << "CPdfzFeatureResult" << "\t" << mapzfeatures.size() << std::endl;

}

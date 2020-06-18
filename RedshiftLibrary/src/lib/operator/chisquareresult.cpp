#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string/predicate.hpp>

using namespace NSEpic;

CChisquareResult::CChisquareResult()
{

}

CChisquareResult::~CChisquareResult()
{

}

void CChisquareResult::Init(UInt32 n , Int32 nISM, Int32 nIGM)
{
    ChiSquare.resize( n );
    FitAmplitude.resize( n );
    FitAmplitudeError.resize( n );
    FitAmplitudeNegative.resize( n);
    FitDustCoeff.resize( n );
    FitMeiksinIdx.resize( n );
    FitDtM.resize( n );
    FitMtM.resize( n );
    LogPrior.resize( n );
    Redshifts.resize( n );
    Overlap.resize( n );
    Status.resize( n );

    ChiSquareIntermediate.clear();
    IsmDustCoeffIntermediate.clear();
    IgmMeiksinIdxIntermediate.clear();
    for(Int32 k=0; k<n; k++)
    {
        std::vector<TFloat64List> _ChiSquareISMList;
        std::vector<TFloat64List> _IsmDustCoeffISMList;
        std::vector<TInt32List> _IgmMeiksinIdxISMList;
        for(Int32 kism=0; kism<nISM; kism++)
        {

            TFloat64List _chi2IGMList(nIGM, DBL_MAX);
            _ChiSquareISMList.push_back(_chi2IGMList);

            TFloat64List _dustCoeffIGMList(nIGM, -1.0);
            _IsmDustCoeffISMList.push_back(_dustCoeffIGMList);

            TInt32List _meiksinIdxIGMList(nIGM, -1);
            _IgmMeiksinIdxISMList.push_back(_meiksinIdxIGMList);

        }
        ChiSquareIntermediate.push_back(_ChiSquareISMList);
        IsmDustCoeffIntermediate.push_back(_IsmDustCoeffISMList);
        IgmMeiksinIdxIntermediate.push_back(_IgmMeiksinIdxISMList);
    }
}

void CChisquareResult::Load( std::istream& stream )
{
    // Clear current lines list
    Redshifts.clear();
    ChiSquare.clear();
    Overlap.clear();
    Status.clear();

    std::string line;

    // Read file line by line
    while( std::getline( stream, line ) )
    {
        if( !boost::starts_with( line, "#" ) )
        {
            boost::char_separator<char> sep(" \t");

            // Tokenize each line
            typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
            ttokenizer tok( line, sep );

            // Check if it's not a comment
            ttokenizer::iterator it = tok.begin();
            if( it != tok.end() )
            {
                // Parse redshift
                Float64 r = -1.0;
                try
                {
                    r = boost::lexical_cast<Float64>(*it);
                }
                catch (boost::bad_lexical_cast&)
                {
                    return;
                }

                // Parse merit
                ++it;
                Float64 c = DBL_MAX;
                if( it != tok.end() )
                {
                    try
                    {
                        c = boost::lexical_cast<Float64>(*it);
                    }
                    catch (boost::bad_lexical_cast&)
                    {
                        return;
                    }
                }
                else
                {
                    return;
                }

                // Parse overlap
                ++it;
                Float64 o = -1.0;
                if( it != tok.end() )
                {
                    try
                    {
                        o = boost::lexical_cast<Float64>(*it);
                    }
                    catch (boost::bad_lexical_cast&)
                    {
                        o = -1.0;
                    }
                }
                else
                {
                    o = -1.0;
                }

                Redshifts.push_back( r );
                ChiSquare.push_back( c );
                Overlap.push_back( o );
                Status.push_back( COperator::nStatus_OK );
            }
        }
    }
}

void CChisquareResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream <<  "#Redshifts\tChiSquare\tOverlap"<< std::endl;
    for ( int i=0; i<Redshifts.size(); i++)
    {
        stream <<  Redshifts[i] << std::setprecision(16) << "\t" << std::scientific << ChiSquare[i] << std::fixed << "\t" << Overlap[i] << std::endl;
    }

    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;

        if(ChiSquare.size()>0){
            stream <<  "#ExtremaMerit = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(16) << "\t" << std::scientific <<   ChiSquare[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitAmplitude.size()>0){
            stream <<  "#Extrema FitAmplitudes = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(16) << "\t" << std::scientific <<   FitAmplitude[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitAmplitudeError.size()>0){
            stream <<  "#Extrema FitAmplitudeErrors = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(16) << "\t" << std::scientific <<   FitAmplitudeError[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitAmplitudeNegative.size()>0){
            stream <<  "#Extrema FitAmplitudesNegative = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << "\t" <<   FitAmplitudeNegative[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitDustCoeff.size()>0){
            stream <<  "#Extrema FitDustCoeff = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << std::setprecision(4) << "\t" <<   FitDustCoeff[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(FitMeiksinIdx.size()>0){
            stream <<  "#Extrema FitMeiksinIdx = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << "\t" <<   FitMeiksinIdx[idx] << "\t";
            }
            stream << "}" << std::endl;
        }

        if(LogPrior.size()>0){
            stream <<  "#Extrema LogPrior = {";
            for ( int i=0; i<Extrema.size(); i++)
            {
                Int32 idx=0;
                for ( UInt32 i2=0; i2<Redshifts.size(); i2++)
                {
                    if(Redshifts[i2] == Extrema[i]){
                        idx = i2;
                        break;
                    }
                }

                stream << "\t" <<   LogPrior[idx] << "\t";
            }
            stream << "}" << std::endl;
        }
    }
}

void CChisquareResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "ChisquareResult" << "\t" << Redshifts.size() << std::endl;
}

bool CChisquareResult::CallFind() {
  Int32 extremumCount = 10;
  if (Redshifts.size() > extremumCount) {
    TPointList extremumList;
    TFloat64Range redshiftsRange( Redshifts[0], Redshifts[Redshifts.size() - 1]);
    CExtremum extremum(redshiftsRange, extremumCount, true);
    extremum.Find(Redshifts, ChiSquare, extremumList);
    // Refine Extremum with a second maximum search around the z candidates:
    Float64 radius = 0.001;
    for (Int32 i = 0; i < extremumList.size(); i++) {
      Float64 x = extremumList[i].X;
      Float64 left_border = std::max(redshiftsRange.GetBegin(), x - radius);
      Float64 right_border = std::min(redshiftsRange.GetEnd(), x + radius);

      TPointList extremumListFine;
      TFloat64Range rangeFine = TFloat64Range(left_border, right_border);
      CExtremum extremumFine(rangeFine, 1, true);
      extremumFine.Find(Redshifts, ChiSquare,
        extremumListFine);
      if (extremumListFine.size() > 0) {
        extremumList[i] = extremumListFine[0];
      }
    }
    // store extrema results
    Extrema.resize(extremumCount);
    for (Int32 i = 0; i < extremumList.size(); i++) {

      Extrema[i] = extremumList[i].X;
    }

  } else {
    // store extrema results
    Extrema.resize(Redshifts.size());
    TFloat64List tmpX;
    TFloat64List tmpY;
    for (Int32 i = 0; i < Redshifts.size(); i++) {
      tmpX.push_back(Redshifts[i]);
      tmpY.push_back(ChiSquare[i]);
    }
    // sort the results by merit
    NSEpic::CQuickSort < Float64 > sort;
    vector < Int32 > sortedIndexes(Redshifts.size());
    sort.SortIndexes(tmpY.data(), sortedIndexes.data(), sortedIndexes.size());
    for (Int32 i = 0; i < Redshifts.size(); i++) {
      Extrema[i] = tmpX[sortedIndexes[i]];
    }
  }
  return true;
}
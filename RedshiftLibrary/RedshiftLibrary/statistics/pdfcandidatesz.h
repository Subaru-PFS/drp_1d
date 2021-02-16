#ifndef _REDSHIFT_STATISTICS_PDFCANDIDATESZ_
#define _REDSHIFT_STATISTICS_PDFCANDIDATESZ_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/common/range.h>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <utility>

namespace NSEpic
{

struct TCandidateZ {
    Float64           	  Redshift = NAN;
    Float64               ValProba = NAN;
    Float64           	  ValSumProba = 0.;
    Float64               Deltaz = 0.;
    Float64           	  ValSumProbaZmin = NAN;
    Float64           	  ValSumProbaZmax = NAN;

    //opt 1: direct integration
    //
    //opt 2: gaussian fit
    Float64           		GaussAmp = NAN;
    Float64           		GaussAmpErr = NAN;
    Float64           		GaussSigma = NAN;
    Float64           		GaussSigmaErr = NAN;

};
/*
struct TCandidateZ: TCandidateZ 
{  
  std::string Id;

  TCandidateZRanked(const TCandidateZ & _candidate, 
                     cosnt std::string> & _Id):
      TCandidateZ(_candidate),
      Id(_Id)
  {}
*/


typedef std::map<std::string, TCandidateZ> TCandidateZbyID;
typedef std::vector<std::pair<std::string, TCandidateZ>> TCandidateZbyRank;
typedef std::map<std::string, TFloat64Range> TCandidateZRangebyID;

class CPdfCandidateszResult;

class CPdfCandidatesZ : public COperator
{

public:

    CPdfCandidatesZ(const TCandidateZbyID & candidates);
    CPdfCandidatesZ(const TFloat64List & redshifts);
   
    std::shared_ptr<CPdfCandidateszResult> Compute(TRedshiftList const & PdfRedshifts, 
                                  TFloat64List const & PdfProbaLog);
    
    TStringList SetIntegrationWindows(const TFloat64Range PdfZRange, TCandidateZRangebyID & ranges);

    Int32     m_optMethod = 0; //0: direct integration, 1:gaussian fit
    Float64   m_dzDefault = 1e-3; // default value in case deltaz couldnt be calculted, should be instrument dependant (parameter ?)

    TCandidateZbyID m_candidates;
  
private:

    Bool getCandidateSumTrapez(const TRedshiftList &redshifts,
                                  const TFloat64List &valprobalog,
                                  const TFloat64Range  &zrange,
                                  TCandidateZ & candidate) const;//default: zwidth_left = zwidth_right

    Bool getCandidateRobustGaussFit(const TRedshiftList &redshifts,
                                    const TFloat64List &valprobalog,
                                    const TFloat64Range &zrange,
                                    TCandidateZ & candidate) const;
    
    Bool getCandidateGaussFit(const TRedshiftList &redshifts,
                              const TFloat64List &valprobalog,
                              const TFloat64Range &zrange,
                              TCandidateZ & candidate) const;

    void SortByValSumProbaInt(TCandidateZbyRank & ranked_candidates) const; 
};



}

#endif

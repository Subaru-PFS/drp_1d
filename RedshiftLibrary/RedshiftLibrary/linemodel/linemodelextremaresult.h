#ifndef _REDSHIFT_OPERATOR_LINEMODELEXTREMARESULT_
#define _REDSHIFT_OPERATOR_LINEMODELEXTREMARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <RedshiftLibrary/linemodel/linemodelsolution.h>

namespace NSEpic
{

class CLineModelExtremaResult
{

public:

    CLineModelExtremaResult();
    ~CLineModelExtremaResult();

    Void Resize(Int32 size);

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    //Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    //Void Load( std::istream& stream );

    //Extrema results
    TFloat64List            Extrema;    // z extrema
    TFloat64List            ExtremaMerit;    // extrema merit
    TFloat64List            DeltaZ;    // extrema delta z
    TFloat64List            mTransposeM;    // extrema model norm
    TFloat64List            ExtremaLastPass; //z extrema with more precision
    TFloat64List            lmfitPass;// z found with lmfit

    //Deprecated?
    TFloat64List            ExtremaExtendedRedshifts;    // z range around extrema
    TFloat64List            Posterior;    // z extrema
    TFloat64List            LogArea;    // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;    //corrected z for each extrema
    TFloat64List            SigmaZ; //sigmaz for each extrema

    //
    TFloat64List            StrongELSNR;
    TFloat64List            bic;    // bayesian information criterion for each extrema
    std::vector<CContinuumIndexes::TContinuumIndexList> ContinuumIndexes; //continuum indexes for each extrema
    std::vector<CMask>      OutsideLinesMask;   //Mask with 0 under the lines and 1 anywhere else
    std::vector<std::string>      FittedTplName;    //Name of the best template fitted for continuum
    TFloat64List            FittedTplAmplitude;     //Amplitude for the best template fitted for continuum
    TFloat64List            FittedTplDustCoeff;     //Calzetti dustcoeff for the best template fitted for continuum
    std::vector<Int32>      FittedTplMeiksinIdx;    //Meiksin igm index for the best template fitted for continuum
    std::vector<std::string>      FittedTplshapeName;   //Name of the best template fitted for tplcorr/tplshape

};


}

#endif

class CPdfMargZLogResult : public COperatorResult
{

  public:
  CPdfMargZLogResult();
    ~CPdfMargZLogResult() = default;
    CPdfMargZLogResult(const TFloat64List & redshifts);

    Int32 getIndex( Float64 z ) const;

    TFloat64List          Redshifts;
    TFloat64List          valProbaLog;
    Float64               valEvidenceLog;
    UInt32 				  countTPL;
};

class TCandidateZ : public COperatorResult{
public:

   TCandidateZ(const TCandidateZ&) = default;
    TCandidateZ(TCandidateZ&&) = default;
    TCandidateZ& operator=(const TCandidateZ&) = default;
    TCandidateZ& operator=(TCandidateZ&&) = default;
    virtual ~TCandidateZ() = default;
    TCandidateZ() = default;
    
    Float64           	  Redshift = NAN;
    Float64               ValProba = NAN;
    Float64           	  ValSumProba = 0.;
    Float64               Deltaz = 0.;
    Float64           	  ValSumProbaZmin = NAN;
    Float64           	  ValSumProbaZmax = NAN;
               
    std::string           ParentId = "";

    //opt 1: direct integration
    //
    //opt 2: gaussian fit
    Float64           		GaussAmp = NAN;
    Float64           		GaussAmpErr = NAN;
    Float64           		GaussSigma = NAN;
    Float64           		GaussSigmaErr = NAN;

};

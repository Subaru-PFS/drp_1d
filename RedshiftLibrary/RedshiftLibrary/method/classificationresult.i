class CClassificationResult : public CSolveResult
{

public:

    CClassificationResult();

    void SetTypeLabel( std::string lbl );
    void SetG(Float64 evidence, Float64 prob);
    void SetS(Float64 evidence, Float64 prob);
    void SetQ(Float64 evidence, Float64 prob);


    std::string m_TypeLabel="-1";
    Float64 m_evidence_galaxy=-1.0;
    Float64 m_evidence_star=-1.0;
    Float64 m_evidence_qso=-1.0;
    Float64 m_prob_galaxy=-1.0;
    Float64 m_prob_star=-1.0;
    Float64 m_prob_qso=-1.0;
};

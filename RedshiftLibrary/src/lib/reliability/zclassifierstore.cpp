#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>
#include <gsl/gsl_matrix.h>
#include <vector>

#include <RedshiftLibrary/reliability/zclassifierstore.h>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <string>
#include <fstream>
#include <iostream>

//#include <boost/tokenizer.hpp>
//#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace bfs = boost::filesystem;
namespace bft = boost::this_thread;

using namespace boost;
using namespace std;
using namespace NSEpic;


/*  ---------------------------------------------------------------------
 * 					>>CLASS C_CLASSIFIER_STORE
/*  --------------------------------------------------------------------- */
CClassifierStore::CClassifierStore()
{
    m_isInitialized = false;
    m_classifier_option = 2; //1: svm gauss, 2:svm lin
    m_file_format_version = -1;
}

CClassifierStore::~CClassifierStore()
{

}


/* ---------------------------------------------------------------------
 * 					>>CLASS C_LEARNER
/* --------------------------------------------------------------------- */

CClassifierStore::CLearner::CLearner()
{

}

CClassifierStore::CLearner::~CLearner()
{

}

/* ---------------------------------------------------------------------
 * 					>>	LOAD - PARAMETERS FILES (.dat or .txt)
/* --------------------------------------------------------------------- */
Bool CClassifierStore::Load ( const char* directoryPath )
{
	Bool disp_time = false;
	boost::posix_time::ptime  startTime;
	boost::posix_time::time_duration diff;

	startTime = boost::posix_time::microsec_clock::local_time();

	SetTypeClassifier("SVM" );
	SetTypeCoding( "ECOC_OVA" );


    Log.LogDetail("  ZClassifier: Loading with classifier option: %d", m_classifier_option);
    Log.LogInfo("  ZClassifier: Loading from zclassifier directory: %s", directoryPath);

    m_file_format_version = Load_version( directoryPath );
    Log.LogDetail("  ZClassifier: Found file format version = %.3f", m_file_format_version);

    Bool ret = Load_params( directoryPath );

	if ( ret )
	{
		FILE* f;
		Int32 rows, cols;
		bfs::path dirPath;

		// LOAD CODING MATRIX
        Log.LogDetail("  ZClassifier: Loading codingMatrix");
		dirPath = directoryPath;
		dirPath=dirPath/"zClassifier_codingMatrix.dat";
		rows = GetNbClasses();
		cols = GetNbLearners();
		gsl_matrix* codMat = gsl_matrix_alloc(rows,cols);
		f = fopen( dirPath.c_str(), "r");
		if( f==NULL)	{
            Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
			return false;
		}
		gsl_matrix_fscanf( f, codMat );

		SetCodingMatrix(codMat);
		SetCodingMatrixPos();
		SetCodingMatrixNeg();
        //DisplayQ(GetCodingMatrix());


		// LOAD SIGMOID PARAMS & BIAS
        Log.LogDetail("  ZClassifier: Loading learnersParams");
		dirPath = directoryPath;
		dirPath=dirPath/"zClassifier_learnersParams.dat";
        if ( !boost::filesystem::exists( dirPath ) )
        {
            dirPath = directoryPath;
            dirPath=dirPath/"zClassifier_learnersParams.csv";
            if ( !boost::filesystem::exists( dirPath ) )
            {
                Log.LogError("  ZClassifier: ERROR File does not exist: %s", dirPath.string().c_str());
            }
        }

		cols = 3; // SLOPE ; INTERCEPT; BIAS
		rows=GetNbLearners();
		params_L = gsl_matrix_alloc(rows, cols);
		f = fopen( dirPath.c_str(), "r");
		if ( f==NULL) 	{
            Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
			return false;
		}
		gsl_matrix_fscanf( f, params_L );

		// LOAD LEARNERS
		Int32 key_learner;
		std::string name, nbLearner;
		gsl_matrix* sv_params;
		gsl_matrix* sv_vectors;
		for (Int32 i = 0; i<GetNbLearners(); i++)
		{
            key_learner = i;
            nbLearner = boost::lexical_cast<std::string>(key_learner+1); //= std::to_string(key_learner);


            if(m_classifier_option==1)
            {
                // LOAD VECTORS
                Log.LogDetail("  ZClassifier: Loading vectors #%d", i);
                dirPath = directoryPath;
                name = "sv"+nbLearner+"/sv"+nbLearner+"_vectors.dat";
                dirPath=dirPath/name;
                if ( !boost::filesystem::exists( dirPath ) )
                {
                    dirPath = directoryPath;
                    name = "sv"+nbLearner+"/sv"+nbLearner+"_vectors.csv";
                    dirPath=dirPath/name;
                    if ( !boost::filesystem::exists( dirPath ) )
                    {
                        Log.LogError("  ZClassifier: ERROR File does not exist: %s", dirPath.string().c_str());
                    }
                }
                rows = temp_sv[i];
                cols = GetNbFeatures()+ 1+ 1; // the matrix contains the vector Alpha, the vector SVlabels and the matrix SVectors
                Log.LogDetail("  ZClassifier: rows:%d cols:%d", rows, cols);
                sv_vectors = gsl_matrix_alloc(rows,cols);
                f = fopen( dirPath.c_str(), "r");
                if ( f==NULL){
                    Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
                    return false;
                }
                gsl_matrix_fscanf( f, sv_vectors );
            }

            if(m_classifier_option==2)
            {
                // LOAD PARAMS
                Log.LogDetail("  ZClassifier: Loading params #%d", i);
                dirPath = directoryPath;
                name = "sv"+nbLearner+"/sv"+nbLearner+"_params.dat";
                dirPath=dirPath/name;
                if ( !boost::filesystem::exists( dirPath ) )
                {
                    dirPath = directoryPath;
                    name = "sv"+nbLearner+"/sv"+nbLearner+"_params.csv";
                    dirPath=dirPath/name;
                    if ( !boost::filesystem::exists( dirPath ) )
                    {
                        Log.LogError("  ZClassifier: ERROR File does not exist: %s", dirPath.string().c_str());
                    }
                }
                rows = GetNbFeatures();
                cols = 3;
                sv_params = gsl_matrix_alloc(rows,cols);
                f = fopen( dirPath.c_str(), "r");
                if ( f==NULL) {
                    Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
                    return false;
                }
                gsl_matrix_fscanf( f, sv_params );
            }


            Log.LogDetail("  ZClassifier: creating learner #%d", i);
			auto learner = new CLearner();
			learner->m_SVmu = gsl_vector_alloc(GetNbFeatures());
            learner->m_SVsigma = gsl_vector_alloc(GetNbFeatures());
			learner->m_SVsigmoiid.resize(2);
            if(m_classifier_option==2)
            {
                learner->m_SVbeta = gsl_vector_alloc(GetNbFeatures());
                learner->LoadParams (sv_params);
            }
            if(m_classifier_option==1)
            {
                learner->m_SValpha = gsl_vector_alloc(temp_sv[i]);
                learner->m_SVectorLabels = gsl_vector_alloc(temp_sv[i]);
                learner->m_SVectors = gsl_matrix_alloc(temp_sv[i],GetNbFeatures());

                learner->LoadVectors ( sv_vectors);

                learner->m_nbSVectors = temp_sv[i];
            }
			learner->m_nbDesriptors = GetNbFeatures() ;
			learner->m_SVbias = gsl_matrix_get( params_L, key_learner, 2);

			learner->m_SVsigmoiid[0] = gsl_matrix_get( params_L, key_learner, 0); // SLOPE
			learner->m_SVsigmoiid[1] = gsl_matrix_get( params_L, key_learner, 1); // INTERCEPT

			m_learners[key_learner + 1] = std::shared_ptr<CLearner>(learner);  // keys in MAP_LEARNERS start at 1

		}

		// FREE ALLOX_MEMORY
		gsl_matrix_free ( codMat );
        if(m_classifier_option==1)
        {
            gsl_matrix_free ( sv_vectors );
        }
        if(m_classifier_option==2)
        {
            gsl_matrix_free ( sv_params );
        }

        gsl_matrix_free( params_L );
		fclose(f);

		// CHECK COMPUTATIONAL TIME
		diff = boost::posix_time::microsec_clock::local_time() - startTime;
		if (disp_time) {std::cout <<"\t>> [ZQUAL] load : " << to_simple_string(diff) <<std::endl;}


		/*
		Float64 bias;
		std::cout << "############### boost uMAP_Learners - acces 1 ###############" << std::endl;
		MapLearners::const_iterator it=(GetLearners()).begin();
		for (it=(GetLearners()).begin();it!=(GetLearners()).end();++it) {
			bias = (it->second)->m_SVbias;
			std::cout << "m_learners  "<< (it->first)<<" ; bias = " << bias << std::endl;
		}
		std::cout << "############### boost uMAP_Learners - acces 2 ###############" << std::endl;
		it=(m_learners).begin();
		for (it=(m_learners).begin();it!=(m_learners).end();++it) {
			auto Lrn = std::dynamic_pointer_cast<CLearner>( (it->second) );
			bias = Lrn->m_SVbias;
			std::cout << "m_learners  "<< (it->first)<<" ; bias = " << bias << std::endl;
		}
		std::cout << "############### boost uMAP_Learners - acces 3 ###############" << std::endl;
		for (Int32 j=0; j<GetNbLearners(); j++) {
			Int32 key_id = j+1;
			std::cout << "m_learners "<<key_id<<" ;  bias "<< m_learners[key_id]->m_SVbias << std::endl;
		}
		 */
        m_isInitialized = true;
        Log.LogDetail("  ZClassifier: Successfully initialized");
    }else{
        Log.LogError("  ZClassifier: Unable to load the learner params from directory %s", directoryPath);
        m_isInitialized = false;
        return false;
    }


	return true;
}

/**
 * @brief CClassifierStore::Load_version
 * Load version number of the zclassifier dir calibration data
 * output <1 : load using SJamal file format
 * output = 1.0 : new format designed with M. Gray
 * @param directoryPath
 * @return
 */
Float64 CClassifierStore::Load_version ( const char* directoryPath)
{
    Float64 ver = 0.0;
    bfs::path filePath;
    filePath = directoryPath;
    filePath=filePath/"zClassifier_version.csv";
    if ( !boost::filesystem::exists( filePath ) )
    {
        return 0.0;
    }
    ifstream file;
    file.open( filePath.string(), ifstream::in );
    if( file.rdstate() & ios_base::failbit )
    {
        return -1;
    }

    std::string line;

    // Read file line by line
    while( std::getline( file, line ) )
    {
        if( boost::starts_with( line, "#" ) )
        {
            continue;
        }

        ver = boost::lexical_cast<Float64>(line);
        break;
    }
    return ver;
}

Bool CClassifierStore::Load_params ( const char* directoryPath)
{
	FILE* f;
    Int32 rows;
    bfs::path dirPath;
    gsl_vector* params_global_info;
    gsl_vector* params;


	// LOAD GLOBAL PARAMS
	dirPath = directoryPath;
    dirPath=dirPath/"zClassifier_params_ova.dat";
    if ( !boost::filesystem::exists( dirPath ) )
    {
        dirPath = directoryPath;
        dirPath=dirPath/"zClassifier_params_ova.csv";
        if ( !boost::filesystem::exists( dirPath ) )
        {
            Log.LogError("  ZClassifier: ERROR File does not exist: %s", dirPath.string().c_str());
        }
    }
    //first read 3 first values
    rows=3;
    params_global_info = gsl_vector_alloc(rows);
    f = fopen( dirPath.c_str(), "r");
    if ( f==NULL) 	{
        Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
        return false;
    }
    gsl_vector_fscanf( f, params_global_info );

    // UPDATE GLOBAL INFO
    Int32 nbFeat = (Int32)  params_global_info->data[0];
    Log.LogDetail("  ZClassifier: Nb Features read = %d", nbFeat);
    SetNbFeatures( nbFeat );
    Int32 nbClasses = (Int32) params_global_info->data[1];
    Log.LogDetail("  ZClassifier: Nb Classes read = %d", nbClasses);
    SetNbClasses( nbClasses );
    Int32 nbLearners = (Int32) params_global_info->data[2];
    Log.LogDetail("  ZClassifier: Nb Learners read = %d", nbLearners);
    SetNbLearners( nbLearners );
    m_Labels.resize(GetNbClasses());
    for (Int32 j = 0; j< GetNbClasses(); j++) {
        m_Labels[j] = "C"+std::to_string(j+1);
    }

    if(m_file_format_version<1.0) //read nums separated by EOL
    {
        rows = nbClasses+3+1; //ROWS: nbFeatures; nbClass; nbLearners; size each SVectors ; NAN (to dismiss)
    }else{
        rows = nbClasses+3; //COLS: nbFeatures; nbClass; nbLearners; size each SVectors ;
    }
	params = gsl_vector_alloc(rows);
    f = fopen( dirPath.c_str(), "r");
    if ( f==NULL) 	{
        Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
        return false;
    }
    gsl_vector_fscanf( f, params );

    // STORE TEMP INFO
    temp_sv.resize(GetNbClasses());
    for (Int32 j = 0; j< GetNbClasses(); j++) {
        temp_sv[j] = params->data[3+j];
    }

	// LOAD LEARNERS_WEIGHT
	dirPath = directoryPath;
    dirPath=dirPath/"zClassifier_learnersWeight.dat";
    if ( !boost::filesystem::exists( dirPath ) )
    {
        dirPath = directoryPath;
        dirPath=dirPath/"zClassifier_learnersWeight.csv";
        if ( !boost::filesystem::exists( dirPath ) )
        {
            Log.LogError("  ZClassifier: ERROR File does not exist: %s", dirPath.string().c_str());
        }
    }
    if(m_file_format_version<1.0) //read nums separated by EOL
    {
        rows = GetNbLearners()+1;
    }else{
        rows = GetNbLearners();
    }

	params = gsl_vector_alloc(rows);
	f = fopen( dirPath.c_str(), "r");
	if ( f==NULL) 	{
        Log.LogError("  ZClassifier: ERROR in reading %s", dirPath.string().c_str());
		return false;
	}
	gsl_vector_fscanf( f, params );

	// UPDATE WEIGHTS
	m_learnersWeight = gsl_vector_alloc(GetNbLearners());
	SetLearnerWeight(params);

	// FREE ALLOC_MEMORY
	gsl_vector_free( params );
	fclose(f);

	return true;

}



Bool CClassifierStore::CLearner::LoadParams ( gsl_matrix* params)
{
    //gsl_matrix_get_col (gsl_vector * v, const gsl_matrix * m, size_t j) 	:	copies the elements of the j-th column of the matrix m into the vector v.
    gsl_matrix_get_col ( m_SVmu, params, 0);
    gsl_matrix_get_col ( m_SVsigma, params, 1);
    gsl_matrix_get_col ( m_SVbeta, params, 2);
    return true;

}
Bool CClassifierStore::CLearner::LoadVectors ( gsl_matrix* vectors)
{

    // gsl_matrix_set_col (gsl_matrix * m, size_t j, const gsl_vector * v)	:  copies the elements of the vector v into the j-th column of the matrix m.
    gsl_matrix_get_col ( m_SValpha, vectors, 0);
    gsl_matrix_get_col ( m_SVectorLabels, vectors, 1);

    gsl_vector* temp_v = gsl_vector_alloc(vectors->size1);
    for ( Int32 col = 0; col<vectors->size2-2; col++ )
    {
        gsl_matrix_get_col ( temp_v, vectors, col+2);
        gsl_matrix_set_col ( m_SVectors, col, temp_v);
    }


    gsl_vector_free(temp_v);

    return true;

}


void CClassifierStore::DisplayQ ( const gsl_matrix* m )
{
	std::cout << "---------------------------------------------------------------------------------------------------------"<<"\n"
			<<">> ["<<GetTypeCoding()<<"] Coding Matrix :  \t\t"
			<< m->size1 << " classes (rows) x "
			<< m->size2 << " learners (columns) "<<std ::endl;
	for (Int32 i = 0; i<m->size1; i++){
		std::cout <<"\t Class ["<<(i+1)<<"] \t";
		for (Int32 j = 0; j<m->size2; j++){
			std::string val;
			if (gsl_matrix_get( m, i, j) == 1) {val = "+1";}
			else {val = std::to_string((int)gsl_matrix_get( m, i, j));}
			std::cout << val << " ";
		}
		std::cout <<"\n";
	}

}



/* ---------------------------------------------------------------------
 * 					>>	GET - METHODS
/* --------------------------------------------------------------------- */
const CClassifierStore::MapLearners& CClassifierStore::GetLearners() const
{
	return m_learners;
}

const Int32 CClassifierStore::GetNbClasses() const
{
	return m_nbClasses;
}

const Int32 CClassifierStore::GetNbLearners() const
{
	return m_nbLearners;
}

const Int32 CClassifierStore::GetNbFeatures() const
{
	return m_nbFeatures;
}

const std::string CClassifierStore::GetTypeCoding() const
{
	return m_typeCoding;
}

const std::string CClassifierStore::GetTypeClassifier() const
{
    return m_typeClassifier;
}

const Int32 CClassifierStore::GetOptionClassifier() const
{
    return m_classifier_option;
}

const std::string CClassifierStore::GetLabel( Int32 idLabel ) const
{
	return m_Labels[idLabel];
}

const TStringList CClassifierStore::GetLabels( ) const
{
	return m_Labels;
}

const gsl_vector* CClassifierStore::GetLearnersWeight() const
{
	return m_learnersWeight;
}

const gsl_matrix* CClassifierStore::GetCodingMatrix() const
{
	return m_codingMatrix;
}

const gsl_matrix* CClassifierStore::GetCodingMatrixPos() const
{
	return m_codingMatrixPos;
}

const gsl_matrix* CClassifierStore::GetCodingMatrixNeg() const
{
	return m_codingMatrixNeg;
}

/* ---------------------------------------------------------------------
 * 					>>	SET - METHODS
/* --------------------------------------------------------------------- */

void CClassifierStore::SetLearners ( CClassifierStore::MapLearners& learners )
{
	m_learners = learners;
}

void CClassifierStore::SetNbClasses( Int32 nbclasses )
{
	m_nbClasses = (Int32) nbclasses;
}

void CClassifierStore::SetNbLearners( Int32 nblearners )
{
	m_nbLearners = (Int32) nblearners;
}

void CClassifierStore::SetNbFeatures( Int32 nbfeatures )
{
	m_nbFeatures = (Int32) nbfeatures;
}

void CClassifierStore::SetTypeCoding( std::string typecoding )
{
	m_typeCoding = typecoding;
}

void CClassifierStore::SetTypeClassifier( std::string typeclassifier )
{
    m_typeClassifier = typeclassifier;
}

void CClassifierStore::SetLearnerWeight( gsl_vector* w )
{
	for (Int32 i = 0; i<w->size-1; i++ ){
		m_learnersWeight->data[i] = w->data[i];
	}
}

void CClassifierStore::SetCodingMatrix( gsl_matrix* m )
{
	m_codingMatrix = gsl_matrix_alloc(GetNbClasses(), GetNbLearners());
	gsl_matrix_memcpy(m_codingMatrix, m);
}

void CClassifierStore::SetCodingMatrixPos( )
{
	m_codingMatrixPos = gsl_matrix_alloc( GetNbClasses(), GetNbLearners() );
	gsl_matrix_memcpy( m_codingMatrixPos, GetCodingMatrix() );
	for ( Int32 i = 0; i< m_codingMatrixPos->size1; i++) {
		for ( Int32 j = 0; j< m_codingMatrixPos->size2; j++) {
			if ( gsl_matrix_get(m_codingMatrixPos, i, j ) == -1) {
				gsl_matrix_set(m_codingMatrixPos, i, j, 0);
			}
		}
	}
}

void CClassifierStore::SetCodingMatrixNeg( )
{
	m_codingMatrixNeg = gsl_matrix_alloc( GetNbClasses(), GetNbLearners());
	gsl_matrix_memcpy( m_codingMatrixNeg, GetCodingMatrix() );
	for ( Int32 i = 0; i< m_codingMatrixNeg->size1; i++) {
		for ( Int32 j = 0; j< m_codingMatrixNeg->size2; j++) {
			if ( gsl_matrix_get(m_codingMatrixNeg, i, j ) == 1) {
				gsl_matrix_set(m_codingMatrixNeg, i, j, 0);
			}
		}
	}
}


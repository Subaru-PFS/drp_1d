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
#ifndef _REDSHIFT_EXCEPTION_
#define _REDSHIFT_EXCEPTION_

#include <exception>
#include <iostream>
#include <string>
//#include <boost/stacktrace.hpp>

namespace NSEpic {
  typedef enum ErrorCode
    {
      INTERNAL_ERROR=0,
      EXTERNAL_LIB_ERROR,
      INVALID_SPECTRA_FLUX,
      INVALID_NOISE	,
      SMALL_WAVELENGTH_RANGE ,
      NEGATIVE_CONTINUUMFIT	,
      BAD_CONTINUUMFIT	,
      NULL_AMPLITUDES	,
      PEAK_NOT_FOUND_PDF	,
      MAX_AT_BORDER_PDF	,
      UNKNOWN_PARAMETER  ,
      BAD_PARAMETER_VALUE,
      UNKNOWN_ATTRIBUTE ,
      BAD_LINECATALOG,
      BAD_LOGSAMPLEDSPECTRUM,
      BAD_COUNTMATCH,
      BAD_TEMPLATECATALOG,
      INVALID_SPECTRUM,
      OVERLAPRATE_NOTACCEPTABLE,
      DZ_NOT_COMPUTABLE
    } ErrorCode;

  class AmzException : public std::exception 
  {
  protected:
    std::string _msg;
    std::string stacktrace;
    ErrorCode code;
  public:

    AmzException(ErrorCode ec,std::string message);

    AmzException(const AmzException& e);
    
    virtual ~AmzException();
    virtual const char* what() const noexcept override
    {
      return _msg.c_str();
    }

    const std::string & getMessage()
    {
      return _msg;
    }

    ErrorCode getErrorCode(){return code;}
    void setErrorCode(ErrorCode ec){code = ec;}

    const char* getStackTrace() const ;
  };


  // A solve exception stops the whole pipeline
  // This exception should be caught only from pylibamazed or a client
  class GlobalException: public AmzException
  {
  public:
    GlobalException(ErrorCode ec,std::string message):
      AmzException(ec,message)
    {

    }
    GlobalException(const GlobalException& e):
      AmzException(e){
    }

    virtual ~GlobalException();

  };

  // A solve exception stops a solve method (which is considered as failed), but not the whole pipeline
  // This exception should be caught only from pylibamazed or a client
  class SolveException: public AmzException
  {
  public:
    SolveException(ErrorCode ec,std::string message):
      AmzException(ec,message)
    {      
    }
    SolveException(const SolveException& e):
      AmzException(e){}

    virtual ~SolveException();
  };

    // This exception is for now reserved for the unknown parameter exception
  // It should be handled at initialization only, not the case as parameters are retrieved in the process flow
  // This exception class should also be used for bad parameters values (e.g. a negative value when a positive value is required)
  class ParameterException: public AmzException
  {
  public:
    ParameterException(ErrorCode ec,std::string message):
      AmzException(ec,message)
    {
    }
    ParameterException(const ParameterException& e):
      AmzException(e){}

    virtual ~ParameterException();
    const std::string & getMessage()
    {
      return _msg;
    }
  };

  //only this exception class should be caught from the library
  class InternalException: public AmzException
  {
  public:
    InternalException(ErrorCode ec,std::string message):
      AmzException(ec,message)
    {      
    }
    InternalException(const InternalException& e):
      AmzException(e){}

    virtual ~InternalException();

  };

  
 
}

#endif

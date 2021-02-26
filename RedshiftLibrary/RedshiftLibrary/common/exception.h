#ifndef _EXCEPTION_
#define _EXCEPTION_

#include <exception>
#include <iostream>
#include <string>
//#include <boost/stacktrace.hpp>

namespace NSEpic {
  typedef enum ErrorCode
    {
      INVALID_SPECTRA_FLUX=	0,
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
      BAD_LINECATALOG
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
    
    ~AmzException(){}
    virtual const char* what() const noexcept override
    {
      return _msg.c_str();
    }

    void getMessage(std::string &msg)
    {
      msg = _msg;
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

    ~GlobalException(){}

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

    ~SolveException(){}
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

    ~ParameterException(){}
    void getMessage(std::string &msg)
    {
      msg = _msg;
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

    ~InternalException(){}

  };

  
 
}

#endif

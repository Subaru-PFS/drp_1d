#ifndef _EXCEPTION_
#define _EXCEPTION_

#include <exception>
#include <iostream>
#include <string>

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
  class Exception : public std::exception 
  {
  protected:
    std::string _msg;
    ErrorCode code;
  public:

    Exception(ErrorCode ec,std::string message):
      _msg(message),
      code(ec)
    {

    }

    Exception(const Exception& e):
      _msg(e._msg),
      code(e.code){}
      
    ~Exception(){}
    virtual const char* what() const noexcept override
    {
      return _msg.c_str();
    }

    void getMessage(std::string &msg)
    {
      msg = _msg;
    }

    ErrorCode getErrorCode(){return code;}
  };


  // A solve exception stops the whole pipeline
  // This exception should be caught only from pylibamazed or a client
  class GlobalException: public Exception
  {
  public:
    GlobalException(ErrorCode ec,std::string message):
      Exception(ec,message)
    {
    }
    GlobalException(const Exception& e):
      Exception(e){}

    ~GlobalException(){}

  };

  // A solve exception stops a solve method (which is considered as failed), but not the whole pipeline
  // This exception should be caught only from pylibamazed or a client
  class SolveException: public Exception
  {
  public:
    SolveException(ErrorCode ec,std::string message):
      Exception(ec,message)
    {
      
    }
    SolveException(const Exception& e):
      Exception(e){}

    ~SolveException(){}
  };

    // This exception is for now reserved for the unknown parameter exception
  // It should be handled at initialization only, not the case as parameters are retrieved in the process flow
  // This exception class should also be used for bad parameters values (e.g. a negative value when a positive value is required)
  class ParameterException: public Exception
  {
  public:
    ParameterException(ErrorCode ec,std::string message):
      Exception(ec,message)
    {
    }
    ParameterException(const Exception& e):
      Exception(e){}

    ~ParameterException(){}
    void getMessage(std::string &msg)
    {
      msg = _msg;
    }
  };

  //only this exception class should be caught from the library
  class InternalException: public Exception
  {
  public:
    InternalException(ErrorCode ec,std::string message):
      Exception(ec,message)
    {
      
    }
    ~InternalException(){}

  };

  
 
}

#endif

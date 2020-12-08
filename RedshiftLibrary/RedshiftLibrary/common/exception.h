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
    /*
    template<typename ... Args>
    Exception(ErrorCode ec,std::string& format, Args&& ... args) {
      size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
      //      if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
      std::unique_ptr<char[]> buf( new char[ size ] ); 
      snprintf( buf.get(), size, format.c_str(), args ... );
      _msg = std::string( buf.get(), buf.get() + size - 1 );
    }

    template<typename ... Args>
    Exception(ErrorCode ec,const char* format, Args&& ... args) {
      size_t size = snprintf( nullptr, 0, format, args ... ) + 1; // Extra space for '\0'
      //      if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
      std::unique_ptr<char[]> buf( new char[ size ] ); 
      snprintf( buf.get(), size, format, args ... );
      _msg = std::string( buf.get(), buf.get() + size - 1 );
    }
    //  Exception(const std::string& msg) : _msg(msg){}
  */
    virtual const char* what() const noexcept override
    {
      return _msg.c_str();
    }  
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
  };

  //only this exception class should be caught from the library
  class InternalException: public Exception
  {
  public:
    InternalException(ErrorCode ec,std::string message):
      Exception(ec,message)
    {
      
    }
  };

  
 
}

#endif

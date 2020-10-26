#ifndef _EXCEPTION_
#define _EXCEPTION_

#include <exception>
#include <iostream>
#include <string>

namespace NSEpic {

  class Exception : public std::exception 
  {
    std::string _msg;

  public:
    template<typename ... Args>
      Exception(std::string& format, Args&& ... args) {
      size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
      //      if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
      std::unique_ptr<char[]> buf( new char[ size ] ); 
      snprintf( buf.get(), size, format.c_str(), args ... );
      _msg = std::string( buf.get(), buf.get() + size - 1 );
    }

    template<typename ... Args>
      Exception(const char* format, Args&& ... args) {
      size_t size = snprintf( nullptr, 0, format, args ... ) + 1; // Extra space for '\0'
      //      if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
      std::unique_ptr<char[]> buf( new char[ size ] ); 
      snprintf( buf.get(), size, format, args ... );
      _msg = std::string( buf.get(), buf.get() + size - 1 );
    }
      //  Exception(const std::string& msg) : _msg(msg){}

    virtual const char* what() const noexcept override
    {
      return _msg.c_str();
    }
  
  };
 
}

#endif

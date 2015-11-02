#include <epic/core/log/filehandler.h>

#include <epic/core/log/handler.h>

#include <stdio.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT_NOT_INSTANCIABLE( CLogFileHandler )

CLogFileHandler::CLogFileHandler( CLog& logger, const char* filePath ) :
    CLogHandler( logger )
{
    m_OutputStream.open( filePath, std::fstream::out );
}

CLogFileHandler::~CLogFileHandler()
{
    m_OutputStream.close();
}

Void CLogFileHandler::LogEntry( UInt32 lvl, const char* header, const char* msg )
{
    if( header )
    {
        m_OutputStream << header;
    }

    m_OutputStream << msg << std::endl;

}

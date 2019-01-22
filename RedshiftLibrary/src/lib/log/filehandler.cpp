#include <RedshiftLibrary/log/filehandler.h>

#include <RedshiftLibrary/log/handler.h>

#include <stdio.h>

using namespace NSEpic;

CLogFileHandler::CLogFileHandler( CLog& logger, const char* filePath ) :
    CLogHandler( logger )
{
    m_OutputStream.open( filePath, std::fstream::out );
}

CLogFileHandler::~CLogFileHandler()
{
    m_OutputStream.close();
}

void CLogFileHandler::LogEntry( UInt32 lvl, const char* header, const char* msg )
{
    if( header )
    {
        m_OutputStream << header;
    }

    m_OutputStream << msg << std::endl;

}

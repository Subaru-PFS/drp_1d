#include <RedshiftLibrary/log/handler.h>
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;


CLogHandler::CLogHandler( CLog& logger )
{
    m_Logger = &logger;
    m_LevelMask = CLog::nLevel_Warning;
    m_Logger->AddHandler( *this );
}

CLogHandler::~CLogHandler()
{
    m_Logger->RemoveHandler( *this );
}

void CLogHandler::SetLevelMask( UInt32 mask )
{
    m_LevelMask = mask;
}

UInt32 CLogHandler::GetLevelMask() const
{
    return m_LevelMask;
}

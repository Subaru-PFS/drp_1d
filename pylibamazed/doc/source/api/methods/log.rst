Log
===

.. contents::
   :local:
   :depth: 2
   :backlinks: entry

Logging
-------

.. currentmodule:: pylibamazed.redshift

Class ``CLog`` can be used to log informations.
To use it, first get its instance ( it is a singleton), then use the desired method to log at the desired log level.

.. autosummary:: 

  CLog
  CLog.GetInstance

To log information, use:

.. autosummary:: 

  CLog.LogInfo
  CLog.LogDetail
  CLog.LogDebug

Its available log levels are : 

.. autosummary:: 
  
  CLog.nLevel_Critical
  CLog.nLevel_Error
  CLog.nLevel_Warning
  CLog.nLevel_Info
  CLog.nLevel_Detail
  CLog.nLevel_Debug
  CLog.nLevel_None


Minimal example :

.. code:: python

    from pylibamazed.redshift import CLog

    zlog = CLog.GetInstance()
    zlog.LogDebug("Here is the error message")


Handling logs
-------------

There are two classes depending on the way you want to handle logs:

.. currentmodule:: pylibamazed.redshift

.. autosummary:: 

  CLogConsoleHandler
  CLogFileHandler


Both inherits from a ``CLogHandler`` base class and have the following methods :

.. autosummary:: 

  CLogConsoleHandler.SetLevelMask

.. autosummary:: 
  
  CLogFileHandler.SetLevelMask

Minimal example:

.. code:: python

    from pylibamazed.redshift import CLog, CLogConsoleHandler

    logConsoleHandler = CLogConsoleHandler()
    logConsoleHandler.SetLevelMask(CLog.nLevel_Info)

Reference
=========

.. currentmodule:: pylibamazed.redshift

.. autofunction:: get_version() -> str
  
  Gets lib version


AbstractExternalStorage
-----------------------

.. currentmodule::  pylibamazed.AbstractExternalStorage

.. autoclass:: AbstractExternalStorage
  :members: 
  :undoc-members:
  :show-inheritance:


.. automodule:: pylibamazed.AbstractExternalStorage
   :members:


AbstractSpectrumReader
----------------------

.. currentmodule::  pylibamazed.AbstractSpectrumReader

.. autoclass:: AbstractSpectrumReader
  :members: 
  :undoc-members:
  :show-inheritance:


.. automodule:: pylibamazed.AbstractSpectrumReader
   :members:

BuilderPdfHandler
-----------------
.. currentmodule::  pylibamazed.PdfHandler

.. autoclass:: BuilderPdfHandler
  :members: 
  :undoc-members:
  :show-inheritance:

CalibrationLibrary
------------------

.. currentmodule:: pylibamazed.CalibrationLibrary

.. autoclass:: CalibrationLibrary
  :members: 
  :undoc-members:
  :show-inheritance:

CLog
----

.. currentmodule:: pylibamazed.redshift

.. autoclass:: CLog
   :members: nLevel_Critical, nLevel_Error, nLevel_Warning, nLevel_Info, nLevel_Detail, nLevel_Debug, nLevel_None
   :undoc-members:
   :show-inheritance:

   .. automethod:: LogInfo(message: str)

    Log an information at info level

   .. automethod:: LogDetail(message: str)

    Log an information at detail level

   .. automethod:: LogDebug(message: str)

    Log an information at debug level
  
   .. automethod:: GetInstance()

    Access the ``CLog`` instance

CLogConsoleHandler
------------------

.. currentmodule:: pylibamazed.redshift

.. autoclass:: CLogConsoleHandler
   :members:
   :undoc-members:

   Write logs to the console
    
   .. automethod:: SetLevelMask
     
     For log levels below mask, log will not be handled.

CLogFileHandler
---------------

.. currentmodule:: pylibamazed.redshift

.. autoclass:: CLogFileHandler
   :members:
   :undoc-members:

   Write logs to a file. Must be initialized with a file name.

   .. automethod:: SetLevelMask
     
     For log levels below mask, log will not be handled.

Container
---------

.. currentmodule:: pylibamazed.Container

.. autoclass:: Container
  :members:
  :undoc-members:
  :show-inheritance:

ErrorCode
---------

Enum containing all regirstered error codes. See :doc:`/api/errorswarnings/errors`. 

H5Writer
--------

.. currentmodule:: pylibamazed.H5Writer

.. autoclass:: H5Writer
  :members:
  :undoc-members:
  :show-inheritance:

Parameters
----------

.. currentmodule:: pylibamazed.Parameters

.. autoclass:: Parameters
  :members:
  :undoc-members:
  :show-inheritance:
  :inherited-members:

PdfHandler
----------

.. currentmodule:: pylibamazed.PdfHandler

.. autoclass:: PdfHandler
  :members:
  :undoc-members:
  :show-inheritance:

.. automodule:: pylibamazed.PdfHandler
  :members:


ProcessFlow
-----------
   
.. currentmodule:: pylibamazed.ProcessFlow

.. autoclass:: ProcessFlow
  :members:
  :undoc-members:
  :show-inheritance:

ResultStoreOutput
-----------------

.. currentmodule:: pylibamazed.ResultStoreOutput

.. autoclass:: ResultStoreOutput
  :members:
  :undoc-members:
  :show-inheritance:
  :inherited-members:

Spectrum
---------

.. currentmodule:: pylibamazed.Spectrum

.. autoclass:: Spectrum
  :members:
  :undoc-members:
  :show-inheritance:


WarningCode
-----------

Enum containing all regirstered warning codes. See :doc:`/api/errorswarnings/warnings`. 
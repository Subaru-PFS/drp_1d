Running the ProcessFlow
========================

ProcessFlow
-----------

.. currentmodule:: pylibamazed.ProcessFlow

The ``ProcessFlow`` class is used to launch the processing of a spetrum.

.. autosummary:: 
  ProcessFlow

The processing is launched with:

.. autosummary:: 
  ProcessFlow.run
     
Running the ProcessFlow consists in two steps :

1. Building the ProcessFlow :

   .. code:: python

    process_flow = ProcessFlow(
      {
        "calibration_dir":"absolute_path_to_calib_directory",
        "extended_results":True
      },
      Parameters(params_dict)
    )

2. Running the processFlow

   .. code:: python

     result = process_flow.run()

Its public attributes are:

.. autosummary:: 
  ProcessFlow.parameters
  ProcessFlow.calibration_library


CalibrationLibrary
--------------------

``CalibrationLibrary`` allows to access the calibration data needed for the process flow.

.. currentmodule:: pylibamazed.CalibrationLibrary

Its public attributes are:

.. autosummary:: 
  CalibrationLibrary.line_catalogs_df
  CalibrationLibrary.parameters

ResultStoreOutput
------------------

.. currentmodule:: pylibamazed.ResultStoreOutput

``ResultStoreOutput`` is the output of the ``ProcessFlow.run()`` method.

.. autosummary:: 
  ResultStoreOutput

Its public methods are:

.. autosummary:: 
  ResultStoreOutput.has_error
  ResultStoreOutput.get_error
  ResultStoreOutput.has_attribute
  ResultStoreOutput.get_attribute
  ResultStoreOutput.get_candidate_data
  ResultStoreOutput.get_nb_candidates
  ResultStoreOutput.get_dataset

And its public attributes are:

.. autosummary::
  ResultStoreOutput.parameters
  ResultStoreOutput.object_results
Parameters
==========

The parameters class
--------------------

The ``Parameters`` class is used to store the parameters of the amazed pipeline. It is initialized with a dictionary structured as defined in :doc:`/json-schema/general`.

.. currentmodule::  pylibamazed.Parameters

.. autosummary::
  Parameters

``Parameters`` contains several getters:

.. autosummary::
    Parameters.get_spectrum_models
    Parameters.get_solve_methods_str
    Parameters.stage_enabled
    Parameters.get_redshift_solver_method
    Parameters.get_linemeas_method
    Parameters.get_json_schema_version
    Parameters.get_objects_linemeas_methods
    Parameters.get_spectrum_model_section

  
Exceptions
----------
Constructor and getters can throw an ``APIException``. Throwing another exception should be considered as a bug (report to amazed-support@lam.fr)

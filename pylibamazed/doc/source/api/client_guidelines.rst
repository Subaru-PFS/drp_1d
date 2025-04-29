Client Developer Guidelines
===========================

Running Amazed consists in:
 * Establishing a **context**
 * Running three main steps

  1. Build and load a **spectrum reader** class
  2. Run the **process flow**
  3. **Write** the results

The context
-----------

You need a **calibration** directory.
This calibration must be compatible with the parameters given to the reader and process flow.

This parametrization is described in detail in :doc:`../../parameters` page. It is a json file later uploaded in a python dictionnary.

API allows to **check** the parameters.

The 3 main steps
----------------
 
Load the spectrum
~~~~~~~~~~~~~~~~~

You need at least three vectors of same size:
    - wavelength
    - flux
    - error (or variance)

Depending on your parametrization, you can also load (same size):
    - lsf
    - photometry
    - other columns used to filter spectrum

You need to define a :doc:`spectrum reader class <methods/input>` that will allow to load this data, and return the spectrum objects to be processed.

If the wavelengths are in air and not in vacuum, you must specify it in the parameters, choosing a `conversion method <../json-schema/general.html#property-general-airvacuummethod>`_.

Multiple spectrum can be processed at the same time. This is called **multiobs** and will be detailed in a separate section.

Run the process flow
~~~~~~~~~~~~~~~~~~~~~~~~

You build a process flow from:
 - the parameters
 - a calibration directory path

.. todo add a link to ResultStoreOutput ref

The run can then be launched. It takes the spectrum object as input, and returns a ``ResultStoreOutput`` object. This object contains the results of the processing, including datasets, attributes and errors.

More details in :doc:`methods/run`.

Writing the results
~~~~~~~~~~~~~~~~~~~

.. todo improve this section, and add : remove the links

A **dataset** is a set of attributes of same dimension. They are all defined in result specifications. (link). The output API is defined in more details in Output API (to be written)



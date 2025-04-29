Output
======

hdf5 output
-----------

.. currentmodule::  pylibamazed.H5Writer

To write a hdf5 file output, use ``H5Writer``

.. autosummary:: 
  H5Writer

Its public methods are:

.. autosummary:: 
  H5Writer.write_hdf5

Public attribute is :

.. autosummary:: 
  H5Writer.excluded_datasets


Minimal example, with ``output`` the output of a process flow run.

.. code:: python

    import h5py
    from pylibamazed.H5Writer import H5Writer
    

    with h5py.File("/some/path/output.hdf5", "w") as output_file:
        writer = H5Writer(output)
        writer.excluded_datasets = config.excluded_datasets
        writer.write_hdf5(output_file, "some group name")


PdfHandler
----------

.. currentmodule:: pylibamazed.PdfHandler

If you want to write a PDF file output, use ``PdfHandler``.
First build the ``PdfHandler`` object with a ``BuilderPdfHandler`` object.

.. autosummary:: 
  PdfHandler

Public methods:

.. autosummary:: 
  PdfHandler.convertToRegular

Public attributes:

.. autosummary:: 
  PdfHandler.redshifts
  PdfHandler.valProbaLog

Module methods:

.. autosummary:: 
  get_final_regular_z_grid


BuilderPdfHandler
------------------

.. autosummary::
  BuilderPdfHandler

Public methods: 

.. autosummary::
  BuilderPdfHandler.add_params
  BuilderPdfHandler.build
  
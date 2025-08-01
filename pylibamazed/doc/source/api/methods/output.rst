Output
======

.. contents::
   :local:
   :depth: 2
   :backlinks: entry

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
  
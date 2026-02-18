Input
=====

.. contents::
   :local:
   :depth: 2
   :backlinks: entry

Implementing the Reader
-----------------------

In order to manage input spectrum, the developer has to implement a derived class from ``AbstractSpectrumReader``.

.. currentmodule::  pylibamazed.AbstractSpectrumReader

.. autosummary::
    AbstractSpectrumReader

Following methods need to be implemented:

.. autosummary::
    AbstractSpectrumReader.load_wave
    AbstractSpectrumReader.load_flux
    AbstractSpectrumReader.load_error
    AbstractSpectrumReader.load_lsf
    AbstractSpectrumReader.load_photometry
    AbstractSpectrumReader.load_others
    AbstractSpectrumReader.set_air_or_vaccum

If necessary, it is possible to overload the constructor :

.. autosummary::
    AbstractSpectrumReader.__init__

The other public methods are:

.. autosummary::
    AbstractSpectrumReader.load_and_get_spectrum
    AbstractSpectrumReader.get_spectrum
    AbstractSpectrumReader.load_all


And its public attributes are:

.. autosummary::
    AbstractSpectrumReader.source_id

Once the reader is implemented, it can be registered using the ``register_reader`` method.

.. autosummary:: 
    register_reader

Loading spectrum
----------------
The spectrum can be loaded using the load methods above. The ``load_all`` method can only be used if there is only one resource.

Once loaded, ``Spectrum`` object can be retrieved using ``get_spectrum``.

**NB**: ``AbstractSpectrumReader`` is a context manager, so you can use the ``with``` python statement.


Example
-------

A simple example is given below:

.. code:: python
	  
  from pylibamazed.AbstractSpectrumReader import AbstractSpectrumReader

  class AsciiReader(AbstractSpectrumReader):

      def load_wave(self, df, obs_id=""):
          """
          Load the spectral axis in self.wave , units are in Angstrom by default
          :param df: df of the resource where the wave can be found
          """
          self.waves.append(df.Lambda)

      def load_flux(self, df, obs_id=""):
          """
          Load the spectral axis in self.flux , units are in erg.cm-2 by default
          :param df: df of the resource where the wave can be found
          """
          self.fluxes.append(df.Flux)

      def load_error(self, df, obs_id=""):
          self.errors.append(df.Error)

      def load_lsf(self, df, obs_id=""):
          self.lsf_type = "no_lsf"



.. code:: python

    from pylibamazed.Parameters import Parameters
    from pylibamazed.CalibrationLibrary import CalibrationLibrary
    from pylibamazed.ProcessFlow import ProcessFlow
    from pylibamazed.Spectrum import Spectrum

    import pandas as pd

    def get_spectrum(path: str, params_dict: dict, calibration_dir :str)->Spectrum:
        spectrum_df = pd.read_csv(path)
        parameters = Parameters(params_dict)
        calib = CalibrationLibrary(parameters,calibration_dir)
        reader = AsciiReader(parameters, calib, "x")
        reader.load_all(spectrum_df)
        return reader.get_spectrum()


The Spectrum object
-------------------

.. currentmodule::  pylibamazed.Spectrum

.. autosummary::
    Spectrum

Its public methods are:

.. autosummary::
    Spectrum.get_flux
    Spectrum.get_wave
    Spectrum.get_error
    Spectrum.get_others
    Spectrum.get_lsf
    Spectrum.get_photometric_data
    Spectrum.get_spectrum_infos 
    Spectrum.init


Its public attributes are:

.. autosummary::
    Spectrum.source_id
    Spectrum.observation_ids

Loading spectrum
----------------

.. currentmodule::  pylibamazed.AbstractExternalStorage

In order to open spectrum files and return their data for the reader, an ``AbstractExternalStorage`` class is provided.

.. autosummary::
    AbstractExternalStorage

This abstract class needs to be implemented with the following methods:

.. autosummary::
    AbstractExternalStorage.read
    AbstractExternalStorage.close

You must then use the ``register_storage`` method to register your implementation of ``AbstractExternalStorage``.

.. autosummary::
    register_storage


Utils
-----

Container
~~~~~~~~~

.. currentmodule::  pylibamazed.Container


An object type ``Container`` is defined.
Its aim is to be a multi-dimensional list with a key-value system.
It is mainly used for multi obs spectrum.

.. autosummary:: 
    Container

``Container`` public methods are, for any data type ``T``:
 
.. autosummary:: 
    Container.append
    Container.get
    Container.keys
    Container.size

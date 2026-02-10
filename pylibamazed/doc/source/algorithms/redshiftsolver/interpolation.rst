.. _interpolation_section:

Interpolation
=============

For rebinning a spectrum or a template, various methods have been implemented :

- :ref:`linear_interpolation`
- :ref:`full_linear_interpolation`
- :ref:`nearest_grid_point_interpolation`
- :ref:`spline_interpolation`
- :ref:`rebin_fine_grid`

Interpolation methods
---------------------

.. _linear_interpolation:

Linear interpolation
~~~~~~~~~~~~~~~~~~~~

A simple linear interpolation is made on the selected range. For all elements not in range, mask is set to 0 (and to 1 inside range).

Flux is interpolated as follows:

.. math::
    y = y_k + (y_{k+1} - y_k) \frac{x - x_k}{x_{k+1} - x_k}

Noting :math:`t = \frac{x - x_k}{x_{k+1} - x_k}`, we deduce the variance of the interpolated flux.

.. math::
   y = y_k \cdot (1 - t) + t \cdot y_{k+1}
   Var(y) = (1-t)^2 Var(y_k) + t^2 Var(y_{k+1})

A compensation factor is then applied (see :ref:`compensation-factor`).


.. _full_linear_interpolation:

Masked linear interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The masked linear interpolation method is a variation of the linear interpolation.

It is to be used on an input spectrum containing masked pixels defined with an input mask. With this method the mask is propagated through the interpolation: a wavelength having a neighbor pixel masked is not interpolated and the interpolated spectrum will have the corresponding sample flagged as masked.

This rebin method is used for rebinning input spectrum if parameter `fftProcessing` is set to `true`.

.. comment : this masked linear interpolation corresponds to full linear interpolation in the code


.. _nearest_grid_point_interpolation:

Nearest grid point
~~~~~~~~~~~~~~~~~~

A nearest grid point interpolation is applied on the selected range, for both the flux and the error.

A compensation factor is then applied on the error computation (see :ref:`compensation-factor`).

.. _spline_interpolation:

Spline
~~~~~~

A spline interpolation is applied on the selected range, for each element.

Note that noise rebinning is not implemented for spline interpolation.

.. _rebin_fine_grid:

Rebin fine grid
~~~~~~~~~~~~~~~

The template is finely sampled (with :math:`\text{sub}\_{\text{step}} = 0.1 \, \text{initial}\_{\text{step}}`) and interpolated using cubic spline once for all, using the `gsl` library.
When rebinning, for each element of the subgrid, the closest interpolated value is used (nearest grid point), and no more calculation is needed.

Note that as for spline interpolation, noise rebinning is not implemented for rebin fine grid.



.. _compensation-factor:

The compensation factor
-----------------------

The steps sizes can vary, especially in the case of log-sampling, and this impacts the variance values. In order to compensate this, the variance is ponderated by a compensation factor :math:`c = \frac{\Delta x_{target}}{\Delta x_{source}}` defined for each target step.

This compensation factor is used for each of these methods for which the variance is interpolated. In practice this compensation factor is only used in the case of FFT processing, when the input spectrum needs to be resampled because its sampling is not logarithmic.


.. _continuum:

Continuum estimation
~~~~~~~~~~~~~~~~~~~~

Depending on the chosen continuum processing method (parameter `continuumComponent`), the continuum estimation can be done in different ways :

- noContinuum
- fromSpectrum
- tplFit
- powerLaw
- auto : tplFitAuto / powerLawAuto

.. _nocontinuum:

No continuum
------------

The continuum is considered non significant, and we approximate :math:`\mathbf{y_c} \approx \mathbf{0}` leading to :math:`\mathbf{d} \approx \mathbf{y_l} + \boldsymbol{\eta}`.


Note that in this case absorption lines are not fittable as they are proportional to the estimated continuum, except if the linefitting method enables the fit of a polynomial under each line.


From spectrum
--------------

The continuum is estimated directly from the input spectrum through an iterative median smoothing algorithm described below. As the continuum is not fitted, it isn't taken into account in the chi2 computation: the total chi2 corresponds to the chi2 of the fit between the linemodel and the spectrum with continuum substracted.

First, in order to avoid border effects, the spectrum is extended by reflection. Two reflection methods are available : odd reflection and even reflection (chosen through the parameter :math:`medianEvenReflection`). Even reflection extends the signal by pure mirror symmetry. Odd reflection also reflects the signal, but about a fitted boundary value, not about the data itself. This allows a better continuity at the boundary, especially when the signal has a strong slope at the borders, but is more complex as it needs to estimate a smoothed boundary value.

Then, the algorithm uses several median smoothing cycles (5 of kernel size  :math:`2 \, w + 1` and 5 of kernel size :math:`w` mith :math:`w` the median kernel width) and a mean smoothing (kernel size :math:`0.25 \, w`) to estimate the continuum. This allows to first remove spikes, and then remove Gaussian noise and residual small scale fluctuations. 

Template fitting
----------------

For template fitting continuum estimation, we use a pre registered catalog of continuum templates, and find the best template and amplitude using the same method than described in :ref:`template fitting solver <templatefitting>`.

However, to jointly fit the continuum amplitude and the lines amplitudes one should in principle do a matrix inversion with an SVD for instance. This is rather costly and does not benefit from the fact that the line model contains many zeros outside the line support. The alternative is to fit interdependently the continuum and the lines by trying to decouple the two.

Two options are available to differentiate the continuum from the linemodel : 

- fit the continuum using only the pixels which are **outside** the lines supports
- **orthogonalization** : use a modified version of the continuum template, orthogonal to the line model (null cross-product). To do so we fit the lines from the line catalog on the continuum template. We then substract these fitted lines from the initial continuum template, and call the result the orthogonalized template. The obtained orthogonalized templates will be the continuum templates used for continuum fitting.


Power law
---------

Power law is a fitting method made for the quasar model.
The continuum is estimated by fitting two power laws on each side of a cut wavelength (5400 Å). See :ref:`power law section <powerlaw>` for more details.

.. _auto:

Auto
------

This auto option applies to template fitting and power law methods. After fitting the continuum at all redshifts, the quality of the fit is assessed using three criteria, described below, in this order:

- **reduced chi2** : if the lowest reducecd chi2 of the continuum fits (looping on redshift & continuums, if several) is above the threshold set by parameter `continuumFit.badChi2Threshold`, the continuum fit is considered bad, and the continuum estimation method is forced to `fromSpectrum`.
- **negative amplitude** : if all continuum fits at all redshifts have amplitudes lower than :math:`n_{neg}` times the computed amplitude uncertainty (:math:`n_{neg}` being set to parameter `continuumFit.negativeThreshold`), the continuum amplitude is considered negative, and the continuum estimation method is forced to `fromSpectrum`. Applies only if reduced chi2 is below the threshold.
- **not significant amplitude** : if all continuum fits have amplitudes lower than :math:`n_{null}` times the computed amplitude uncertainty (:math:`n_{null}` being this time set to parameter `continuumFit.nullThreshold`), the continuum amplitude is considered not significant, and the continuum estimation method is forced to `noContinuum`. It applies only if reduced chi2 is below the threshold, and if amplitude is above negative amplitude threshold, which leads to the condition :math:`snr_{neg} < snr < snr_{null}`.

Note that in all cases, these comparisons are made for template fitting and power law : if they fail and auto is not enabled, an error is thrown.

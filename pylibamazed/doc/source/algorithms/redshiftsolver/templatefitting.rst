.. _templatefitting:

Template fitting
================

Mathematical description
------------------------

Template fitting is a common operation in classical spectroscopic redshift methods (see :cite:`tonry_survey_1979`). Here, the fitting of the spectra consists in estimating :

- :math:`a(z, \boldsymbol{\theta})` the optimal scale parameter
- :math:`\chi^2` the least-squares metric

With :

- :math:`z` the redshift
- :math:`\boldsymbol{\theta} = \{k_{template},EBV_{ISM}, k_{IGM}\}` the template parameters, containing the template index (i.e.the chosen template among a set of templates), the :term:`ISM` dust attenuation and the :term:`IGM` extinction correction terms.

We use the  :ref:`linear least square fitting method <leastsquare>`.

In practice
-----------

Simple algorithm description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- List the templates corresponding to the selected :term:`spectrum model` (from `templateDir` parameter).  
- **For each template** of the template catalog, loop on all redshifts.The redshifts are defined by the `lambdaRange` and `redshiftStep` parameters.
   
   - **For each redshift**: 
   
      - Check that the template wavelength range is able to cover a large enough fraction of the blueshifted input spectrum spectral axis. The minimal required fraction is set via the parameter `overlapThreshold`. If not, throw and error. This common wavelength range will be used for the next computations.
      - **Rebin the template** on this intersected range accordingly to the method selected in the `interpolation` parameter (see :ref:`interpolation_section`).
      - **Find the best igm / ism combination** (if activated), and store template computed amplitude and chi2 informations. For each igm / ism (if activated through parameters `igmFit` and `ismFit`):
   
         - Apply igm / ism on the rebined template
         - Compute optimal template amplitude and corresponding chi2 according to equations above. Store these data.
      - Store the igm / ism combination with the lowest chi2 as the best result (for this template, for this redshift). Compute and add quality fit informations to this result (see section fit quality indicators), as well as `cstLog` (see :ref:`pdf` section) for later pdf computations.
- Once all templates and redshifts have been processed, compute pdf. For this:

   - Build chi square array, a 2d array of first dimension :math:`n_{templates} \cdot n_{igm} \cdot n_{ism}`, and second dimension :math:`n_{redshifts}`. Fill each cell with the chi2 value computed above.
   - From this chi2 array, compute **pdf** as described in dedicated section, depending on the `pdfCombination` parameter, and retrieve the **best candidates** (see :ref:`pdf` section), which number correspond to the `extremaCount` parameter. Each selected candidate must be distant of at least the parameters' `extremaRedshiftSeparation` value from each other. Store these results.
   - For each candidate solution retrieve the best model (template, ism, igm) as the one with the lowest chi2. Recompute the model (rebin the template, apply ism, igm and scale with the estimated amplitude), and save it.
   

Two pass computation
~~~~~~~~~~~~~~~~~~~~~

In order to save computation time, this fitting process can be made in two pass: a first pass with a coarser grid to get a rough estimate of the best redshifts, and a second pass with a thiner grid around the selected candidates to refine the best results.

- Compute the **first pass** as described above, except that: 

   - The redshifts steps are multiplied by the parameter `firstPass.largeGridStepRatio`
   - The number of candidates selected and stored during first pass is `firstPass.extremaCount`
   - Each candidate must be distant of at least twice the value of `secondPass.m_opt_secondpass_halfwindowsize` parameter

- Compute the **second pass** around the first pass candidates

   - Create the **new redshift grid**. For each first pass candidate, create a subgrid centered on its redshift spanning across +/- parameter `secondPass.halfWindowSize`, with step `redshiftStep`.
   - Loop on first pass candidates, and recompute the model on each corresponding subgrid. For each first pass candidate:

      - Use same best template, igm, ism than registered candidate.
      - As before, compute optimal template amplitude and corresponding chi2 according to equations above.

- Compute **chi square array** as before (on the full grid containing both the coarse grid and the subgrids). Note that this time, in each redshift subgrid, we didn't compute the model and chi2 for each igm and ism but only for the ones obtained from the first pass (for each candidate). We make here an approximation and interpolate the chi2 values between the redshifts indexes of the coarse grid for the igm / ism which were not pre-computed.
- Compute **pdf**, this time without constraint on peak separation, but retain only one solution (the maximum of the pdf) per candidate window. Note that if the maximum of the pdf is on the border of the window, this candidate is rejected. Retrieve the number of candidates defined in `extremaCount` (at most).


FFT
~~~

The redshift grid is computed as a log-sampled grid : 

.. math::
   z_i = (z_0 + 1)^{i\,\delta} - 1

With:

   - :math:`i` the index of the redshift
   - :math:`z_0` the first redshift of the grid (min of parameter `lambdaRange`)
   - :math:`\delta` the step of the grid (parameter `redshiftStep`)

   


.. comment
   Photometry
   ~~~~~~~~~~

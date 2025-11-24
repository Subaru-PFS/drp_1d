.. _linemodelsolve:

Line model solve
=================

The line model solve algorithm works in two steps :

- First the continuum of the observed spectrum is estimated, and substracted from the input spectrum
- A set of line profiles is then fitted on the resulting spectrum.

From the resulting model (continuum + lines), we can then compute pdf from the merit values, and obtain the best candidates.


We note :

.. math::
   \mathbf{d} = \mathbf{y_c}(z, \boldsymbol{\theta}) + \mathbf{y_l}(z, \boldsymbol{\omega}) + \boldsymbol{\eta}

Where :

- :math:`\mathbf{d}` is the observed data
- :math:`\mathbf{y_c}` is the modeled continuum component
- :math:`\mathbf{y_l}` is the modeled lines component
- :math:`\boldsymbol{\eta}` is additive Gaussian noise
- :math:`z` is the redshift
- :math:`\boldsymbol{\theta}` are the other model parameters used for continuum fitting
- :math:`\boldsymbol{\omega}` are the other model parameters used for line fitting

And: 

.. math::
   \mathbf{y_l} = \mathbf{y_{em}} + \mathbf{y_{abs}}

Where:

- :math:`\mathbf{y_{em}}` is the modeled emission lines component
- :math:`\mathbf{y_{abs}}` is the modeled absorption lines component

With:

.. math::
   \mathbf{y_em} = \sum_{i \in lines_{em}} e_i\, \mathbf{p_em}_i

.. math::
   \mathbf{y_abs} = \sum_{i \in lines_{abs}} a_i\, \mathbf{p_{abs}}_i \cdot \boldsymbol{y_c}

- :math:`e_i` / :math:`a_i` is the emission/absorption line amplitude
- :math:`\mathbf{p_em}_i` / :math:`\mathbf{p_{abs}}_i` is the emission/absorption line profile



.. toctree::

   linemodelsolve/continuum
   linemodelsolve/linefitting




Whole linemodelsolve algorithm, in practice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- First, set redshift grid according to parameters `redshiftStep` and `redshiftRange`
- Then load, the necessary templates and catalogs corresponding to our spectrum model
- If the continuum fit method is `tplFit` or `powerLaw` (or their `auto` counterparts), compute the continuum:

  - if lines support should be ignored (from parameter `continuumFit.ignoreLineSupport`), create a mask builder which will mask all lines supports for continuum computation
  - if the continuum fit method is powerlaw: 

    - compute "best" power law coefficients for each redshift, and store the corresponding igm / ism indexes, using  :ref:`power law method<powerlaw>`

  - if the continuum fit method is tplFit:

    - orthogonalize templates
    - for each template, compute "best" amplitude for each redshift, and store the corresponding igm / ism indexes, using  :ref:`template fitting method<templatefitting>`

  - Update the continuum fitting methods as described in :ref:`continuum auto section <auto>` if necessary

- Loop on redshifts

  - set lya profile according to parameter `lya.profile`, and fit the igm (see :ref:`lya`)
  - prepare continuum
  - fit according to the selected method line ratio and line fittinf methods
  - compute merit
  - store these results

- As in template fitting method, compute the chi2 array, corresponding pdf, and retrieve the best candidates. Retrieve the best model, and save it.

.. continue : describe how chi2 array is built for linemodelsolve

.. note
   We do not document the possibility nContinuum > 1 (continuumFit.count) for the moment

Like in the template fitting solver, in order to save computation time, this fitting process can be made in two pass.

.. comment - to be developed
   The principle is the same in line model solve than in template fitting solve, with some added steps:

   - Compute first pass : 

      - The redshift grid has wide steps, each redshift step is defined as the product of parameters `redshiftStep` and `firstPass.largeGridStepRatio`
      - The number of candidates selected and stored during first pass is `firstPass.extremaCount`
      - Each candidate must be distant of at least twice the value of `secondPass.m_opt_secondpass_halfwindowsize` parameter
   

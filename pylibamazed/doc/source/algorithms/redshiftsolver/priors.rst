.. _priors:

Priors
======

We define priors as :

.. math::
    \Pi(z, \boldsymbol{\theta}) = p(z, \boldsymbol{\theta} \mid \text{'gal'}) = p(z \mid \text{'gal'}) \times p(\boldsymbol{\theta} \mid z, \text{'gal'})

We will call:

- **z-priors** : :math:`p(z \mid \text{'gal'})`
- **model priors** : :math:`p(\boldsymbol{\theta} \mid z, \text{'gal'})`

Note that after computation, all priors are then normalized: i.e. divided by their integral (evaluated using the trapezoidal rule) so that the integral of the prior is equal to 1.

z-priors
--------

z-priors only depend on the redshift and not on the model parameters. They can be used to favor or disfavor some redshiftsranges. The different available z-priors are :

- :ref:`uniform_prior`
- :ref:`halpha_prior`
- :ref:`pozzetti_prior`
- :ref:`strong_line_prior`

.. _uniform_prior:

Uniform prior
^^^^^^^^^^^^^

When no other z-prior is specified the uniform prior is used. 

.. math::
    p(z \mid \text{'gal'}) = \frac{1}{n_{z}}

With :math:`n_{z}` the number of redshift steps.

.. _halpha_prior:

Halpha prior
^^^^^^^^^^^^

The Halpha prior favors redshifts where the Halpha emission is the strongest observed line.
If Halpha is the strongest line, set prior to 1. Otherwise, set prior to penalization factor from `hAlphaPrior` parameter.

.. _pozzetti_prior:

Pozzetti prior
^^^^^^^^^^^^^^^

The Pozzetti prior, from :cite:`Pozzetti_2016`, provides an estimate of Nz (source distribution as a function of z) for Ha emitters. It is activated with parameter `nOfZPriorStrength`.

.. _strong_line_prior:

Strong line prior
^^^^^^^^^^^^^^^^^

The strong line prior favors redshifts where at least one strong emission line is detected.
If there is at least one strong line, set prior to 1. Otherwise, set prior to penalization factor from `strongLinesPrior` parameter.


.. comment : sections to add :  Model priors & Combining priors

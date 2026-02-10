Overview
==========

The current redshift estimation algorithms described in this document are all based on a least-square minimization of a set of redshift-dependent templates/models on the observed spectrum.

Likelihood
--------------

The determination of the redshift is based on the use of redshift-dependent
models that are compared to the observed spectrum through least-square
fitting. For each model, a likelihood can be obtained:

.. math::

   -2 \log \mathcal{L}(z)
   = \min_{\theta} \sum_{\lambda}
   \frac{\left( S_{\lambda} - M_{\theta}\!\left( \frac{\lambda}{1+z} \right) \right)^{2}}
        {\sigma_{\lambda}^{2}}

Where :

- :math:`\lambda` is the observed wavelength,
- :math:`S_{\lambda}` is the observed spectrum,
- :math:`\sigma_{\lambda}^{2}` is the known noise variance, and
- :math:`M_{\theta}\!\left(\frac{\lambda}{1+z}\right)` is the redshifted model spectrum, with :math:`\theta` the vector of parameters.

Main steps
------------

The typical steps are:

#. **Chi2 computation** : for each redshift, fit all available models, and retain their chi2
#. **PDF computation** : combine the chi2 of all models on a zPDF
#. **Candidates selection** : determine the maxima of the zPDFs. Then for each candidate retrieve the best fitted model i.e. the one with the lowest chi2

Two-pass computation
---------------------

In order to save computation time, a two-pass computation is proposed :

* a first pass with the 3 steps described above, but with a coarser redshift grid, and less free parameters (less models)
* from the set of first pass candidates do a 2nd pass on a fine grid arround the candidate, with more models.

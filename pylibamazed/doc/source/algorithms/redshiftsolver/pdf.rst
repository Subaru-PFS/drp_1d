.. _pdf:

Probability density function
============================

We want to compute the log-pdf :math:`p(z, \boldsymbol{\theta} | \mathbf{d}, \text{'gal'})` for all parameters, which coresponds to the a posteriori probability since we know the realization :math:`\mathbf{d}` of the data.

Mathematical description
------------------------

Notations
^^^^^^^^^

We have the possibility to use one of several models ("galaxy", "star", "QSO", ...). Noting :math:`\mathcal{M}` the model variable and :math:`m` its value, we write :math:`p(\mathcal{X}=x \mid \mathcal{M}=m) = p(x \mid m) = p(x \mid \text{'gal'})`, for clarity.

.. _bayesian-framework:

Bayesian framework
^^^^^^^^^^^^^^^^^^

From Bayes theorem : 

.. math::
  p(z, \boldsymbol{\theta} | \mathbf{d}) = \frac{p(\mathbf{d} | z, \boldsymbol{\theta})p(z, \boldsymbol{\theta})}{p(\mathbf{d})}

Adding a dependency on the model, we obtain the a posteriori probability:

.. math::

  p(z, \boldsymbol{\theta} | \mathbf{d}, \text{'gal'}) = \frac{p(\mathbf{d} | z, \boldsymbol{\theta}, \text{'gal'})p(z, \boldsymbol{\theta} | \text{'gal'})}{p(\mathbf{d} | \text{'gal'})}

Where:

- :math:`p(\mathbf{d} | z, \boldsymbol{\theta}, \text{'gal'}) = \mathcal{L}(z, \boldsymbol{\theta})` is the likelihood of the parameters
- :math:`p(z, \boldsymbol{\theta}| \text{'gal'}) = \Pi(z, \boldsymbol{\theta})` is the prior (see the priors section :ref:`priors <priors>`).
- :math:`p(\mathbf{d}∣\text{'gal'})` is the evidence of the Galaxy model

Which leads to the log-pdf equation :

.. math::
  :label: first-log-pdf

  \log(p(z, \boldsymbol{\theta} | \mathbf{d}, \text{'gal'})) = \log(\mathcal{L}(z, \boldsymbol{\theta})) + \log(\Pi(z, \boldsymbol{\theta})) - \log(p(\mathbf{d} | \text{'gal'}))


Computing log-likelihood
~~~~~~~~~~~~~~~~~~~~~~~~

From the definition of the likelihood of a gaussian noise, we get : 

.. math::
    \mathcal{L}(z, \boldsymbol{\theta})  = \prod_{i \in pixels} \frac{1}{\sqrt{2\pi}\sigma_i} \exp{(-\frac{(d_i-m_i(z, \boldsymbol{\theta}))^2}{2\sigma_i^2})}

Which leads to log-likelihood:

.. math::
  :label: log-likelihood

    \log{\mathcal{L}(z, \boldsymbol{\theta})} = - \frac{1}{2} \chi^2(z, \boldsymbol{\theta}) + \sum_{i \in pixels}\log{\frac{1}{\sqrt{2\pi}\sigma_i}}
    
We note :

.. math::

  cstLog = \sum_{i \in pixels} \log{\frac{1}{\sqrt{2\pi}\sigma_i}}


Computing the evidence
~~~~~~~~~~~~~~~~~~~~~~

We compute the evidence, beginning from marginalization definition :

.. math::
  p(\mathbf{d} | \text{'gal'}) = \int_z \int_{\boldsymbol{\theta}} p(\mathbf{d}, z, \boldsymbol{\theta} | \text{'gal'})\, dz\, d\boldsymbol{\theta}

Using bayes theorem :

.. math::
  p(\mathbf{d} | \text{'gal'}) = \int_z \int_{\boldsymbol{\theta}} p(\mathbf{d} | z, \boldsymbol{\theta}, \text{'gal'})\, p(z, \boldsymbol{\theta} | \text{'gal'}) \, dz\, d\boldsymbol{\theta}

Using gaussian noise definition, we finally obtain the evidence : 

.. math::
  :label: evidence

  p(\mathbf{d} | \text{'gal'}) = \int_z \int_{\boldsymbol{\theta}} \Pi(z, \boldsymbol{\theta}) e^{cstLog} e^{-\frac{1}{2}\chi^2(z, \boldsymbol{\theta})}\, dz\, d\boldsymbol{\theta}


Full log-pdf equation
~~~~~~~~~~~~~~~~~~~~~

From :eq:`first-log-pdf`, :eq:`log-likelihood` and :eq:`evidence`, we finally obtain the full log-pdf equation :

.. math::
  :label: log-pdf

  \log(p(z, \boldsymbol{\theta} | \mathbf{d}, \text{'gal'})) = -\frac{1}{2}\chi^2(z, \boldsymbol{\theta}) + \log(\Pi(z, \boldsymbol{\theta})) - \int_z \int_{\boldsymbol{\theta}} \Pi(z, \boldsymbol{\theta}) e^{-\frac{1}{2}\chi^2(z, \boldsymbol{\theta})}\, dz\, d\boldsymbol{\theta}


Combining models
^^^^^^^^^^^^^^^^

Combining the information corresponding to all models created with different parameters can be
achieved in two ways depending on the parameter `pdfCombination` : 

- **Best proba solution** for each redshift:

  .. math::
      :label: na-best

      p(z, \hat{\boldsymbol{\theta}} | \mathbf{d}, \text{'gal'})
        = \max_{\boldsymbol{\theta}}
          \frac{p(\mathbf{d} \mid z, \boldsymbol{\theta}, \text{'gal'})\,\Pi(z,\boldsymbol{\theta})}
               {p(\mathbf{d} \mid \text{'gal'})}

The best model parameters at a given redshift being:

.. math::

  \hat{\boldsymbol{\theta}}(z) = \arg\max_{\boldsymbol{\theta}} p(z, \boldsymbol{\theta} \mid \mathbf{d}, \text{'gal'})

When priors are uniform (i.e. no prior), it corresponds to `pdfCombination` set to `bestChi2`.

- **Marginalized solution**:

  .. math::
      :label: na-marg

      p(z \mid \mathbf{d}, \text{'gal'})
        = \int_{\boldsymbol{\theta}}
          \frac{p(\mathbf{d} \mid z, \boldsymbol{\theta}, \text{'gal'})\,\Pi(z,\boldsymbol{\theta})}
               {p(\mathbf{d} \mid \text{'gal'})} d\boldsymbol{\theta}

This corresponds to `pdfCombination` set to `marg`.

In the case of a discrete parameter in :math:`\boldsymbol{\theta}`, for instance the choice of a template among a set of templates, the integral on this parameter becomes a discrete sum.


Marginalization gets the best redshift probability and accepts priors, whereas the best chi2 maximizes the probability of the couple :math:`(z,\boldsymbol{\theta})` and does not accept priors. The `bestChi2` pdf combination is to be preferred when the parameters :math:`{\boldsymbol{\theta}}` may be completely uncorrelated between models, for instance if the velocity dispersion can be very different from one template to the other (for instance for stars). The marginalization method should be used in all other cases.



Determining and ranking the redshift solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The redshift candidates are then determined on the combined pdf above, using the :math:`n` first local maxima (`extremaCount` parameter).

For each candidate, a redshift uncertainty is obtained by fitting a Gaussian on the local maxima. Then all the candidates are ranked by the pdf integrated around the solution :math:`\pm 3 \, \sigma` of the fitted Gaussian.

  .. math::
      :label: zbest-maxint

      z_{\text{best}}
        = \arg\max_{z_p}
          \sum_{z_p - \delta \le z \le z_p + \delta}
          p(z \mid \mathbf{d})

With :

- :math:`z_p` the redshifts of the local maxima
- :math:`\delta = 3\,\sigma` the window on which to make the integration

Fitting the gaussian
~~~~~~~~~~~~~~~~~~~~~~

We analytically fit a parabole on each of the local maxima. 

The log-pdf is approximated locally by a parabola constrained to have:

- its maximum at :math:`z = z_0` the redshift of the local maximum
- its value fixed at :math:`m_0` the log-pdf value at :math:`z = z_0`

The only unknown parameter is the curvature coefficient :math:`c_0` : :math:`m_0 - m(z) \;\approx\; c_0 (z - z_0)^2`

Using a least square fitting method, we obtain : 

.. math::
  c_0 = \frac{\sum_i (z_i - z_0)^2 \, (m_0 - m_i)}
       {\sum_i (z_i - z_0)^4}

With :

- :math:`z_i` the redshifts of surrounding pixels
- :math:`m_i` the log-pdf values at these redshifts

From this curvature coefficient we obtain the gaussian width :math:`\sigma =  \sqrt{\frac{1}{2 \, c_0}}`.

We arbitrarily choose to use 5 pixels on each side of the maximum for the fit, then 2 if the fit is not successful, then if the fit is still not successful, :math:`\sigma` is forced to :math:`1e^{-3}`.

The fit is not successful if the curvature coefficient :math:`c_0 \le 0` (i.e. the parabola is not concave).


Integration around the solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that the integration range is defined as :math:`[z_{\text{candidate}} - 3 \, \sigma, z_{\text{candidate}} + 3 \, \sigma]`,we compute the integrated pdf for each candidate. For that we use a simple trapezoidal integration.




In practice
-----------

Compute log pdf
^^^^^^^^^^^^^^^

- Using a precomputed chi2 array (computed during the redshift solver operations), **compute log-pdf for each model** i.e. for each parameters combination using formula :eq:`log-pdf`:

  - For each model, use chi2 values and z priors, compute log-evidence  (see :eq:`evidence`), and compute log pdf.

  - Set model priors. If specified, use provided prior (see section :ref:`priors <priors>`). Otherwise, set constant prior. Include this prior information in the evidence : considered evidence will now be obtained through the sum of log-evidence and log-prior.

  - Sum these new evidences from all models : it will be used for classification
   
- **Combine the pdfs** depending on `pdfCombination` parameter. It can be set to `marg` or `bestChi2`. 

    - If set to **marg**, marginalize on :math:`\boldsymbol{\theta}`, using the equation :eq:`na-marg` above. Note that we are using the log-values, so we use a numerical stability trick to compute the marginalization based on the formula :math:`\log(e^x+e^y)=\max(x,y)+\log(e^{x - \max} + e^{y-\max})`

    - If set to **bestChi2**, for each redshift, find the model which minimized the chi2. Re-compute the log-pdf as above (:eq:`log-pdf`) using these merit values in particular. 

Select candidates
^^^^^^^^^^^^^^^^^

The candidates are selected as the maxima of the pdf. The range on which we look for maxima is :
  
  * The whole range in first pass
  * Several smaller windows in second pass. Second pass windows are centered on first pass candidates redshifts and of half width defined in parameter `secondPass.halfWindowSize`


Looping over the redshift range, we find all the **pdf peaks** (see the different parameters customziations). The pdf peaks are selected depending on various parameters :

   - **Allow extrema at border** : a maxima at the border is considered as valid only in first pass. In second pass, we ignore maxima at the border.
   - **Peak separation** : two maxima should be at least distant of this value, set in parameter `extremaRedshiftSeparation` for first pass.
   - **Max candidates** : maximal number of candidates to keep, corresponds to parameter `extremaCount`. It can be set to different values for first and second passes.
   - **Max peak per window** : maximal number of candidates to retain per window. For first pass, the value is set to `firstPass.extremaCount` (the first pass window corresponds to the full z range) whereas for second pass the value is set to one (only one second pass candidate can be retrieved per first pass candidate window). Note that `firstPass.extremaCount` should be higher than `extremaCount` : the idea is to retain more candidates in first pass, then eventually truncate the number of final candidates based on their probability (integrated pdf under peak).
   - **Merit cut** : a peak should have a merit value higher than this value, set in parameter `extremaCutProbaThreshold`. This parameter is defined only for line model solver, and set to 0 for all other solvers.

Once found, the peaks are ranked by integrated pdf under the peak (see below). They will be called **candidates**.

Compute integrated pdf
^^^^^^^^^^^^^^^^^^^^^^
For each candidate, fit a Gaussian and compute pdf integrated around the solution :math:`\pm 3 \, \sigma` of the fitted Gaussian.

- **Fit a gaussian** on the local maxima (see description above). Retrieve its standard deviation :math:`\sigma`. Set the integration range to :math:`z \pm 3 \, \sigma`.

- Look for **overlapping ranges** :

  - If overlap is less than 30%, split equally the overlapping part between both candidates. 
  - Otherwise, remove the candidate with the lowest pdf value.

- Compute the **integrated pdf** on these ranges.

- Keep only the "best" candidates i.e. **truncate** the number of candidates to the ones with the highest integrated pdf (see max candidates describes above).

- In case the integration range of some of the selected candidates goes beyond the defined window ranges, issue a warning.

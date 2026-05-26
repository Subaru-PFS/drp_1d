.. _leastsquare:

Linear least Square Fitting
===========================

Across the different algorithms, we use the linear least square fitting method for fitting amplitudes. We will describe more precisely the equations and notations here. 

We wish to estimate :

- :math:`a(z, \boldsymbol{\theta})` the optimal scale parameter
- :math:`\chi^2` the least-squares metric

With :

- :math:`z` the redshift
- :math:`\boldsymbol{\theta}` the other parameters of the model

The optimal scale parameter is estimated by considering the following model:

.. math::

   \mathbf{d} = a(z, \boldsymbol{\theta})\, \mathbf{y}(z, \boldsymbol{\theta}) + \boldsymbol{\eta}

where  :

- :math:`\mathbf{y}(z, \boldsymbol{\theta})` the template modeled flux (without amplitude)
- :math:`\mathbf{d}` is the observed data
- :math:`\boldsymbol{\eta}` is additive Gaussian noise


We note :math:`N` the noise covariance matrix. Supposing that the pixels are uncorrelated, we only consider the variance. :math:`N` is a diagonal matrix, and we note :math:`N = \mathrm{diag}(\sigma_1^2, \sigma_2^2, \ldots, \sigma_n^2)` with :math:`\sigma_i` the flux standard deviation of pixel :math:`i`.


Using the least-square fitting method, the estimated scale parameter is obtained by:

.. math::
   :label: eq:scale_param

   \hat{a}(z, \boldsymbol{\theta}) =
   \frac{\mathbf{d}^T N^{-1} \mathbf{y}}{\mathbf{y}^T N^{-1} \mathbf{y}}

the uncertainty on the estimated scale, is obtained by:

.. math::
   :label: eq:scale_uncertainty

   \text{VAR}(\hat{a}(z, \boldsymbol{\theta})) = \frac{1}{\mathbf{y}^T N^{-1} \mathbf{y}}

and the :math:`\chi^2` metric by:

.. math::
   :label: eq:chi2

   \chi^2 = (\mathbf{d} - \hat{a}\, \mathbf{y})^T N^{-1} (\mathbf{d} - \hat{a}\, \mathbf{y}) =
   \mathbf{d}^T N^{-1} \mathbf{d} +
   \hat{a}^2\, \mathbf{y}^T N^{-1}\mathbf{y}
   - 2 \hat{a}\, \mathbf{d}^T N^{-1} \mathbf{y}

with

.. math::

   \begin{aligned}
   \mathbf{d}^T N^{-1} \mathbf{y}
   &= \sum_{i \in pixels} \sigma_i^{-2}\, d_i\, y_i(z, \theta)

   \mathbf{y}^T N^{-1} \mathbf{y}
   &= \sum_{i \in pixels} \sigma_i^{-2}\, y_i^2(z, \theta)

   \mathbf{d}^T N^{-1}  \mathbf{d}
   &= \sum_{i \in pixels} \sigma_i^{-2}\, d_i^2
   \end{aligned}

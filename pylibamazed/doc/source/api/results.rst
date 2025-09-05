Results
=======

Fit quality indicators
----------------------

The fit quality indicators implemented in the library are listed below. They are computed for both continuum and full model. The names listed below are for full model, for continuum just add the keyword ``Continuum`` before. ex : ``LeastSquare`` / ``ContinuumLeastSquare``.

We use here : 

* :math:`n` number of pixels
* :math:`s_i` input spectrum measured flux at pixel :math:`i`
* :math:`m_i` model flux at pixel :math:`i`
* :math:`\sigma_i` input spectrum flux standard deviation of pixel :math:`i`
* :math:`r_i = \frac{s_i - m_i}{\sigma_i}` residual at pixel :math:`i`

**LeastSquare**
: Merit, sum of squared residuals. :math:`\chi^2 = \sum_{i=0}^{n} r_i^2`

**ReducedLeastSquare**
: Reduced merit, sum of squared residuals divided by number of pixels. :math:`\chi^2_{red} = \frac{\chi^2}{n}`

**PValue**
: p-value of the hypothesis "residuals follow a normal distribution" :math:`\mathcal{N} (0,1)`

**ResidualsMean**
: Mean :math:`\bar{\mu}` of the residuals

**ResidualsStd**
: Standard deviation of the residuals :math:`\bar{\sigma}^2 = \frac{1}{n-1}\sum(r_i - \bar{\mu})^2`

**ResidualsSkewness**
: Skewness of the residuals :math:`\gamma_1 = \frac{1}{n}\sum \left(\frac{r_i - \bar{\mu}}{\bar{\sigma}}\right)^3`

**ResidualsKurtosis**
: Kurtosis excess of the residuals :math:`\gamma_2 = \frac{1}{n}\sum \left(\frac{r_i - \bar{\mu}}{\bar{\sigma}}\right)^4 - 3`

**Ks**
: Kolmogorov-Smirnov statistic of the residuals :math:`D = \sup_i | F_n(r_i) - F(r_i) |` where :math:`F_n` is the empirical cumulative distribution function of the residuals and :math:`F` is the cumulative distribution function of the normal distribution :math:`\mathcal{N}(0,1)`

**KsStd**
: Same as above but the standard deviation of the residuals is used instead of 1. :math:`F` is the cumulative distribution function of the normal distribution :math:`\mathcal{N}(0,\bar{\sigma}^2)`

**KsStdMean**
: Same as above but the mean and standard deviation of the residuals are used instead of 0 and 1. :math:`F` is the cumulative distribution function of the normal distribution :math:`\mathcal{N}(\bar{\mu},\bar{\sigma}^2)`

**Anderson**
: Anderson-Darling statistic of the residuals :math:`A^2 = n \sum \frac{(F_n(r_i) - F(r_i))^2F'(r_i)}{F(r_i)(1-F(r_i))}` where :math:`r_{i}` are the ordered residuals and :math:`F(x)` is the cumulative distribution function of the normal distribution :math:`\mathcal{N}(\bar{\mu},\bar{\sigma}^2)`
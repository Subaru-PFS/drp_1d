Chi2 and Probability Density Function
=====================================

What is ChisquareArray.cstLog ?
-------------------------------

Likelihood is defined as 

.. math::
    L(z)  = \prod_i \frac{1}{\sqrt{2\pi}\sigma_i} \exp{(-\frac{(d_i-M_i(z))^2}{2\sigma_i^2})}

Which leads to log likelihood:

.. math::
    \log{L(z)} = -\frac{N_{pix}}{2}\log{2\pi} - \sum_i\log{\sigma_i} - \frac{\chi^2}{2} = cstLog - \frac{\chi^2}{2}

With :math:`cstLog` - constant with respect to z - defined as :

.. math::
    cstLog = -\frac{N_{pix}}{2}\log{2\pi} - \sum_i\log{\sigma_i}




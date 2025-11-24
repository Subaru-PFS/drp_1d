.. _lineprofiles:

Line profiles
=============

Depending on how they are described in the calibration, lines can be modeled using different profiles:


- `SYM` :  simple gaussian profile
- `ASYM` : asymmetric gaussian profile, with :math:`\sigma = 1`, :math:`\alpha = 4.5` and :math:`\delta = 0`
- `ASYMFIT` : asymmetric gaussian profile with a fit on the tails, with :math:`\sigma = 2`, :math:`\alpha = 2` and :math:`\delta = 0`, with a "mean" centering


Simple Gaussian profile
------------------------

Follows the equation : 

.. math::

   f(x) = e^{-\frac{(x - x_0)^2}{2 \sigma^2}}


Asymmetric gaussian profile
---------------------------

It corresponds to the pdf of the skewed normal distribution, which expression is :

.. math::

   f(x) = e^{- \frac{1}{2} x_{\rm surc}^2} \left[ 1 + \operatorname{erf} \left( \frac{\alpha}{\sqrt{2}} x_{\rm surc} \right) \right]


The variable :math:`x_{\rm surc}` represents the scaled and shifted coordinate in the
skew-normal line profile. It is computed from the original coordinate :math:`x`
and the line center :math:`x_0` as follows:

.. math::

   x_{\rm c} &= x - x_0 \\[2mm]
   x_{\rm surc} &= \frac{x_{\rm c} + \Delta + \sigma m_0}{\sigma}

where:

- :math:`\sigma` is the standard deviation of the line (width parameter).  
- :math:`m_0` is a centering offset that depends on the chosen centering method (see below)
  ("none", "mean", or "mode").  
- :math:`\Delta` is an additional shift parameter

The value of :math:`x_{\rm surc}` ensures that the line profile is properly scaled
and optionally centered depending on the level of asymmetry given by the parameter :math:`\alpha`.

Centering methods
~~~~~~~~~~~~~~~~~~~~


The centering method determines how the line profile is positioned relatively to the central wavelength :math:`x_c`. Two options are "mean" and "none":

None
^^^^^

The line profile is **not shifted**; the center of the profile remains at the nominal line center :math:`x_0`. The resulting profile may appear slightly off-centered if it is strongly
skewed.

Mean
^^^^^
In order to center the profile, we would like to position the maximum of the profile at the central wavelength :math:`x_c`, which corresponds to the skewed normal distribution mode. Since the mode is more complicated and costly to calculate, we use the mean instead, which is close to the mode.

The line profile is **shifted so that its mean matches the line center** :math:`x_0`. A centering offset :math:`m_0` is applied, computed as  

.. math::

    m_0 = \mu_z = \delta \sqrt{\frac{2}{\pi}}, \quad
    \delta = \frac{\alpha}{\sqrt{1 + \alpha^2}}
  
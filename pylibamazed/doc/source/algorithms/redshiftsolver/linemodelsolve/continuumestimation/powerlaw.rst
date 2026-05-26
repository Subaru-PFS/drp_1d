.. _powerlaw:

Power law
=========

.. role:: raw-latex(raw)
     :format: latex html


A fitting method to use for QSOs.
The idea is to fit the flux with two power laws (a first power law for wavelength < 5400 Å, a second one for wavelength > 5400 Å).
We make the fit using the least square fitting method on a log/log scale.


Mathematical description
------------------------


Applying least square fitting method on log/log scale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
     \begin{cases}
          y = a_1 x^{b_1} & \text{for } x < x_c \\
          y = a_2 x^{b_2} & \text{for } x > x_c \\
          a_1 x_c^{b_1} = a_2 x_c^{b_2}
     \end{cases}

With :math:`Y = \ln{y}`,  :math:`X = \ln{x}`, :math:`A = \ln{a}` :

.. math::
     \Leftrightarrow
     \begin{cases}
     Y = A_1 +  b_1 X & \text{for } X < X_c \\
     Y = A_2 + b_2 X & \text{for } X > X_c \\
     A_2 = A_1 + b_1 X_c - b_2 X_c
     \end{cases}

.. math::

     \Leftrightarrow
     \begin{cases}
          Y = A_1 +  b_1  X & \text{for } X \leqslant X_c \\
          Y = A_1 + b_1 X_c + b_2 (- X_c + X) & \text{for } X > X_c \\
     \end{cases}
     \quad (1)

     

To solve this, we use least square fitting method.

With :math:`d` the measured data,  :math:`M` the power law, :math:`\eta` the noise, and :math:`P` the probability:

.. math::
     d = M(\theta) + \eta \Rightarrow P(d - M(\theta)) = P(\eta)

The noise :math:`\eta` is supposedly following a normal law :math:`P(\eta) = \alpha \exp(-\frac{1}{2} \sum (\frac{\eta_i}{\sigma})^2)`.

We want to maximize likelihood :math:`\mathcal{L}(\theta) = P(d \mid \theta) = P(\eta)`, which leads to maximizing the log-likelihood, and therefore minimizing :math:`-2 \ln{\mathcal{L}(\theta)}`.

.. math::
     -2 \ln{\mathcal{L}(\theta)} = -2 \ln{P(\eta)} = -2 \ln{\alpha} + \sum \left( \frac{\eta_i}{\sigma_i} \right)^2

Therefore, maximizing log-likelihood leads to minimizing :math:`\sum (\frac{\eta_i}{\sigma_i})^2`.

.. math::
     \sum \left( \frac{\eta_i}{\sigma_i} \right)^2 = \sum \left( \frac{d_i - M_i}{\sigma_i} \right)^2


Using matrix notation, minimze 

.. math::
     (d-M \cdot \theta)^T N^{-1} (d - M \cdot \theta)
     \quad (2)

With :math:`N` the covariance matrix,  a diagonal matrix with :math:`n_{i,i}=\sigma_i^2`.

Supposing :math:`M(\theta) =  M \cdot \theta` :

.. math::
     \frac{\partial}{\partial{\theta}} \mathcal{L} = 0 \Rightarrow \frac{\partial}{\partial{\theta}} ((d-M \theta)^T N^{-1} (d - M \theta)) = 0 \Rightarrow ... \Rightarrow M^T N^{-1} d = M^T N^{-1} M \theta

Finally:

.. math::
     \boxed{
          \theta =  (M^T N^{-1} M)^{-1} \cdot M^T N^{-1} d
     }


With :

.. math::
     \theta = \left( \begin{array}{c} A_1 \\ b_1 \\ A_2 \\ b_2 \end{array} \right)
     \text{, and }
     M = \left( \begin{array}{cccc} 1 & X_i & 0 & 0 \\ \vdots & \vdots & \vdots & \vdots \\ 0 & 0 & 1 & X_i \\ \vdots & \vdots & \vdots & \vdots \end{array} \right)

Which can be reduced with the continuity constraint of the 2 power laws in **(1)** :

.. math::
     \theta = \left( \begin{array}{c} A_1 \\ b_1 \\ b_2 \end{array} \right)
     \text{, and }
     M = \left( \begin{array}{ccc} 1 & X_i & 0 \\ \vdots & \vdots & \vdots \\ 1 & X_c & -X_c + X_i \\ \vdots & \vdots & \vdots \end{array} \right)


Matrix calculations then allow to find an analytic solution to this equation.

Case were b1 / b2 is fixed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Constraints on min / max value of b1 and b2 can be set in parameters.
If b1 or b2 reaches these limits, we redo the fit with b1 or b2 fixed to the limit value.
We then have a linear least square fitting with 2 parameters.

If :math:`b_1` or :math:`b_2` is fixed to :math:`\bar{b_1}` (resp. :math:`\bar{b_2}`), we use :math:`M(\theta) = M \cdot \theta + \gamma` in **(2)**, which leads to: 

.. math::
     \boxed{
          \theta =  (M^T N^{-1} M)^{-1} \cdot M^T N^{-1} (d - \gamma)
     }

b1 fixed
^^^^^^^^

From **(1)** , we obtain:

.. math:: 
     \theta = \left( \begin{array}{c} A_2 \\ b_2 \end{array} \right)
     \text{, }
     M = \left( \begin{array}{cc} 1 & X_c \\ \vdots & \vdots \\ 1 & X_i \\ \vdots & \vdots \end{array} \right)
     \text{, and }
     \gamma = \left( \begin{array}{c} \bar{b_1}(X_i - X_c) \\ \vdots \\ 0 \\ \vdots \end{array} \right)


b2 fixed
^^^^^^^^
From **(1)** , we obtain:

.. math:: 
     \theta = \left( \begin{array}{c} A_1 \\ b_1 \end{array} \right)
     \text{, }
     M = \left( \begin{array}{cc} 1 & X_i \\ \vdots & \vdots \\ 1 & X_c \\ \vdots & \vdots \end{array} \right)
     \text{, and }
     \gamma = \left( \begin{array}{c} 0 \\ \vdots \\ \bar{b_2}(X_i - X_c) \\ \vdots \end{array} \right)

b1 and b2 fixed
^^^^^^^^^^^^^^^
From **(1)** , we obtain:

.. math::
     \theta = \left( \begin{array}{c} A_1 \end{array} \right)
     \text{, }
     M = \left( \begin{array}{c} 1 \\ \vdots \end{array} \right)
     \text{, and }
     \gamma = \left( \begin{array}{c} \bar{b_1}X_i \\ \vdots \\ \bar{b_1}X_c + \bar{b_2}(X_i - X_c) \\ \vdots \end{array} \right)
   
Calculating coefs standard deviations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variances and covariances of A1, b1, b2 are the terms of :math:`M^{-1}`.
We use the approximation :math:`\text{Var}(a_1) = \text{Var}(\exp(A_1)) \approx a_1^2 var(A1)` 

For A2, based on :math:`A2 = A1 + (b1-b2) X_c` we find:

.. math::
     \text{Var}(A_2) = \text{Var}(A_1) + X_c^2 \left( \text{Var}(b_1) + \text{Var}(b_2) - 2 \text{Cov}(b_1,b_2) \right) + 2 X_c \left(\text{Cov}(A_1, b_1) - \text{Cov}(A_1, b_2)\right)

And using the same approximation than for :math:`a_1`, we deduce the formula :math:`\text{Var}(a_2)` from :math:`\text{Var}(A_2)`

Estimating continuum SNR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To calculate continuum SNR, we use :math:`\text{max}(\frac{a_1}{\sigma_{a_1}}, \frac{a_2}{\sigma_{a_2}})`

Propagating the noise in log scale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Knowing the in input standard deviation :math:`\sigma`` on each sample y, we need an expression of the standard deviation :math:`\sigma_{log}` of :math:`\log(y)`.
Using the transformation :math:`\ln(y+ \partial y) = \ln(y) + \frac{\partial y}{y} + o(\partial y^2)`,
we approximate :math:`\ln(y \pm \sigma) \approx \ln(y) \pm \frac{\sigma}{y}`

We will use:

.. math::
     \boxed{
          \sigma_{log} = \frac{\sigma}{y}
     }

Note that this has an impact on the calculations : the standard deviation to use in matrix :math:`N` is not the simple standard deviation of pixel i anymore, but it is dependent on the model flux, which we don't know before fitting. Therefore, we will compute the power law coefficients in two steps : first step without ponderation, and once a first estimate of the power law coefficients have been made, we compute a second step using this power law estimate to ponderate the standard deviation of the flux, and retrieve our final power law coefficients.


.. comment
     :math:`ln` is also applied to the x-axis. We could empirically compensate with :math:`\sigma_{loglog} = x \sigma_{log} = \frac{x}{y} \sigma`

     For the moment, it did not show a big difference in the tests so we will not apply this ponderation yet.






In practice
--------------


Simple case 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each redshift : 

#. Intialize flux curve : apply a blue shift on the input flux so lambda axis corresponds to lambda rest
#. If the number of unmasked samples is lower than thershold defined in parameter `nbSamplesMinForContinuumFit`, we force coefficients to zero and issue a warning.
#. Compute SNR-compliant pixels : pixels for which :math:`\frac{flux}{error} > \text{threshold}` with :math:`\text{threshold}` defined in parameter `continuumFit.nullThreshold`. If the number of valid pixel i.e. pixels which are both not masked (from input spectrum) and SNR-compliant is lower than the threshold defined in parameter `nbSamplesMinForContinuumFit`, we compute constant law : :math:`y = a` with :math:`a = \frac{\sum d_i \, w_i}{\sum w_i}` with :math:`w_i = \frac{1}{\sigma_i^2}`, taking into account all unmasked pixels (event the non-snr compliant ones). Igm and ism are forced to zero.
#. Compute the emitted curve: we create one flux curve per igm / ism correction.
#. For each of these emitted curves, compute the power law coefs, according to the equations above. If one of either sides of lambda cut does not have enough samples (compared to threshold `nbSamplesMinForContinuumFit`), the coefficients are computed using the side with enough samples: we use a simple least square fitting method. They are then extended to the side with too little samples before chi2 calculation. If there is no side with enough samples, coefficient are forced to zero, and a warning is issued.
#. At each power law computation, if :math:`a` is too small (i.e. smaller than :math:`\text{__DBL_MIN__}` i.e. the smallest positive value representable by a double on the machine), all coefficients are forced to zero, and their standard deviation to infinity.
#. Compute chi2 for each obtained power law coefficients, and find igm / ism indexes which minimize this chi2.

Power law
=========

.. role:: raw-latex(raw)
     :format: latex html


A fitting method to use for QSOs.
The idea is to fit the flux with two power laws (first power law for wavelength < 5400 Å, second one for wavelength > 5400 Å).
We make the fit using the least square fitting method on a log/log scale.


Applying least square fitting method on log/log scale
-----------------------------------------------------

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


Using matrix notation, minimze :math:`(d-M \cdot \theta)^T N^{-1} (d - M \cdot \theta)` with :math:`N` the covariance matrix,  a diagonal matrix with :math:`n_{i,i}=\sigma_i^2`.

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

Which can be reduced with the continuity constraint of the 2 power laws in (1) :

.. math::
     \theta = \left( \begin{array}{c} A_1 \\ b_1 \\ b_2 \end{array} \right)
     \text{, and }
     M = \left( \begin{array}{ccc} 1 & X_i & 0 \\ \vdots & \vdots & \vdots \\ 1 & X_c & -X_c + X_i \\ \vdots & \vdots & \vdots \end{array} \right)


Matrix calculations then allow to find an analytic solution to this equation.


Calculating coefs standard deviations
-------------------------------------

Variances and covariances of A1, b1, b2 are the terms of :math:`(M^T N^{-1} M)^{-1}`.
We use the approximation :math:`\text{Var}(a_1) = \text{Var}(\exp(A_1)) \approx a_1^2 var(A1)` 

For A2, based on :math:`A2 = A1 + (b1-b2) X_c` we find:

.. math::
     \text{Var}(A_2) = \text{Var}(A_1) + X_c^2 \left( \text{Var}(b_1) + \text{Var}(b_2) - 2 \text{Cov}(b_1,b_2) \right) + 2 X_c \left(\text{Cov}(A_1, b_1) - \text{Cov}(A_1, b_2)\right)

And using the same approximation than for :math:`a_1`, we deduce the formula :math:`\text{Var}(a_2)` from :math:`\text{Var}(A_2)`

Estimating continuum amplitude
------------------------------

To calculate continuum SNR, we use :math:`\text{max}(\frac{a_1}{\sigma_{a_1}}, \frac{a_2}{\sigma_{a_2}})`
If SNR is lwoer than threshold defined in parameters, we apply same behavior than in template fitting.


Behaviour for too little sample
-------------------------------

Two main cases of too little cases are possible:

* The total number of samples is too low : coefficient are forced to zero, and a warning is issued.
* One of either sides of lambda cut does not have enough samples:
    - The coefficients are calculated using the side with enough samples
    - They are extended to the side with too little samples before chi2 calculation


Calculating the noise
---------------------

Using the transformation :math:`\ln(y+ \partial y) = \ln(y) + \frac{\partial y}{y} + o(\partial y^2)`,
we approximate :math:`ln(y \pm \sigma) \approx \ln(y) \pm \frac{\sigma}{y}`

We will use:

.. math::
     \boxed{
          \sigma_{log} = \frac{\sigma}{y}
     }

:math:`ln` is also applied to the x-axis. We could empirically compensate with :math:`\sigma_{loglog} = x \sigma_{log} = \frac{x}{y} \sigma`

For the moment, it did not show a big difference in the tests so we will not apply this ponderation yet.







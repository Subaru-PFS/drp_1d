.. _linefitting:



Line Fitting
===================

Line fitting can be done and customized in a few different ways: different constraints can be applied to the lines, and different methods can be used to fit the lines amplitudes.

Lines can be fitted independently (linemodel free), jointly (lineratio), or independently but with some constraints (rules).

Different algorithms can be used to perform the line fitting (parameter `lineFittingMethod`) :

- individual
- svd
- hybrid
- lbfgsb

Note: in the following we name an **element** a set of lines physically linked together (for instance the OII doublet). This means that the ratios of all lines inside an element are fixed and only one amplitude parameter will be fitted for te whole element, whatever the selected fitting method. The status of lines to be linked or not is defined in the input line catalog (parameter `lineCatalog`).

Line profiles
--------------

Different line profiles are defined and set in the calibration. See :ref:`lineprofiles`.

Initializing lines supports
-----------------------------

Is a line visible ?
^^^^^^^^^^^^^^^^^^^^

For each line (knowing its central wavelength defined in the line catalog), for each redshift, we estimate the wavelength window possibly impacted by the line, named the **line support**.
The wavelength window is defined as :math:`\mu \pm ( \sigma \, n_{\sigma} + o)`, where :

- :math:`\mu` is the line wavelength in observed frame
- :math:`\sigma` is the line gaussian width. It is computed as the quadratic sum of the sigma of the LSF (for Line Spread Function which includes the source size in the case of slitless spectroscopy) and the sigma of the velocity dispersion :math:`\sigma = \sqrt{\sigma_{LSF}^2 + \sigma_{velocity}^2}`, with :math:`\sigma_{velocity} = \frac{v}{c} \, \mu`, :math:`c` being the speed of light and :math:`v` the observed object's velocity dispersion.
- :math:`n_{\sigma}` is defined in pararmeter `nSigmaSupport`
- :math:`o`  is the maximal allowed line offset, when the wavelength offset fitting is enabled.

Then we study the overlap between this computed line window and lambda range. The line is considered as visible only if :

- there are enough spectral samples inside the wavelength window defined above
- there is at least one sample near the line center. : noting :math:`d_{min}` the minimal distance between the line center and each sample in the window, we require :math:`d_{min} \leq d_{max}`, with :math:`d_{max} = n_{\sigma max} \, \sigma + o`, with :math:`n_{\sigma max}` defined in parameter `nSigmaMax`


.. _lya:

Fitting lya profile 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that the lya emission line can be fitted using two methods, chosen with parameter `lya.profile` :

- **asym** : fitting a profile which corresponds to a skew normal distribution (see :ref:`line profiles<lineprofiles>`). Three parameters are fitted: 

  - :math:`alpha`, the skewness parameter of the skewed gaussian
  - :math:`delta`, the line center offset, in wavelength, in Angstrom
  - :math:`sigma`, the line width.

 The fit consists in a simple brute force search algorithm on these three parameters, fitting best amplitude for each combination, minimizing the chi2. :math:`alpha` range is determined through parameters `asymProfile.min`, `asymProfile.max` and `asymProfile.step`.

- **igm** : the lya asymetry is reproduced using igm. Since lines below Lya (1216 Angstrom restframe) are impacted by the IGM, their profile can be substantially distorted. Thus we need to compute this distorted profile by fitting the IGM on these lines.

   - The lines possibly impacted by igm are inventoried i.e. the lines which have been considered as visible (see above), and which rest wavelength is below 1216 Å (lya).
   - We loop on possible igm profiles, fit lines amplitudes using least square fitting, and compute merit.
   - The igm minimizing the merit is selected as the "best" igm.

Here and below, fitting the lines is made through :ref:`least square fitting <leastsquare>` method.

Performing the fit
---------------------

Depending on the chosen line fitting method, the lines parameters are fitted differently:

- the individual fitting method is only able to fit independently the amplitude of each element
- the svd, hybrid and lbfgsb fitting method can fit jointly the elements amplitudes and in additon can fit a second degree polynomial under each element, and an offset in wavelength of each element (with parameter `lbdaOffsetFit`) 
- the lbfgsb fitting method can in addition fit the line width jointly with all the other free parameters

Here is a brief description of each method.

Individual
^^^^^^^^^^^^^^

Elements are fitted independently, using for each element a :ref:`least square fitting <leastsquare>` method.

Svd
^^^^^^^^^^

With this method, superposed line profiles on some of the pixels are taken into account by fitting all elements simultaneously, using a linear least square fitting method (using singular value decomposition of the design matrix). 
It is also possible to simultaneously fit a polynomial of degree 2, in order to take into account possible small scale continuum residuals (set with parameter `ampOffsetFit`).
Emission and absorption lines are fitted independently. A shift in wavelength of each element can also be fitted (parameter `lbdaOffsetFit`) by just doing a loop on different values of the offset.

Mathematical model
""""""""""""""""""

For each observed spectral sample :math:`i`, the data are modeled as a linear
combination of line profiles plus an optional polynomial offset:

.. math::

   d_i =
   \sum_{k=1}^{N_{\text{lines}}} a_k \, y_k(\lambda_i, z)
   \;+\;
   \sum_{p=0}^{d} b_p \, \lambda_i^p
   \;+\;
   \eta_i

where:

- :math:`d_i` is the observed (continuum-subtracted) flux at wavelength
  :math:`\lambda_i`,
- :math:`y_k(\lambda_i, z)` is the model value of the :math:`k`-th spectral line
  profile evaluated at :math:`\lambda_i` and redshift :math:`z`,
- :math:`a_k` are the line amplitudes to be fitted,
- :math:`b_p` are the coefficients of the polynomial offset of degree
  :math:`d`,
- :math:`\eta_i` is additive Gaussian noise

In matrix form, the model can be written as:

.. math::

   \mathbf{d} = X\,\mathbf{c} + \boldsymbol{\eta}

where:

- :math:`\mathbf{X}` is the design matrix whose columns are the line profiles
  :math:`y_k(\lambda_i, z)` and, if enabled, the polynomial basis functions
  :math:`1, \lambda_i, \lambda_i^2`,
- :math:`\mathbf{c}` is the parameter vector containing the line amplitudes
  :math:`a_k` and the polynomial coefficients :math:`b_p`,
- :math:`\boldsymbol{\eta}` is the noise vector.

The parameters are estimated by minimizing the weighted least-squares criterion:

.. math::

   \left(\mathbf{d} - \mathbf{X}\,\mathbf{c}\right)^{\mathrm{T}}
   \mathbf{W}
   \left(\mathbf{d} - \mathbf{X}\,\mathbf{c}\right)

with the weight matrix :math:`\mathbf{W} = \mathrm{diag}(1 / \sigma_i^2)`.

The solution of this linear least square minimization problem is given by:

.. math::

   \hat{\mathbf{c}} = \left[ \mathbf{X}^\mathrm{T}\mathbf{WX} \right]^{-1} \mathbf{X}^\mathrm{T}\mathbf{W} \mathbf{d}

and the covariance matrix of the estimated parameters :math:`\hat{\mathbf{c}}` is:

.. math::

   \text{COVAR}(\hat{\mathbf{c}}) = \left[ \mathbf{X}^\mathrm{T}\mathbf{WX} \right]^{-1}

The inversion involved in these computations is performed using the singular value decomposition of the matrix X.

In practice
"""""""""""""

- Builds the list of elements to fit
- Compare the number of degrees of freedom (number of elements amplitudes to fit + number of polynomial coefficients - 2) to the number of available spectral samples. If there are not enough samples, the fit is skipped, a warning is issued, and all amplitudes are set to NAN.
- For numerical conditioning, normalize flux and line profiles by the maximal absolute flux value over the fitting range.
- Design the matrix X containing the modeled line profiles (samples, lines)
- Remove columns of zeros i.e. null line prxofiles
- Use linear least square fitting to estimate amplitudes + polynomial coefficients (method `gsl_multifit_wlinear` from gsl library)
- If some amplitudes are negative, force them to zero and refit the positive ones.
- use the square root of the diagonal of the covariance matrix to get the uncertainty on the lines amplitude and the covariance matrix of the polynomial coefficient to get the uncertainty on the estimated continuum under the lines.

hybrid
^^^^^^^^^^^^^^

The hybrid fitter method combines the individual and svd methods. It consists in grouping all the line elements that are sharing pixels. Each group can then be fitted independently from the others using the above SVD method, or even the individual method when an element is alone in a group.

- Build the groups of independant elements
- loop on the groups:

  - if the group contains one single element without the polynomial use the individual fitting method
  - otherwise use the svd fitting method on the elements of the group. Note that if amplitude offsets are enabled (parameter `ampOffsetFit`), the hybrid method always uses the svd method.

- If Balmer improvement is enabled (parameter `improveBalmerFit`), then a joint fit of emission and absorption Balmer lines is performed. The solution minimizing the model RMS is kept.

Description of the Balmer improvement algorithm
""""""""""""""""""""""""""""""""""""""""""""""""""
- Enumerate predefined Balmer line pairs (emission + absorption).

- For each pair:

    - Check whether the pair is valid and physically meaningful.
    - Optionally include nearby related emission lines.
    - Test whether a joint refit is warranted : the observed absorption line width must be significantly broader than the emission line :math:`sigmaA \ge 2 \, sigmaE`
    - Attempt a local linear refit of amplitudes.
    - Keep the refit only if it improves the residuals.

.. _lbfgsb:

lbfgsb
^^^^^^^^^^^^^^

lbfgsb method works like the hybrid method, but uses a general optimizer with box constraint (lbfgsb algorithm) instead of the SVD to solve the non-linear least-square case by including the position and width of the line profiles as free parameters in addition to the amplitude.

Below, we explain only the part which replaces the svd method in the hybrid fitter. The rest is the same.

Like for the svd method, we enumerate the number of degrees of freedom, and check that the number of samples is enough to fit all the degrees of freedom. The degrees of freedom are :

- elements amplitudes
- the dispersion velocities (if enabled with parameter `velocityFit`) : one for emission lines, and one for absorption lines (if absorption and emission are mixed)
- the wavelength offset (if enabled with parameter `lbdaOffsetFit`)
- polynomials coefficients (if enabled with parameter `ampOffsetFit`).

Like in svd, there is a normalization step for all the parameters (flux, dispersion velocities, wavelength offsets), so that the ranges of variation are consistent.

We then define parameters bounds to constrain the solution in these defined boundaries.

- velocities are bounded between min and max values defined in parameters (`emVelocityFitMin` / `emVelocityFitMax`, `absVelocityFitMin`, `absVelocityFitMax`)
- wavelength offsets is bounded between min and max values defined using parameter `lbdaOffsetMin` / `lbdaOffsetMax`
- polynomial coefficients are left unconstrained
- An initial guess is computed for all parameters, using the svd fitter, with a lambda offset of 0, and velocities set to the parameters `velocityEmission` and `velocityAbsorption`.
- We then compute the resulting snr :math:`\frac{amplitude}{amplitudeError}` for each element. If the maximum snr is too low i.e. lower than 1, we stop there and keep this initial guess

Now that we have the initial guess, we can run the lbfgsb algorithm, minimizing on the chi2. 
The uncertainties on each parameter are computed using the square root of the diagonal of the approximated inverse Hessian given by the lbfgsb algorithm. 


Using a line ratio catalog
------------------------------

Rather than fitting the individual elements, we can choose to fix the ratios between the different lines amplitudes (activated with parameter `lineRatioType=tplRatio`). This limit the number of degrees of freedom of the model. Several lines ratios can be defined, they are defined in a line ratio catalog from the calibration directory.

In this case, the term **element** refers to the whole ratio lines. To be more precise, we define one element for all emission lines, and one element for all absorption lines (if activated). Therefore, only two amplitudes are fitted.


Applying rules
------------------------------

Another way to constrain the fitted lines is to apply some rules (activated with parameter `rules`) to the fitted amplitude.
The different defined rules are :

- `strongWeak`
- `balmerSingle`
- `ratioRange`

One or several of these rules can be jointly selected.


Strong weak
^^^^^^^^^^^^^^

Ensures that the all fitted amplitude of weak lines are lower than all fitted amplitudes of strong lines. If weak lines are too strong, their amplitudes are reduced. Weak / strong lines are defined in the line catalog.
Emission / absorption lines amplitudes are treated separately.


Balmer Single
^^^^^^^^^^^^^^

Ensures that the fitted amplitude of a weaker line does not exceed a fixed multiple of a stronger reference line : if the weaker line is too strong, its amplitude is reduced. 

The implemented balmerSingle rules are : 

- :math:`H_\beta \le \frac{1.1}{2.86} \times H_\alpha`
- :math:`H_\gamma \le 1.1 \times 0.47 \times H_\beta`
- :math:`H_\gamma \ge H_\delta \ge H_8 \ge H_9 \ge H_{10} \ge H_{11}` with 10% slack for emission lines
- :math:`H_\beta \ge H_\gamma \ge H_\delta \ge H_8 \ge H_9 \ge H_{10} \ge H_{11}` with 10% slack for absorption lines

Ratio range
^^^^^^^^^^^^^^

This method ensures that the fitted amplitudes of two lines
:math:`A` and :math:`B` satisfy an approximate ratio constraint:

.. math::

   \frac{A}{B} \approx R \quad \text{with} \quad R = m_{\mathrm{Coefficient}}

If the fitted amplitudes violate this ratio beyond tolerance, the method adjusts **both amplitudes together**, using uncertainty-weighted least squares.


Let : 

- :math:`R = m_{\mathrm{Coefficient}}`
- :math:`(a_1, \sigma_1)` the dominant line, and
- :math:`(a_2, \sigma_2)` the weaker line

The weights are defined as:

.. math::

   w_1 = \frac{1}{\sigma_1^2}, \quad
   w_2 = \frac{1}{\sigma_2^2 R^2}

The corrected amplitudes are then computed as:

.. math::

   a_{1,\mathrm{corr}} = \frac{a_1 w_1 + a_2 w_2 R}{w_1 + w_2}, \quad
   a_{2,\mathrm{corr}} = \frac{a_{1,\mathrm{corr}}}{R}

In our case, the ratio range is applied to lines `cIII1907_em` and `cIII1909_em` with a coefficient of 2.

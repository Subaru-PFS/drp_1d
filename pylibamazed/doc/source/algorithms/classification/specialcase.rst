Classification special case : fromSpectrum switch
=================================================

When lineModelSolve (LMS) continuum fit is initially in auto mode (parameter `continuumComponent` is `tplFitAuto`/`powerLawAuto`), 
if the fitted continuum in first pass is considered as bad the continuum fitting is discarded and the continuum 
is estimated by median filtering the input spectrum (equivalent to have the parameter `continuumComponent` set to `fromSpectrum`). 
In such a case the lineModelSolve evidence is computed using continuum subtracted data whereas the templateFittingSolve method, 
for instance used for the Star model is computed on the full data. 
This leads to an unfair comparison, since the fits are not based on the same input data. 
It will essentially favor models using continuum subtracted data with lineModelSolve.

Case one model with lineModelSolve continuumFit set to fromSpectrum
---------------------------------------------------------------------

To reduce this artifact when there is only one model using lineModelSolve in `fromSpectrum` mode 
(either by requesting that mode or by switching automatically to it) 
we use the continuum evidence calculated in the lineModelSolve first pass for the lineModelSolve model.

Note: since the continuum only evidence is used, it means the fit can be degraded by the presence of lines 
in the spectrum that are not in the model.

Case several models set to fromSpectrum
------------------------------------------------------------------------

We first classify the models using the continuum evidence. If two LMS models are side by side in the classification, we refine their evidence usign a ratio combinig their continuum evidence and their "full" evidence (i.e. the one obtained after switching to spectrumModel)

Example :

Noting : 
 - :math:`E_{TFS}` the evidence of spectrum model using template fitting solve

 - :math:`E{i}^c_{LMS}` the continuum evidence of spectrum model i using line model solve

 - :math:`E{i}^f_{LMS}` the "full" evidence of spectrum model i using line model solve

 - :math:`E{i}_{LMS}` the evidence of spectrum model i using line model solve we finally use for classification

If :

   .. math:: 

      E1^c_{LMS} = 3 > E2^c_{LMS} = 2 > E_{TFS} = 1 \\
      E1^f_{LMS} = 4 < E2^f_{LMS} = 6

We use the ratio :

   .. math:: 

      r = \frac{E1^c_{LMS} + E2^c_{LMS}}{E1^f_{LMS} + E2^f_{LMS}} \text{  and  } E{i}_{LMS}= r \times E{i}^f_{LMS}


And finally use the the evidences :

   .. math:: 

      E2_{LMS} = 3 > E1_{LMS} = 2 > E_{TFS} = 1
  



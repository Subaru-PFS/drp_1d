Classification
==============

We compare the marginalized evidence of every spectrum model. The spectrum model with the best evidence is the selected one.

Case lineModelSolve continuumFit switched to fromSpectrum
---------------------------------------------------------------------

If lineModelSolve (LMS) continuum fit is initially in auto mode (`tplFitAuto`/`powerLawAuto`), its continuum calculated in first pass is considered as bad, continuum fit is forced to `fromSpectrum`. Line model solve evidence is now calculated using modified spectrum : its continuum is filtered out. 
However, template fitting solve still calculates the evidence using the full spectrum.
This leads to artifically classifying more objects as models using line model solve.

To reduce this behaviour, we use the continuum evidence calculated in line model solve first pass.


Case several models switched to fromSpectrum
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
  



LineMeas solver
===============

For each model (Galaxy, Star, QSO...), once the redshifSolver stage has been processed, 
the spectral feature measurement can optionaly be triggered (lineMeasSolver stage) 
and is using the estimated redshift of the previous redshifSolver stage.

The processing is very similar to the :ref:`linemodelsolve` applied on one redshift, without going through the zPDF computation.

The goal here being to fit each line profile with all three parameters (amplitude, position and width), the fitting method should be
:ref:`lbfgsb` (the other line fitting methods are available but should not be used).

The handling of the continuum is also the same than for the redshifSolver :ref:`linemodelsolve` (see :ref:`continuum`), with all possible choices 
but the goal being to get a good fit of the lines, the prefered way is to:

- set the continuum as null (parameter `continuum` set as "noContinuum")
- activate the fitting of a polynomial under each line (parameter `ampOffsetFit` set at "True")
  
Then

- activate the velocity fit (parameter `velocityFit` set at "True")
- activate the lambda offset fit to enable the position of each line to be shiffted relative to the fixed estimated redshift (parameter `lbdaOffsetFit` set at "True")
- set the `lineRatioType` to "rules" and `rules` to "no" (ie not using fixed template ratios, nor line ratio rules)
 
When the preceding resdhift estimation solver is lineModelSolve the initial velocity dispersion (used to compute the line width) 
used is the one estimated from the redshiftSolver. In the other cases (templateFittingSolve), 
the initial velocity dispersion is determined by the parameters (`velocityEmission` and `velocityAbsorption`) 


Note:

- The lineMeas Solver can be requested for each model, even when the redshifSolver is not lineModelSolve.  
- The lineMeas stage, when requested, can be skipped if the model is not selected by the classification (whith parameters `lineMeasRunMode` set to "classif" instead of "all")
Classification
==============

The classification between the models Galaxy, Sar, Qso... is based on the Bayesian model selection.
The best selected model is the one having the largest probability :math:`p(\mathcal{M}=m | \mathbf{d})` 
where :math:`m` can be :math:`\text{'gal'}`, :math:`\text{'star'}`, :math:`\text{'qso'}` 
(depending on the predefined model category).

We have:

.. math::

   p(m | \mathbf{d}) = \frac{p(\mathbf{d} | m) p(m)}{p(\mathbf{d})}

where

- :math:`p(\mathbf{d} | m)` is the evidence of the model :math:`m` (see section :ref:`bayesian-framework` )
- :math:`p(m)` is the prior probability of the model
- :math:`p(\mathbf{d}) = \sum_i p(\mathbf{d} | m_i)` is the sum of the evidences of all the models

Currently no model prior can be given as input, thus a non-informative uniform prior is used:
:math:`p(m) = \frac{1}{n_{\text{models}}}`

Note: when using 2 pass processing for the redshift solvers of the different models,
there is the possibility to make the classification early after the completion of the first pass for each model 
(see parameter `secondPassAfterClassification`). 
Then the second pass is only performed for the classified model, the others keep their results from the first pass.

.. toctree::
   :maxdepth: 4
   :titlesonly:
   
   specialcase
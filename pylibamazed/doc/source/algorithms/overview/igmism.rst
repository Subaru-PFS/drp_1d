Igm ism
=======

The ISM and IGM flux correction curves allow for physically motivated variability in the model. The currently implemented ISM correction curve is based on :cite:`Calzetti2000`. The IGM extinction curve is based on a set of redshift-dependent tables from :cite:`Meiksin2005`.

For a model :math:`M` without IGM / ISM, which can be for instance a line profile or a continuum template, we note :

.. math::
    M_{\text{igm_ism}} = M(\lambda_{\text{rest}})\, IGM_{z,m}(\lambda_{\text{rest}})\, ISM_{\tau}(\lambda_{\text{rest}})


- IGM applies to the continuum and to emission lines which wavelengths are lower or equal to Lya (both are independent)
- ISM applies to the continuum and to emission lines (only in template ratio, because smooth variation does not modify profiles, but line ratios), and again, both are independent.
- Absorption lines are not affected because we fit a continuous absorption profile that already includes the ISM or IGM factor.


Igm
---

Meiksin curves are used to reproduce IGM extinction.
Different meiksin curves are defined on the redshift bin.
For each redshift bin, 7 extinction curves are defined: one "mean" curve, and +/- 0.5, 1, 1.5 :math:`\sigma`

.. image:: ../_static/images/meiksin-curves.png




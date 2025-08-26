Preliminary concepts
====================

.. todo revoir cette section, et son articulation avec la page client developer guidelines

In pylibamazed, we define several solvers organized in a pipeline. The created pipeline can process different object types (ex: galaxy, star, ...).

It is composed of severel stages:

- **init**: preceeds object types processing, common for all.

.. TODO ici clariifier, et mettre des vrais noms qui correspondent à quelque chose

- **object stages**: It contains ``redshift_solver``, ``linemeas_catalog_load``, ``linemeas_solver``, ``reliability_solver`` and ``sub_classif_solver``. These are the main parts of the pipeline.
- **classification** and **result_store_fill**: one or two stages that follow the object stages, common for all object types.


Pipeline definition
-------------------
All the parameters used to define the solver are in a dictionary. This dictionary is usually described in a json file named ``parameters.json``.

``spectrumModels`` define the names of the object types on which to launch object solvers, for example "galaxy", "star", "qso". Each object has its dedicated section, with the solvers to use, their methods and parametrization (see :doc:`../../parameters` page).

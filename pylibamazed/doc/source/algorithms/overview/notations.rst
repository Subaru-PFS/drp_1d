Notations
=========

Here we describe the main notation rules whoch apply accross this documentation. First we describe some simplifications we made in the notations. Then we list some symbols which have the same signification across the documentation (except explicit mentions).

Main rules
----------

- Bold letters correspond to vectors, lowercase letters to scalars, and uppercase letters to matrices.

- For simplicity, we note :math:`P(x) = P(\mathcal{X} = x)` the probability density that the variable :math:`\mathcal{X}` takes the value :math:`x`.

- We have the possibility to use one of several models ("galaxy", "star", "QSO", ...). Noting :math:`\mathcal{M}` the model variable and :math:`m` its value, we write :math:`p(\mathcal{X}=x \mid \mathcal{M}=m) = p(x \mid m) = p(x \mid \text{'gal'})`, for clarity.


Some symbols
------------

- :math:`a(z, \boldsymbol{\theta})` : scale parameter
- :math:`\chi^2` : least-squares metric
- :math:`z` : redshift
- :math:`\mathbf{y}(z, \boldsymbol{\theta})` : modeled flux (without amplitude)
- :math:`\mathbf{d}` : observed data
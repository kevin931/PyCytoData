######################
Installation Guide
######################

You can easily install ``PyCytoDR`` from either ``conda`` or ``PyPI``. We recommend
``conda`` for its easy environment management, but you can of course use ``PyPI``
if you prefer. Just follow the instructions below and you are good to go!

---------

***********
Conda
***********

You can install our package on ``conda``:

.. code-block::

    conda install pycytodata -c kevin931 -c bioconda

Our ``conda`` package is published `here <https://anaconda.org/kevin931/pycytodata>`_.

----------------

***********
PyPI
***********

You can also install our package on from ``PyPI``:

.. code-block::

    pip install PyCytoData

Our ``PyPI`` package is published `on this page <https://pypi.org/project/PyCytoData/>`_.

----------------

*************
Dependencies
*************

We try to make this package as lightweight as possible. Thus, as of now, we only depend on three packages:

- fcsparser
- pandas
- numpy

``fcsparser`` and ``pandas`` are used for handling benchmark datasets and their associated ``fcs`` files.
We use ``numpy`` for all computations and preprocessing.

Optional Dependency
--------------------

If you wish to perform dimension reduction (DR), you can install ``CytofDR`` and have it integrated into
this package. We know that DR is a pretty common workflow, but we didn't make it mandatory because
``CytofDR`` depends on a large number of DR packages, which makes the dependency unnecessarily complex
if one does not wish to perform DR. To install ``CytofDR``, you can run the following:

.. code-block:: shell

    pip install CytofDR

Or, with ``conda``:

.. code-block:: shell

    conda install -c kevin931 cytofdr -c conda-forge -c bioconda

For detailed documentation on ``CytofDR`` installation, which can get quite tricky if you want to use
its own optional dependencies such as ``SAUCIE`` and ``GrandPrix``, you can visit its
`Installation Guide <https://cytofdr.readthedocs.io/en/latest/installation.html>`_.
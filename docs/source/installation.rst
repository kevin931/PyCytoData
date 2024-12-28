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
- numpy >= 1.20

``fcsparser`` and ``pandas`` are used for handling benchmark datasets and their associated ``fcs`` files.
We use ``numpy`` for all computations and preprocessing.


Python Versions
------------------

We currently are tested on the following versions of Python:

- 3.9
- 3.10
- 3.11
- 3.12

However, we try to make ``PyCytoData`` as widely supported as possible. In fact, we in theory support
any Python version all the way back to `3.7`.

.. note:: Newer Version of Python
    We in theory should support newer versions of Python as well, such as Python `3.13`. However, one of the
    optional dependencies `CytofDR` below depends on `numba`, which does not support `3.13`. Therefore, this
    version does not currently pass our CI pipeline. If you do not plan on using `CytofDR`, feel free to use
    the newest version of Python available.


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
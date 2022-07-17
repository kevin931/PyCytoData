##########
Releases
##########

Here we include our release notes for past releases in sequential order.

------------------

********
v0.1.0
********

This is the first official release of ``PyCytoData`` with LTS.

Bug Fixes
-----------

- Fix the default cofactor of ``PyCytoData.preprocess``
- Fix lineage channels for ``levine32`` and ``samusik``

Changes and New Features
--------------------------

- ``PyCytoData.run_dr_methods`` now runs on lineage channels only
- Add tutorials and references for documentation
- Improve docstrings and general documentations
- Add preprocessing options for loading benchmark datasets
- Relax ``numpy`` dependency to >=1.20
- Add descriptors for PyPI releases


Deprecations
----------------

- Permanently remove ``FileIO.make_dir`` function for safety reasons


********
v0.0.1
********

- This is the first official prerelease of the ``PyCytoData`` package.
- We have proper support for the following workflows, including:
    - Downloading data
    - Using PyCytoData as CyTOF data analysis pipeline
    - FileIO
    - CyTOF DR Integration
- Releases on PyPI and conda

.. warning::

    There is a potential issue of compatibility with ``CytofDR`` on ``conda``. If a problem occurs, try
    using pip instead.
============
Changelog
============

Here are the most recent releases and changes. Currently, we're still under developmet.
Therefore, we don't have any official releases. However, check out our git history to see what we're
doing!

-------------------

*************************
Latest Release: v0.1.0
*************************

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

-------------------------

.. toctree::
    :maxdepth: 1

    releases
    recent

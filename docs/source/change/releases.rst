##########
Releases
##########

Here we include our release notes for past releases in sequential order.

------------------

********
v0.1.2
********

This is a minor release that fixes a critical bug that affects all previous releases.
Update is strongly recommended.

Bug Fixes
-----------

- Fixed a critical issue with ``PyCytoData.load_dataset()`` not working due to old web source down (`#3 <https://github.com/kevin931/PyCytoData/issues/3>`_)
    - This issue stems from upstream web source, thus affecting all previous releases.
    - This update fixes the issue with a different implementation.
    - There is no implementation change.


Changes and New Features
--------------------------

- No new feature added.

-------------------


********
v0.1.1
********

This is a minor release with various bug fixes and documentation improvents.

Bug Fixes
-----------

- Fixed a potential issue with loading benchmark datasets' samples out of order (This behavior is not guaranteed, but a implementation detail)
- Fixed an issue with ``DataLoader.load_dataset`` not recognizing downloaded datasets
- Fixed an issue with ``preprocess.bead_normalization`` having uninitialized array (#9)


Changes and New Features
--------------------------

- No new feature added
- Improved documentations with streamlined front page and updated links
- Added docstrings for `+` and `+=` operators
- Updated references for CytofDR paper publication in Nature Communications

-------------------

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

----------------------

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
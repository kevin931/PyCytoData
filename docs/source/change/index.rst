============
Changelog
============

Here are the most recent releases and changes. Currently, we're still under developmet.
Therefore, we don't have any official releases. However, check out our git history to see what we're
doing!

-------------------

*************************
Latest Release: v0.1.1
*************************

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

-------------------------

.. toctree::
    :maxdepth: 1

    releases
    recent

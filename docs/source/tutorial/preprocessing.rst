#########################
CyTOF Data Preprocessing
#########################

For CyTOF datasets, we oftentimes need to preprocess the expression. We currently support
the follow preprocssing steps:

* Arcsinh transformation
* Gating:
    * Debris Removal
    * Find intact cells
    * Find live cells
    * Center, offset, and residual outlier removal
* Bead normalization

.. note:: Currently, we do not support any builtin cross batch normalization.

We break this tutorial down to two sections. The first focuses on our builtin benchmark
datasets, which requires little normalization. The second part details the full API
and the details of each step.

---------------------

************************
Benchmark Datasets
************************

The only preprocessing step for all three benchmark datasets, ``levine13``,
``levine32``, and ``samusik`` is arcinh transformation. Typically, we use
a cofactor of 5 for such a transformation. You can let the ``DataLoader.load_dataset``
method perform the transformation automatically upon loading the
datasets:

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset="levine13", preprocess=True)

This will automatically apply the arcsinh transformation to the expression matrix.
You can use the data for downstream analyses if needed. The transformation is
applied to the lineage channels only.


If you wish to preprocess later, you can call the ``preprocess`` method
separately as well:

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset="levine13", preprocess=False)
    >>> exprs.preprocess(arcsinh=True, cofactor=5)

This will have the same effect as the previous snippet, but you have the choice
of applying the transformation whenever you wish or with different cofactors. 

---------------------------------------------

****************************
Full Preprocessing API
****************************


To access the API, you can use the ``preprocess`` method on your ``PyCytoData``
method:

.. code-block:: python

    >>> exprs.preprocess(arcsinh=True,
    ...                  gate_debris_removal=True,
    ...                  gate_intact_cells=True,
    ...                  gate_live_cells=True,
    ...                  gate_center_offset_residual=True,
    ...                  bead_normalization=True)
    Runinng Arcsinh transformation...
    Runinng debris remvoal...
    Runinng gating intact cells...
    Runinng gating live cells...
    Runinng gating Center, Offset, and Residual...
    Runinng bead normalization...

Each step can be optionally included. By default, the method automatically
attempts to resolve the necessary instrument channels. However, in case
this step fails, some manual intervention is needed.


Channel Names
-----------------

The preprocessing pipeline some specific channels to work properly. If none
are specified, the ``preprocess`` method tries to guess how such channels are
named. If ``auto_channels`` is set to ``False``, then specifying instrument
channles is required for gating and bead normaliztion; otherwise, any specified
channels will override the ``auto_channels`` function.

Below, we have a list of channels that users may wish to specify for preprocessing:

* ``lineage_channels``: These are protein marker channels used for analyses and arcsinh transformation (e.g. CD4).
* ``bead_channels``: These are the added beads to track time decay. Typically there are multiple bead channels (e.g. Bead1)
* ``time_channel``: The time of each event, typically named "Time".
* ``cor_channel``: This consists of "Center", "Offset", and "Residual".
* ``dead_channel``: Typically named "Dead" or "Live" to indicate cell status.
* ``DNA_channels``: Typically consisting of "DNA1" and "DNA2".

We assume that these channels are present and conventionally named. If they're 
not present, it can be an indication that some preprocessing has been done already.
As for arcsinh transformation, we suggest users plot the distributions of each channel
and verify whether transformation is needed.

Arcsinh Transformation
------------------------

Arcsinh transformation is often used to transform the channels. The transformation has
the following formula:

.. math::

    arcinh(\frac{\cdot}{cofactor})

Typically, people set :math:`cofactor=5` as per convention. If you prefer, you can
change it for your own data analyses.

As an example:

.. code-block:: python

    >>> exprs.preprocess(arcsinh=True,
    ...                  cofactor=2)
    Runinng Arcsinh transformation...


Gating: Debris Removal
------------------------

This is the first step in gating in which we use the bead channels and remove any
cells that are three standard deviations above the mean for each channel.

The ``bead_channels`` parameter must be given without ``auto_channels``. Otherwise,
the bead channels must start with "bead" (case insensitvie) for the method to
automatically detect channels.


.. code-block:: python

    >>> exprs.preprocess(gate_debris_removal=True,
    ...                  bead_channels = ["Bead1", "Bead2", "Bead3", "Bead4"])
    Runinng debris remvoal...


Gating: Finding Intact Cells
---------------------------------

This step uses the DNA channels to gate for intact cells. Specifically, it
trims cells with DNA greater than or smaller than ``cutoff_DNA_sd`` times
the standard deviation of the channels. By defaults, it preserves cells
within two standard deviations of the mean. The DNA channel names must
contain "DNA" or be provided specifically.

.. code-block:: python

    >>> exprs.preprocess(gate_intact_cells=True,
    ...                  DNA_channels = ["DNA1", "DNA2"],
    ...                  cutoff_DNA_sd = 2)
    Runinng gating intact cells...

Gating: Finding Live Cells
---------------------------------

This step uses the Dead channel to gate for live cells. Specifically, it
trims cells from the top percentile of the dead channel as specified by
``cutoff_quantile``. By default, the top 3rd percentile will be trimmed.
The DNA channel name must contain "dead" or be provided specifically.

.. code-block:: python

    >>> exprs.preprocess(gate_live_cells=True,
    ...                  dead_channel=["Dead"],
    ...                  dead_cutoff_quantile=0.03)
    Runinng gating live cells...


Gating: Center, Offset, and Residuals
----------------------------------------

This step gates cells using the Center, Offset, and Residual channels.
Specifically, it trims cells from the top and bottom percentile of the
three channels  as specified by ``cutoff_quantile``. By default, the top
and bottom 3rd percentile will be trimmed. The channels must be named as
such given here or be provided.

.. code-block:: python

    >>> exprs.preprocess(gate_center_offset_residual=True,
    ...                  cor_channels = ["Center", "Offset", "Residual"],
    ...                  cor_cutoff_quantile = 0.03)
    Runinng gating Center, Offset, and Residual...


Bead Normalization
------------------------

Our bead normalization algorithm uses the bead channels and the time
channel to correct signal decay. We uses a two-step process:

1. We remove cells whose bead sigals are in the bottom 5th quantile.
2. We perform the transformation using the most correlated bead channels.

This algorithm is developed and implemented in house. To perform bead
normalization:

.. code-block:: python

    >>> exprs.preprocess(bead_normalization=True,
    ...                  bead_channels = ["Bead1", "Bead2", "Bead3", "Bead4"],
    ...                  time_channel = ["Time"])
    Runinng bead normalization...                
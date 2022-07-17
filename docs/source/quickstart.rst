####################
Quickstart Guide
####################

``PyCytoData`` is a package that allows you to manage, download, and work with your
CyTOF datasets. It aims to create a unified interface, like that of ``Seurat`` in the
single cell universe. This guide walks you through the gist of the package, and with this,
you will be able to use this package without much of an issue at all. For more detailed
tutorials, check out the **Tutorial** section of the docmentation.

-------------------------

****************************
Loading Datasets
****************************

In this package, we provide a nice interface for you to load CyTOF
datasets. We offer two modes: you can either bring your own dataset (BYOD)
or you can use one of the existing benchmark datasets. The former allows
for the most flexibity whereas the latter is great for comparing your
results with thoese who have already walked this path.

Benchmark Datasets
--------------------

Currently, we support three benchmark datasets: 

=============== ==================
Dataset           Literal
--------------- ------------------
Levine 13-dim     ``levine13``
Levine 32-dim     ``levine32``
Samusik           ``samusik``
=============== ==================

As you may recognize, they are famous datasets that have been used numerous
times in the past. In fact, we use the `HDCytoData <https://github.com/lmweber/HDCytoData>`_
interface for downloading such data! To load these datasets, all you need to do
is the following:

.. code-block:: python

    >>> from PyCytoData import DataLoader

    >>> exprs = DataLoader.load_dataset(dataset = "levine13")
    Would you like to download levine13? [y/n]y
    Download in progress...
    This may take quite a while, go grab a coffee or cytomulate it!
    >>> type(exprs)
    <class 'PyCytoData.data.PyCytoData'>
    
Now you have a ``PyCytoData`` object to work with. The good news is that we cache
all datasets once you have downloaded them, meaning that there is no need to
download them again and again, which will save you time and bandwidth. The next
time you load the same dataset, you will no longer see the prompt:

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset = "levine13")

And of course, there are a few customization options you can have. For ``levine32``
and ``samusik`` which have more than one sample, you can choose to load specific
samples instead of loading them all (this will save you some RAM):

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset = "levine32", sample = ["AML08"])

or in the case of ``samusik``:

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset = "samusik", sample = ["01", "05"])
    >>> exprs.n_samples
    2

If you load multiple samples, they will be combined into one ``PyCytoData``, but sample
indices are preserved for you to subset and distinguish later one if you prefer.

.. note::
    
    The only preprocessing that you need to do is ``arcsinh`` transformation. All other steps have been performed.

Bring Your Own Dataset (BYOD)
------------------------------

Just as you can load your benchmark datasets easily, we also allow you to use your own
CyTOF datasets that you like. This offers the best flexibility. To do so, we have a
``FileIO`` class at your disposal:

.. code-block:: python

    >>> from PyCytoData import FileIO
    >>> exprs = FileIO.load_expression("Your_File_Path", col_names = True)

This returns a ``PyCytoData`` object! And you can access your expression matrix and
chennels names with the following attributes:

.. code-block:: python

    >>> exprs.expression_matrix
    >>> exprs.channels

And of course, if you don't have the first row as channel names, you can turn the option off:

.. code-block:: python

    >>> exprs = FileIO.load_expression("Your_File_Path", col_names = False)

In this case, no channel names will be stored. For more in-depth guide on IO and all its
functionalities, please head to the tutorials section and read the
`IO Guide <https://pycytodata.readthedocs.io/en/latest/tutorial/fileio.html>`_.


---------------------------

***************************
The ``PyCytoData`` Object
***************************

As you have seen in the previous section, the ``FileIO.load_expression`` method returns a
``PyCytoData`` object instead of an array. This is intentional: we want to group things
together. The ``PyCytoData`` object is able to store not only the expression matrix, but
also cell types, sample indices, and other metadata! Furthermore, it automatically checks
for errors when you manipulate these metadata. This makes it much less likely that things
go sideways when you work with your CyTOF data in your experiment. This section shows you
a little bit on how this works.

Accessing Attributes
----------------------

This is easy and pythonic:

.. code-block:: python

    >>> exprs.expression_matrix
    >>> exprs.channels
    >>> exprs.cell_types
    >>> exprs.sample_index

And the attributes are self-explanatory as well! By the same token, you can set these
attributes yourself! For example, when you load an expression matrix as a ``PyCytoData``
object, there are no cell types. You can set them accordingly:

.. code-block:: python

    >>> exprs.cell_types = cell_types

The setter method will ensure that dimension matches.

Metadata
---------

The object automatically computes a few metadata and they are automatically updated as well:

.. code-block:: python

    >>> exprs.n_cells
    >>> exprs.n_cell_types
    >>> exprs.n_samples
    >>> exprs.n_channels

These are implemented most for convenience and error checking! You don't have to work with
arrays' shape any more: you can simply refer to these dimensions by name!

Operations
-----------

You can not only store your data with ``PyCytoData``, but you can also do things with them.
You can preprocess your data and then run DR with the same object with the following verbs:

- ``preprocess()``
- ``run_dr_methods()``

Both of them will be further documented in the tutorials section.


Create Your ``PyCytoData`` Object
----------------------------------

The constructor is very easy to use:

.. code-block:: python

    >>> from PyCytoData import PyCytoData
    >>> exprs = PyCytoData(expression_matrix = expression_matrix,
    ...                    channels = channels,
    ...                    cell_types = cell_types,
    ...                    sample_index = sample_index,
    ...                    lineage_channels = lineage_channels)

All the parameters are self-explanatory as well! The only thing that you may be
unfamiliar with is ``lineage_channels``, which delineates actual lineage channels
from other instrument channels, such as Bead and time channel.


------------------------------

*********************
Preprocessing
*********************

We offer a full suite of preprocessing workflows at your disposal. The easiest way
is simply perform it on your ``PyCytoData`` object:

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

These are the six steps if you choose to do everything, but you can of course pick and choose.
It also depends the dataset you have: if your dataset doesn't have a lot of
instrument channels, it's likely been processed already! We detect these channels
automatically. For more details on each preprocessing step, go look at our
`CyTOF Data Preprocessing <https://pycytodata.readthedocs.io/en/latest/tutorial/preprocessing.html>`_ page.

------------------------------

**************************
Integration with CytofDR
**************************

The good news is that ``PyCytoData`` supports a ``CytofDR`` interface as
an optional extension of this package. After loading in your data and
performing all your necessary preprocessing steps, you can run DR methods
by simply calling a wrapper:

.. code-block:: python

    >>> exprs.run_dr_methods(methods = ["PCA", "UMAP", "ICA"])
    Running PCA
    Running ICA
    Running UMAP
    >>> type(exprs.reductions)
    <class 'CytofDR.dr.Reductions'>

And then you can perform any downstream DR workflows supported by ``CytofDR``.
Of course, if you're aware of the ``run_dr_methods`` methods from
``CytofDR``, you know that this is the "easy" mode. For more advanced usage,
you can set the ``exprs.reductions`` attribute with a ``Reductions`` object.
More information on the latter can be found at length on the
`CytofDR Documentation <https://cytofdr.readthedocs.io/en/latest/index.html>`_
website. 


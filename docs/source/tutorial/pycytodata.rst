###########################
The PyCytoData Object
###########################

The ``PyCytoData`` object is the workhorse that carries you through every step
of the way! If you're familiar with ``seurat`` and its signature workflow, you
should feel more or less at ease here: we aim to use one single object for
all your needs! All the actions will be as verbs whereas the metadata can be
accessed via attributes. 

This tutorial walks you through the details of the class and how you can best
utilize all its neat features beyond the most common usage as shown in the
quickstart guide.

---------------------------------

*******************************
Creating ``PyCytoData`` Object
*******************************

In general, there are three ways that you will have a ``PyCytoData`` object:

1. You load a benchmark dataset using the ``DataLoader.load_dataset`` method.
2. You load an expression matrix with the ``FileIO.load_expression`` method.
3. You create your own object using the constructor.

While the first two applications are much more common, we will walk you through
each of the three methods.

Loading Benchmark Dataset
----------------------------

This is by far the easiest and the most common method: 

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset="levine13")

The ``exprs`` object here is a ``PyCytoData`` object. The details of this method works
with benchmark datasets are included in the
`FileIO and Datasets <https://pycytodata.readthedocs.io/en/latest/tutorial/fileio.html>`_ section.

Loading Expression Matrix
---------------------------

If you have an existing expression matrix stored in plain text (not ``fcs``, which
is currently unsupported), you can easily load it into a ``PyCytoData`` object as
well:

.. code-block:: python

    >>> from PyCytoData import FileIO
    >>> exprs = FileIO.load_expression(files="<path>", col_names=True, delim="\t")

Again, ``exprs`` is a ``PyCytoData`` object, and the ``load_expression`` method has
more features, which are documented thoroughly in the 
`FileIO and Datasets <https://pycytodata.readthedocs.io/en/latest/tutorial/fileio.html>`_
section as well.

Using ``PyCytoData`` Constructor
------------------------------------

If you already have your data as arrays, you can of course construct an object:

.. code-block:: python

    >>> from PyCytoData import PyCytoData
    >>> exprs = PyCytoData(expression_matrix=expression_matrix,
    ...                    channels=channels,
    ...                    cell_types=cell_types,
    ...                    sample_index=sample_index,
    ...                    lineage_channels=lineage_channels)

Here, you can supply all your information on your CyTOF data. The only mandatory
component is the ``expression_matrix``. All parameters should be ``ArrayLike``,
but they don't have to be a ``numpy`` array. For the expression matrix, it naturally
should consist of numeric values, whereas others are expected to be strings (or maybe
integers).

In some cases, your data have some channels, such as beads, Time, or other instrument
channels. In this case, you can specify your ``lineage_channels`` to indicate which
are actually protein channels instead of other supplementary so that ``PyCytoData``
can process downstream analyses using appropriate channels accordingly.

-------------------------

********************************
Object Attributes and Metadata
********************************

One of the beauties of working with a ``PyCytoData`` object is that it offers a whole
set of builtin metadata and attributes you can work with. Some of them are automatically
added to the object whereas others you will need to supply. In any case, they will allow
you to get any information you need from your dataset.

Metadata
-----------

The metadata are basically statistics that are calculated from your dataset. You can
easily get the number of cells, channels, cell types, etc. Here, we will use a builtin
benchmark dataset as an example:

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset="levine32")

After loading ``levine32``, you can look at some of its statistics:

.. code-block:: python

    >>> exprs.n_cells
    265627
    >>> exprs.n_cell_types
    15
    >>> exprs.n_samples
    2
    >>> exprs.n_channels
    39

.. note:: There are 39 channels for ``levine32`` because there are instrument channels.

The metadata are automatically updated with changes and operations to the expression
matrix, such as subsetting, concatenating, etc. You can always be sure that the metadata
are up to date. 


Attribute: ``channels``
-------------------------

This attribute stores all the channel names. If no channel names are given during
instantiation, names are automatically generated using the ``"Channel" + int`` convention.
Namely, the first channel will be "Channel0" and all the rest are automatically numbers.
If channels are given or built into benchmark datasets, they're stored as is in this
channel.


Attribute: ``lineage_channels``
---------------------------------

This attribute denotes all the lineage channels, which are a subset of the ``channels``.
Lineage channels are protein channels typically used for analyses. They are stored to
differentiate from instrument channels, which are not used for analyses and transformations


Attribute: ``cell_types``
--------------------------

This attribute stores the cell types for each cell. If no cell types are given or available,
an array of ``None`` is automatically created and stored.

Attribute: ``sample_index``
---------------------------

This attribute stores the sample indices for each cell. They are stroed as strings or integers
within an array. If no sample information is available, all cells are assumed to be from the
same sample and indices of `0` are given.


Attribute: ``reductions``
--------------------------

This stores a ``Reductions`` object for dimension reduction from ``CytofDR`` package. By default,
this is ``None``. If the ``run_dr_methods`` method is called, the results will be automatically
stored. Users can supply and set their own ``Reductions`` object as well.


----------------------------------


***************************
Subsetting and Indexing
***************************

The ``PyCytoData`` object supports both the standard bracket notation, which is mostly consitent
with ``numpy`` indexing, and a custom ``subset`` method to subset by metadata, such as channels,
cell types, and samples. Here is a tutorial on both.

Indexing Using Brackets
-------------------------

We have implemented many of the features from ``numpy``'s behavior. Here are a list of allowable
types in the brackets (which are passed into the magic ``__getitem__`` method):

- slice
- List[int]
- np.ndarray

Indexing is performed on the expression matrix. A new object will be created and returned along
with the appropriate metadata. The indices can be one-dimensional or two-dimensional, allowing you
to index by the zeroth axis without having to specify the other axis. For example:

.. code-block:: python

    >>> new_exprs = exprs[:10]
    >>> new_exprs.n_cells
    10
    >>> new_exprs.n_channels
    39

As expected, this indexes the first 10 cells from the original object. This is equivalent to the
following:

.. code-block:: python

    >>> new_exprs = exprs[:10, :]

And of course, you can index both rows and columns: 

.. code-block:: python

    >>> new_exprs = exprs[5:100, [2,3,4]]

which will index the 5th to the 99th cells along with the 3 given channels. You can also use
``numpy`` arrays to index, which will enable you do something such as:

.. code-block:: python

    >>> new_exprs = exprs[5:100, np.isin(exprs.channels, ["CD3", "CD4"])]

As shown, you can mix and match indices types.

We do not support higher dimensional indexing because we assume that expression matrices are
two dimensional. Higher dimensional arrays will cause confusions.

.. note::

    We do not support indexing with integers such as ``exprs[5]``. This is to avoid the complexities
    introduced. If you wish to index a single obervation, use ``exprs[[5]]`` or ``exprs[5:6]``.


Subsetting
--------------

Instead of using indices and arrays, you can specifically subset based on metadata using the
``subset`` method provided. You can subset based on the following:

- channels
- sample
- cell_types

This usage is quite common and will save you a few seconds from using ``np.where`` or ``np.isin``.
To get started, you can simply do:

.. code-block:: python

    >>> new_exprs = exprs.subset(channels=['CD13(Er168)Di', 'CD3(Er170)Di'], sample=[0], cell_types=["Pro B Cells"], in_place=False)
    >>> print(new_exprs)
    A 'PyCytoData' object with 10 cells, 39 channels, 1 cell types, and 1 samples at 0x7fbf5eb0db80.

This will subset the dataset with:

- ``CD13(Er168)Di``, ``CD3(Er170)Di`` channels
- The 0th sample
- Pro B cells

You can pick and choose which of the metadata to subset. By default, the operation is done in
place. However, if you wish a new object to be returned, you can set ``in_place=False``.

Optionally, you can also negate the selection by setting ``not_in=True``, which will exclude the
given channels, sample, or cell types.


Subset Channels by Name
-------------------------

The standard ``subset`` method returns a ``PyCytoData`` object. However, if you wish to get an
array with the specified channels, you can subset using the ``get_channel_expression`` method:

.. code-block:: python

    >>> expressions, channels = exprs.get_channel_expressions(['CD13(Er168)Di', 'CD3(Er170)Di'])
    >>> expressions
    array([[ 1.92805946e-01, -1.63007990e-01],
    ...    [ 1.01540089e+01, -2.17397958e-01],
    ...    [ 1.07605422e+00,  1.63160920e+00],
    ...    ...,
    ...    [-8.74943915e-04, -1.50887862e-01],
    ...    [ 7.03829479e+00, -7.69027993e-02],
    ...    [ 1.62252128e+00, -2.47358829e-01]])
    >>> channels
    array(['CD13(Er168)Di', 'CD3(Er170)Di'], dtype='<U18')

This will save a step of getting the expression from an object.

-----------------------

*******************
Adding New Samples
*******************

If you have an existing ``PyCytoData`` object and you have a new sample, you can add your new
sample to your existing object. This can be easily achieved by the following:

.. code-block:: python

    >>> exprs = exprs.add_sample(expression_matrix=expression_matrix, sample_index = sample_index)

where ``sample_index`` is an array-like object with sample names, and ``expression_matrix`` is
just a matrix.

-----------------------------

********************
Magic Methods
********************

``PyCytoData`` implements a number of magic methods for convenience. Here, we detail the usage
of such methods.


Method: ``len``
-------------------

This returns the number of cells of the object, which is equivalent to using ``n_cells``:

.. code-block:: python

    >>> len(exprs)
    265627


Method: ``print``
--------------------

This method prints out the general information of the object with information on the number of
cells, channels, cell types, and samples along with its memory address.

.. code-block:: python

    >>> print(exprs)
    A 'PyCytoData' object with 265627 cells, 39 channels, 15 cell types, and 2 samples at 0x7fbf5ea09f10.


Method: ``+`` and ``+=``
-----------------------------

These two methods are used to concatenate two ``PyCytoData`` objects. The ``+`` operator returns a new
object:

.. code-block:: python

    >>> exprs1 = exprs = DataLoader.load_dataset(dataset="levine32", sample="AML08")
    >>> exprs2 = exprs = DataLoader.load_dataset(dataset="levine32", sample="AML09")
    >>> exprs = exprs1 + exprs2
    >>> print(exprs)
    A 'PyCytoData' object with 265627 cells, 39 channels, 15 cell types, and 2 samples at 0x7fbf313f7190.

Or if you prefer to concatenate a new object to an existing one, you can use the ``+=`` operator:

.. code-block:: python

    >>> exprs1 += exprs2
    >>> print(exprs1)
    A 'PyCytoData' object with 265627 cells, 39 channels, 15 cell types, and 2 samples at 0x7fbf5ea6f100.
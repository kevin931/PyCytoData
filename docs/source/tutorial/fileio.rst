###########################
File IO and Datasets
###########################

``PyCytoData`` provides an efficient interface for loading and working with
benchmark datasets. We offer the flexibility of both off-the-shelf benchmark
datasets and the ability loading your own dataset. Here, we detail how to
use each mode.

----------------------------

********************
Benchmark Datasets
********************

As stated in the `Quickstart Guide <https://pycytodata.readthedocs.io/en/latest/quickstart.html>`_,
we currently support three benchmark datasets: 

=============== ==================
Dataset           Literal
--------------- ------------------
Levine 13-dim     ``levine13``
Levine 32-dim     ``levine32``
Samusik           ``samusik``
=============== ==================

They are provided by the open-source `HDCytoData <https://github.com/lmweber/HDCytoData>`
implementation in R. We provide an easy interface for not only downloading
but also downstream analyses of such datasets.

Download Datasets
-------------------

You can load the dataset easily with the following code:

.. code-block:: python

    >>> from PyCytoData import DataLoader

    >>> exprs = DataLoader.load_dataset(dataset = "levine13")
    Would you like to download levine13? [y/n]y
    Download in progress...
    This may take quite a while, go grab a coffee or cytomulate it!
    

The resulting ``exprs`` is a ``PyCytoData`` object. You will have to do this
for only the first time: we cache all datasets once you have downloaded them.
The next time you load the same dataset, you will no longer see the prompt:

.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "levine13")

This will allow us to save some bandwidth and time. If you wish to download
datasets without being prompted for user input, you can using the
``force_download`` option:

.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "levine13", force_download=True)

.. note::

    This ``force_download`` currently will not overwrite any existing cache.
    It is used to silence the message only.


Load Specific Samples
-----------------------------

We have showcased the easiest way to load a dataset. As with the case
of ``levine13``, which has only 1 sample, you can load it with the
default settings:

.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "levine13")

In this case, it is all you need to do. If you load ``levine32`` or
``samusik`` with default settings, it will load all samples:

.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "levine32")
    >>> exprs.n_samples
    2

This is helpful if you wish to load all samples. On the other hand,
if you only need one sample, doing so can be quite useful for datasets
like ``samusik`` which has 10 samples. In this case, you can
specify the sample you would like:


.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "samusik", sample=["01", "02"])
    >>> exprs.n_samples
    2

This will load the first two samples. 

Note that not all datasets have the sample indices. Here are the conventions:

=============== ===============================
Dataset           Sample Index
--------------- -------------------------------
Levine 13-dim     ``0``
Levine 32-dim     ``AML08`` and ``AML09``
Samusik           ``01``, ``02``, ..., ``10``
=============== ===============================

This is done for convention of dataests. So, for ``levine32``, you can
load the specific sample:

.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "levine32", sample=["AML08"])
    >>> exprs.n_samples
    1


Built-in Preprocessing
------------------------

Most of the proprocessing for these three datasets have been done already. All you
need to do is to perform Arcsinh transformation. You can either do so yourself,
which is documented in `this tutorial <https://pycytodata.readthedocs.io/en/latest/tutorial/preprocessing.html>`_.
However, we also have a built-in option to perform proprocessing (arcsinh transformation
with a cofactor of 5) upon loading the dataset: 

.. code-block:: python

    >>> exprs = DataLoader.load_dataset(dataset = "levine32", sample=["AML08"], preprocessing = True)

The resulting ``exprs`` object is ready for downstream analyses.

For ``samusik`` and ``levine32`` which have more than 1 sample, there is an
argument to be made that we need cross-batch normalization. Currently. we don't have
a normalization procedure in place (nor do we have packages in the PyCytoData alliance
that fills the role). If you need this preprocessing step done, some other packages may be
needed! 


--------------------------

***************************
Load Your Own Dataset
***************************

For the best flexibity and compatibility with a wide range of workflows, you can use
``PyCytoData`` to load your own CyTOF datasets as well. The dataset must be in plain
text, not the ``fcs`` format.

.. note::

    Although ``PyCytoData`` depends on ``fcsparser``, meaning that there is the capability
    to read ``fcs`` by definition, we currently do not support this usage. This is
    due to the potential ambiguity in the format.


Load Dataset
----------------

To load your dataset, you can use the ``FileIO`` class from our package. By default,
we assume that file is an ``tsv``:

.. code-block:: python

    >>> from PyCytoData import FileIO
    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = True)


If your file is saved with different delimiters, such as commas, you can eaily
specify the ``delim`` parameter:

.. code-block:: python

    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = True, delim=",")

The ``col_names`` parameters indicates whether the first row of data is channel names. This
is oftentimes the convention, but if your data comes in only expression data with channels
stored elsewhere or no channels, you can read in the columns only:

.. code-block:: python

    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = False)


Load Multiple Datasets at Once
------------------------------------

If you have multiple samples in the same format (i.e. The columns are in the same
configuration after accounting for dropped columns as seen in the next section),
you can load multiple samples at once into one single ``PyCytoData`` object: 

.. code-block:: python

    >>> exprs = FileIO.load_expression(["<path_1>", "<path_2>"], col_names = True)
    >>> exprs.n_samples
    2

All samples will be stored in a single object, but sample indices will be preserved.
In the cases that channels are misaligned or mismatched, preprocessing on the users'
side is needed to ensure that they are the same. Otherwise, a `ValueError` will be
thrown. 

Specifying Columns
-------------------

Sometimes it is necessary to load only some of the channels (e.g. You want to load
only the lineage channels and not instrument channels). In this case, you can drop
channels by indices:

.. code-block:: python

    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = True, drop_columns=[0,1])

which will drop the first two columns of the dataset. 

Setting Lineage Channels
--------------------------

By default, the ``FileIO`` doesn't set lineage channels. You can set your lineage channels
manually. If all your channels are lineage channels:

.. code-block:: python

    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = True)
    >>> exprs.lineage_channels = exprs.channels

Or if only specific channels are lineage channels:

.. code-block:: python

    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = True)
    >>> exprs.lineage_channels = exprs.channels[5:]

or by specifying names:

.. code-block:: python

    >>> exprs = FileIO.load_expression("<Your_File_Path>", col_names = True)
    >>> exprs.lineage_channels = ["CD13", "CD4", "CD5"]


-----------------------------------

**********************
General FileIO
**********************

We also support an interface for loading and save lists and arrays. They're mostly
implemented as a convenient wrapper for internal usage, but these may be helpful for
working with CyTOF data. For most usecases, it is using the ``FileIO.load_expression``
will suffice.

Load an Array from Text
------------------------

This behaves similarly as as ``load_expression``, except that it cannot handle column
names and it returns a ``np.ndarray`` instead. To read an array with floats only, you
can use:

.. code-block:: python

    >>> exprs = FileIO.load_delim("<Your_File_Path>")

Alternatvie, if your first row is the column names, you can skip it:

.. code-block:: python

    >>> exprs = FileIO.load_delim("<Your_File_Path>", skiprows=1)

If your array consists of other datatypes, you can change the ``dtype`` argument according,
which is passed into ``np.loadtxt``.


Save 2D List to CSV
-----------------------

The ``FileIO.save_2d_list_to_csv`` does what exactly what it sounds like: it save a nested
2D list to a csv file. To do so, you can simply run

.. code-block:: python

    >>> data = [["A", "B", "C"], [1,2,3]]
    >>> FileIO.save_2d_list_to_csv(data=data, path="<your_path>")

By default, this does not overwrite existing files. If you would like to avoid the error
and overwrite, you can set

.. code-block:: python

    >>> FileIO.save_2d_list_to_csv(data=data, path="<your_path>", overwrite=True)

This method is as feature-rich as other methods since we don't use this internally.

Save Numpy Array
-------------------

The ``FileIO.save_np_array`` saves a numpy array as a plain text file. You can save it
simpy by the following:

.. code-block:: python

    >>> data = np.array([[1,2,3], [3,4,5]])
    >>> FileIO.save_np_array(array=data, path="<your_path>")

You can also specify the column names to be saved as the first row:

.. code-block:: python

    >>> data = np.array([[1,2,3], [3,4,5]])
    >>> col_names = np.array(["CD13", "CD4", "CD3"])
    >>> FileIO.save_np_array(array=data, path="<your_path>", col_names=col_names)
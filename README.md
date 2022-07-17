# PyCytoData
> An elegant data analysis tool for CyTOF.

This package is an all-in-one CyTOF data analysis package for your experiments. From loading datasets to DR and evaluation, you have a consistent interface and readable code every step along the way. There is also support for some of ``HDCytoData``'s benchmark datasets as originally implemented in R by Weber & Soneson (2019) in [this repository](https://github.com/lmweber/HDCytoData). Why wait? Start your PyCytoData journal right here, right now! 

## Installation

You can install ``PyCytoData`` easily from ``pip``:

```
pip install PyCytoData
```

or from ``conda``:

```
conda install pycytodata -c kevin931 -c bioconda
```

If you wish to use ``CytofDR`` along with PyCytoData, use can optionally install it as well:

```
pip install CytofDR
```

For more information on optional dependencies or installation details, look [here](https://pycytodata.readthedocs.io/en/latest/installation.html).

## Install and Load Benchmark Datasets

You can load the data easily with the following python snippet:

```python
>>> from PyCytoData import DataLoader

>>> exprs = DataLoader.load_dataset(dataset = "levine13")
>>> exprs.expression_matrix # Expression matrix
>>> exprs.cell_types # Cell types
>>> exprs.sample_index # Sample index
>>> exprs.features # The feature/marker names
```

The resulting ``exprs`` is a ``PyCytoData`` object, which is easy to use. The expression matrix, cell types (if available), and sample index are directly accessible with attributes, and they are all stored as **numpy.array**. You can also access some metadata of the object with the following attributes:

```python
>>> exprs.n_cells
>>> exprs.n_cell_types
>>> exprs.n_samples
>>> exprs.n_features
```

All these metadata is automatically set, and there is protection in place for unintended changes. You can also add a sample with the following:

```python
>>> exprs.add_sample(expression_matrix, cell_types, sample_index) # All inputs should be ArrayLike
```

**Note**: The data are downloaded from a server instead of being shipped with this package. Each dataset only needs to be downloaded once, which is automatically managed. During the first-time download of the data, a command-line confirmation is needed.

## Bring Your Own Dataset (BYOD)

Yes, you read it right! You can load your own datasets. Currently, we only support reading in plain text files with saved with delimiters. The data need to have cells as rows and features as columns. To do load them in as a ``PyCytoData`` object, you can simply do the following:

```python
>>> from PyCytoData import FileIO

>>> FileIO.load_delim(files="/path", # Path to file
...                   col_names=True, # Whether the first row is feature (column) names 
...                   delim="\t" # Delimiter
...                  ) 
```

If your experiment has multiple samples, you can simply import them together:

```python
>>> from PyCytoData import FileIO

>>> expression_paths = ["path1", "path2", "path3"]
>>> FileIO.load_delim(files=expression_paths, # Path to file
...                   col_names=True, # Whether the first row is feature (column) names 
...                   delim="\t" # Delimiter
...                  ) 
```

In this case, the expression matrices are concatenated automatically without any normalization. To access particular samples, you can access the ``sample_index`` of the attribute and use the standard ``numpy`` indexing techniques.

**Note:** This technique does not automatically load cell types. In fact, it does **not** not mixed-datatype array, except for column names. You will need to read in cell types and set them using the ``cell_types`` attribute of the object. 

## Preprocessing

Currently, ``levine13``, ``levine32``, and ``samusik`` have all been mostly preprocessed. All you need to do is to perform ``aecsinh`` transformaion. You can simply do this:

```python
>>> from PyCytoData import DataLoader

>>> exprs = DataLoader.load_dataset(dataset = "levine13", preprocess=True)
```

When you perform BYOD, you can have much more flexibility:

```python
>>> from PyCytoData import FileIO

>>> byod = FileIO.load_delim(files="/path", # Path to file
...                          col_names=True, # Whether the first row is feature (column) names 
...                          delim="\t" # Delimiter
...                         )
>>> byod.lineage_channels = ["CD4", "CD8", "FoxP3", "CD15"]
>>> byod.preprocess(arcsinh=True,
...                 gate_debris_removal=True,
...                 gate_intact_cells=True,
...                 gate_live_cells=True,
...                 gate_center_offset_residual=True,
...                 bead_normalization=True)

byod.expression_matrix # This is preprocessed
```
As the example shows, we support five unique preprocessing steps! And of course, you can use a subset of these to suit your own needs! By default, we automatically detect the necessary channels, such as "Bead1" or "Center". However, if your dataset is unconventionally named, our auto-detect algorithm may fail. Thus, we can perform a manual override:

```python
>>> byod.preprocess(arcsinh=True,
...                 gate_debris_removal=True,
...                 gate_intact_cells=True,
...                 gate_live_cells=True,
...                 gate_center_offset_residual=True,
...                 bead_normalization=True,
...                 bead_channels = ["1bead", "2bead"],
...                 time_channel = ["clock"])
```

## Dimension Reduction

If you wish to run DR on your dataset, you can easily do so as well if you have ``CytofDR`` installed (assume you have loaded the dataset and preprocessed it accordingly):

```python
>>> exprs.run_dr_methods(methods = ["PCA", "UMAP", "ICA"])
Running PCA
Running ICA
Running UMAP
>>> type(exprs.reductions)
<class 'CytofDR.dr.Reductions'>
```
The ``reductions`` attribute is a ``Reductions`` object from ``CytofDR``. You can perform all downstream DR workflows as usual.

## Datasets Supported

We only support the following datasets as of now. The *Literal* is the string literal used in this package to refer to the datasets whereas the *Dataset Name* is what these datasets are more commonly known for.

| Dataset Name | Literal |
| --- | --- |
| Levine-13dim | levine13 |
| Levine-32dim | levine32 |
| Samusik | samusik |

More datasets will be added in the future to be fully compatible with HDCytoData and to potentially incorporate other databases.

## Documentation

For detailed documentation along with tutorials and API Reference, please visit our [Official Documentation](https://pycytodata.readthedocs.io/en/latest/). This is automatically updated with each update.

If you prefer to build documentation on your own, refer to [this guide](https://pycytodata.readthedocs.io/en/latest/change/build.html) for more details.

## Latest Release: 0.1.0

This is the first official release of ``PyCytoData`` with LTS.

### Bug Fixes

- Fix the default cofactor of ``PyCytoData.preprocess``
- Fix lineage channels for ``levine32`` and ``samusik``

### Changes and New Features

- ``PyCytoData.run_dr_methods`` now runs on lineage channels only
- Add tutorials and references for documentation
- Improve docstrings and general documentations
- Add preprocessing options for loading benchmark datasets
- Relax ``numpy`` dependency to >=1.20
- Add descriptors for PyPI releases

### Deprecations

- Permanently remove ``FileIO.make_dir`` function for safety reasons


## References

If you use ``PyCytoData`` to perform DR, citing the [our DR Review paper](https://doi.org/10.1101/2022.04.26.489549) is highly appreciated:

```
@article {Wang2022.04.26.489549,
	author = {Wang, Kaiwen and Yang, Yuqiu and Wu, Fangjiang and Song, Bing and Wang, Xinlei and Wang, Tao},
	title = {Comparative Analysis of Dimension Reduction Methods for Cytometry by Time-of-Flight Data},
	elocation-id = {2022.04.26.489549},
	year = {2022},
	doi = {10.1101/2022.04.26.489549},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/06/02/2022.04.26.489549},
	eprint = {https://www.biorxiv.org/content/early/2022/06/02/2022.04.26.489549.full.pdf},
	journal = {bioRxiv}
}
```

If you use ``Cytomulate`` with this package, [our paper](https://doi.org/10.1101/2022.06.14.496200) can be cited here:

```
@article {Yang2022.06.14.496200,
	author = {Yang, Yuqiu and Wang, Kaiwen and Lu, Zeyu and Wang, Tao and Wang, Xinlei},
	title = {Cytomulate: Accurate and Efficient Simulation of CyTOF data},
	elocation-id = {2022.06.14.496200},
	year = {2022},
	doi = {10.1101/2022.06.14.496200},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/06/16/2022.06.14.496200},
	eprint = {https://www.biorxiv.org/content/early/2022/06/16/2022.06.14.496200.full.pdf},
	journal = {bioRxiv}
}
```

If you use the builtin datasets, please visit our [Reference Page](https://pycytodata.readthedocs.io/en/latest/references.html) and cite the papers accordingly.

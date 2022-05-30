# PyCytoData
> An elegant data analysis tool for CyTOF.

This package is an all-in-one CyTOF data analysis package for your experiments. From loading datasets to DR and evaluation, you have a consistent interface and readable code every step along the way. There is also support for some of ``HDCytoData``'s benchmark datasets as originally implemented in R by Weber & Soneson (2019) in [this repository](https://github.com/lmweber/HDCytoData). Why wait? Start your PyCytoData journal right here, right now! 

## Installation

We're currently under development! To install, you can do the following:

```shell
git clone https://github.com/kevin931/PyCytoData
cd PyCytoData
python setup.py deveop
```

This approach will allow you to use the package while developing!

### Dependencies

We need the following dependencies:

- fcsparser
- pandas
- numpy

## Install and Load Benchmark Datasets

You can load the data easily with the following python snippet:

```python
from PyCytoData import DataLoader

exprs = DataLoader.load_dataset(dataset = "levine13")
exprs.expression_matrix # Expression matrix
exprs.cell_types # Cell types
exprs.sample_index # Sample index
exprs.features # The feature/marker names
```

The resulting ``exprs`` is a ``PyCytoData`` object, which is easy to use. The expression matrix, cell types (if available), and sample index are directly accessible with attributes, and they are all stored as **numpy.array**. You can also access some metadata of the object with the following attributes:

```python
exprs.n_cells
exprs.n_cell_types
exprs.n_samples
exprs.n_features
```

All these metadata is automatically set, and there is protection in place for unintended changes. You can also add a sample with the following:

```python
exprs.add_sample(expression_matrix, cell_types, sample_index) # All inputs should be ArrayLike
```

**Note**: The data are downloaded from a server instead of being shipped with this package. Each dataset only needs to be downloaded once, which is automatically managed. During the first-time download of the data, a command-line confirmation is needed. To override this, you can do the following: 

```python 
from PyCytoData import DataLoader

exprs = DataLoader.load_dataset(dataset = "levine13", force_download = True)
```

## Bring Your Own Dataset (BYOD)

Yes, you read it right! You can load your own datasets. Currently, we only support reading in plain text files with saved with delimiters. The data need to have cells as rows and features as columns. To do load them in as a ``PyCytoData`` object, you can simply do the following:

```python
from PyCytoData import FileIO

FileIO.load_delim(files="/path", # Path to file
                  col_names=True, # Whether the first row is feature (column) names 
                  delim="\t" # Delimiter
                 ) 
```

If your experiment has multiple samples, you can simply import them together:

```python
from PyCytoData import FileIO

expression_paths = ["path1", "path2", "path3"]
FileIO.load_delim(files=expression_paths, # Path to file
                  col_names=True, # Whether the first row is feature (column) names 
                  delim="\t" # Delimiter
                 ) 
```

In this case, the expression matrices are concatenated automatically without any normalization. To access particular samples, you can access the ``sample_index`` of the attribute and use the standard ``numpy`` indexing techniques.

**Note:** This technique does not automatically load cell types. In fact, it does **not** not mixed-datatype array, except for column names. You will need to read in cell types and set them using the ``cell_types`` attribute of the object. 

## Preprocessing

Currently, ``levine13``, ``levine32``, and ``samusik`` have all been mostly preprocessed. All you need to do is to perform ``aecsinh`` transformaion. You can simply do this:

```python
from PyCytoData import DataLoader

exprs = DataLoader.load_dataset(dataset = "levine13")
exprs.preprocess(arcsinh=True)
```

When you perform BYOD, you can have much more flexibility:

```python
from PyCytoData import FileIO

byod = FileIO.load_delim(files="/path", # Path to file
                         col_names=True, # Whether the first row is feature (column) names 
                         delim="\t" # Delimiter
                        )
byod.lineage_channels = ["CD4", "CD8", "FoxP3", "CD15"]
byod.preprocess(arcsinh=True,
                gate_debris_removal=True,
                gate_intact_cells=True,
                gate_live_cells=True,
                gate_center_offset_residual=True,
                bead_normalization=True)

byod.expression_matrix # This is preprocessed
```
As the example shows, we support five unique preprocessing steps! And of course, you can use a subset of these to suit your own needs! By default, we automatically detect the necessary channels, such as "Bead1" or "Center". However, if your dataset is unconventionally named, our auto-detect algorithm may fail. Thus, we can perform a manual override:

```python
byod.preprocess(arcsinh=True,
                gate_debris_removal=True,
                gate_intact_cells=True,
                gate_live_cells=True,
                gate_center_offset_residual=True,
                bead_normalization=True,
                bead_channels = ["1bead", "2bead"],
                time_channel = ["clock"])
```

## Datasets Supported

We only support the following datasets as of now. The *Literal* is the string literal used in this package to refer to the datasets whereas the *Dataset Name* is what these datasets are more commonly known for.

| Dataset Name | Literal |
| --- | --- |
| Levine-13dim | levine13 |
| Levine-32dim | levine32 |
| Samusik | samusik |

More datasets will be added in the future to be fully compatible with HDCytoData and to potentially incorporate other databases.

## Documentation

We use ``sphinx`` and ``readthedocs`` for documentation! You will need to install the following packages:

- sphinx
- sphinx-rtd-theme
- sphinx-git
- sphinxcontrib-autoprogram
- sphinx-autodoc-typehints

We currently don't have an online documentation. You will need to build the docs on your own! More detailed docs coming soon!

## Unit Testing

You will need the following packages:

- pytest
- pytest-cov
- pytest-mock
- coverage

## References

[Levine J.H., Simonds E.F. Bendall S.C., Davis KL, Amir el-A.D., Tadmor M.D., Litvin O., Fienberg H.G., Jager A., Zunder E.R., Finck R., Gedman A.L., Radtke I., Downing J.R., & Pe'er D., Nolan G.P. "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis." *Cell*. 2015 Jul 2;162(1):184-97. doi: 10.1016/j.cell.2015.05.047.](https://pubmed.ncbi.nlm.nih.gov/26095251/)

[Samusik et al. (2016), "Automated mapping of phenotype space with single-cell data", *Nature Methods, 13*(6), 493-496](https://www.ncbi.nlm.nih.gov/pubmed/27183440)

[Weber L.M. and Soneson C. (2019). "HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats." *F1000Research, 8*:1459, v2.](https://f1000research.com/articles/8-1459)

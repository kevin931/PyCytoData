# PyCytoData
> Interface to HDCytoData for CyTOF Benchmark Datasets

This package downloads and load some popular CyTOF benchmark datasets as published by previous studies. The original concept of such a package is conceived and implemented in R by Weber & Soneson (2019) in [this repository](https://github.com/lmweber/HDCytoData). PyCytoData brings this to Python while also adding our own flavor this incredibly
useful tool.

## Installation

We're currently under development! To install, you can do the following:

```shell
git clone https://github.com/kevin931/PyCytoData
cd PyCytoData
python setup.py deveop
```

This approach will allow you to use the package while developing!

### A note on cytomulate

Currently, this package depends on [cytomulate](https://github.com/kevin931/cytomulate). Please make sure that the development version of cytomulate is installed as well!

### Other Dependencies

We need the following dependencies:

- fcsparser
- pandas
- numpy

## Install and Load Data

You can load the data easily with the following python codes:

```python
from PyCytoData import DataLoader

exprs = DataLoader.load_data(dataset = "levine13")
```

The resulting ``exprs`` is a list of numpy arrays: the first consists of column names and the second is the expression matrix itself.

**Note**: The data are downloaded from a server instead of being shipped with this package. Each dataset only needs to be downloaded once, which is automatically managed. 

During the first-time download of the data, a command-line confirmation is needed. To override this, you can do the following: 

```python 
from PyCytoData import DataLoader

exprs = DataLoader.load_data(dataset = "levine13", force_download = True)
```


## Datasets Supported

We only support the following datasets as of now:

| Dataset | Name |
| --- | --- |
| Levine-13dim | levine13 |

More datasets will be added in the future to be fully compatible with HDCytoData.

## Documentation

We use ``sphinx`` and ``readthedocs`` for documentation! You will need to install the following packages:

- sphinx
- sphinx-rtd-theme
- sphinx-git
- sphinxcontrib-autoprogram
- sphinx-autodoc-typehints


## Unit Testing

You will need the following packages:

- pytest
- pytest-cov
- pytest-mock
- coverage

## References

[Levine J.H., Simonds E.F. Bendall S.C., Davis KL, Amir el-A.D., Tadmor M.D., Litvin O., Fienberg H.G., Jager A., Zunder E.R., Finck R., Gedman A.L., Radtke I., Downing J.R., & Pe'er D., Nolan G.P. "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis." *Cell*. 2015 Jul 2;162(1):184-97. doi: 10.1016/j.cell.2015.05.047.](https://pubmed.ncbi.nlm.nih.gov/26095251/)

[Weber L.M. and Soneson C. (2019). "HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats." *F1000Research, 8*:1459, v2.](https://f1000research.com/articles/8-1459)

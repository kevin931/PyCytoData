from __future__ import annotations

import pandas as pd
import numpy as np
from PyCytoData import exceptions, preprocess

from zipfile import ZipFile
from urllib.request import urlopen
from io import BytesIO

import fcsparser
import _csv
import csv
import os
import pkg_resources
import glob
import re
from copy import deepcopy

from numpy.typing import ArrayLike
from typing import Optional, List, Dict, Literal, Any, Union

OPT_PCK: Dict[str, bool] = {"CytofDR": True}

try:
    from CytofDR import dr
except ImportError:
    OPT_PCK["CytofDR"] = False


def _verbose(message: str, verbose:bool=True):
    if verbose:
        print(message)


class PyCytoData():
    """The CytoData Class for handling CyTOF data.

    This is an all-purpose data class for handling CyTOF data. It is compatible with
    benchmark datasets downloaded from the ``DataLoader`` class as well as users' own
    CyTOF datasets. It has wideranging functionalities, include preprecessing, DR,
    and much more.

    :param expression_matrix: The expression matrix for the CyTOF sample. Rows are cells
        and columns are channels.
    :type expression_matrix: ArrayLike
    :param channels: The name of the channels, defaults to None
    :type channels: ArrayLike
    :param cell_types: The cell types of the cells, defaults to None
    :type cell_types: ArrayLike
    :param sample_index: The indicies or names to indicate samples of each cell.
        This allows the combination of multiple samples into one class, defaults to None
    :type sample_index: ArrayLike
    :param lineage_channels: The names of lineage channels, defaults to None
    :type lineage_channels: ArrayLike
    
    :raises exceptions.ExpressionMatrixDimensionError: The expression matrix is not
        or cannot be cast into a two dimensional array.
    :raises exceptions.DimensionMismatchError: The number of channel names does not agree
        with the number of columns of the expression matrix.
    :raises exceptions.DimensionMismatchError: The number of cell types for all cells does not agree
        with the number of rows of the expression matrix.
    :raises exceptions.DimensionMismatchError: The number of sample indices does not agree
        with the number of rows of the expression matrix.
    """
    
    def __init__(self,
                 expression_matrix: ArrayLike,
                 channels: Optional[ArrayLike]=None,
                 cell_types: Optional[ArrayLike]=None,
                 sample_index: Optional[ArrayLike]=None,
                 lineage_channels: Optional[ArrayLike]=None):

        self._expression_matrix: np.ndarray = np.array(expression_matrix)
        if len(self._expression_matrix.shape) != 2:
            raise exceptions.ExpressionMatrixDimensionError(shape=self._expression_matrix.shape)
        
        self._n_cells = self._expression_matrix.shape[0]
        self._n_channels = self._expression_matrix.shape[1]
        
        if channels is not None:
            self._channels = np.array(channels)
        else:
            self._channels = np.array(["Channel" + str(a) for a in range(self._n_channels)])
        
        if cell_types is None:
            self._cell_types = np.full(self.n_cells, None)
        else:
            self._cell_types = np.array(cell_types)
            
        if sample_index is None:
            self._sample_index = np.repeat(0, self.n_cells)
        else:
            self._sample_index = np.array(sample_index)
        
        self._n_samples: int = len(set(self._sample_index))
        self._n_cell_types: int = len(set(self._cell_types))
        
        if self._channels.shape[0] != self.n_channels:
            raise exceptions.DimensionMismatchError(n=self.n_channels, var = "channels")
        if self._cell_types.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "cell_types")
        if self._sample_index.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "sample_index")
        if np.unique(self._channels).shape[0] != self._n_channels:
            raise ValueError("Channel names not unique: This can result in ambiguities.")
        
        self._lineage_channels: Optional[np.ndarray] = lineage_channels if lineage_channels is None else np.array(lineage_channels).flatten()
        if self._lineage_channels is not None and not np.all(np.isin(self._lineage_channels, self._channels)):
            raise ValueError("Some lineage channels are not listed in channel names.")
        
        self._reductions: Optional[dr.Reductions] = None
    
    
    def add_sample(self, expression_matrix: ArrayLike, sample_index: ArrayLike, cell_types: Optional[ArrayLike]=None):
        """Add another CyTOF sample from the same experiment.

        This method allows users to combine samples into existing samples.
        The data must be in the same shape. Sample indices must be provided
        so that the class can properly index these samples using names.

        :param expression_matrix: The expression matrix of the new sample.
        :type expression_matrix: ArrayLike
        :param sample_index: The sample indicies to name the sample.
        :type sample_index: ArrayLike
        :param cell_types: The cell types of each cell, defaults to None
        :type cell_types: Optional[ArrayLike], optional
        :raises exceptions.ExpressionMatrixDimensionError: The expression matrix cannot be cast
        :raises exceptions.DimensionMismatchError: The number of sample indices
        
        :raises exceptions.DimensionMismatchError: _description_
        """
        expression_matrix = np.array(expression_matrix)
        sample_index = np.array(sample_index)
        
        if len(expression_matrix.shape) != 2:
            raise exceptions.ExpressionMatrixDimensionError(expression_matrix.shape)
        if expression_matrix.shape[1] != self.n_channels:
            raise exceptions.ExpressionMatrixDimensionError(expression_matrix.shape)
        if sample_index.shape[0] != expression_matrix.shape[0]:
            raise exceptions.DimensionMismatchError(n=expression_matrix.shape[0], var = "sample_index")
        if cell_types is not None and np.array(cell_types).shape[0] != expression_matrix.shape[0]:
            raise exceptions.DimensionMismatchError(n=expression_matrix.shape[0], var = "cell_types")
        
        self.expression_matrix = np.concatenate((self.expression_matrix, expression_matrix))
        self.sample_index = np.concatenate((self.sample_index, sample_index))
        
        if cell_types is None:
            self.cell_types = np.concatenate((self.cell_types, np.full(expression_matrix.shape[0], None)))
        else:
            self.cell_types = np.concatenate((self.cell_types, np.array(cell_types)))
         
         
    def preprocess(self,
                   arcsinh: bool=False,
                   gate_debris_removal: bool=False,
                   gate_intact_cells: bool=False,
                   gate_live_cells: bool=False,
                   gate_center_offset_residual: bool = False,
                   bead_normalization: bool=False,
                   auto_channels: bool=True,
                   bead_channels: Optional[ArrayLike]=None,
                   time_channel: Optional[ArrayLike]=None,
                   cor_channels: Optional[ArrayLike]=None,
                   dead_channel: Optional[ArrayLike]=None,
                   DNA_channels: Optional[ArrayLike]=None,
                   cofactor: int=2,
                   cutoff_DNA_sd: float=2,
                   dead_cutoff_quantile: float=0.03,
                   cor_cutoff_quantile: float=0.03,
                   verbose: bool=True):
        """Preprocess the expression matrix.

        This is a one-size-fits-all method to preprocess the CyTOF sample using the ``preprocess``
        module. The preprocessing consists of the following steps:
        
        1. Arcsinh transformation. 
        2. Gate to remove debris.
        3. Gate for intact cells.
        4. Gate for live cells.
        5. Gate for anomalies using center, offset, and residual channels. 

        :param data: The expression matrix array of two dimensions.
        :param gate_debris_removal: Whether to gate to remove debris, defaults to True.
        :type gate_debris_removal: bool
        :param gate_intact_cells: Whether to gate for intact cells, defaults to True.
        :type gate_intact_cells: bool
        :param gate_live_cells: Whether to gate for live cells, defaults to True.
        :type gate_live_cells: bool
        :param gate_center_offset_residual: Whether to gate using center, offset, and residual channels, defaults to True.
        :type gate_center_offset_residual: bool
        :param bead_normalizations: Whether to perform bead normalization, defaults to True.
        :type bead_normalizations: bool
        :param auto_channels: Allow the method to recognize instrument and other non-lineage channels automatically.
            This can be overwritten by specifying channels in ``bead_channels``, ``time_channel``, ``cor_channels``,
            ``dead_channel``, and ``DNA_channels``, defaults to True.
        :type auto_channels: bool
        :param bead_channels: The bead channels as specify by name, defaults to None
        :type bead_channels: ArrayLike, optional
        :param time_channel: The time channel as specify by name, defaults to None
        :type time_channel: ArrayLike, optional
        :param cor_channels: The Center, Offset, and Residual channels as specify by name, defaults to None
        :type cor_channels: ArrayLike, optional
        :param dead_channel: The dead channels as specify by name, defaults to None
        :type dead_channel: ArrayLike, optional
        :param DNA_channels: The DNA channels as specify by name, defaults to None
        :type DNA_channels: ArrayLike, optional
        :param cutoff_DNA_sd: The standard deviation cutoff for DNA channels. Here, we
            specifically measure how many standard deviations away from the mean, defaults to 2
        :type cutoff_DNA_sd: float
        :param dead_cutoff_quantile: The cutoff quantiles for dead channels. The top specified quantile
            will be excluded, defaults to 0.03
        :type dead_cutoff_quantile: float
        :param cor_cutoff_quantile: The cutoff quantiles for Center, Offset, and Residual channels. Both the top
            and bottom specified quantiles will be excluded, defaults to 0.03
        :type cor_cutoff_quantile: float
        :param verbose: Whether to print out progress.
        :type verbose: bool

        :return: The gated expression matrix.
        :rtype: np.ndarray
        """
        
        expression_processed: np.ndarray = deepcopy(self.expression_matrix)
        indices: np.ndarray = np.arange(0, self.n_cells)

        channels = self.channels.tolist()
        if auto_channels:
            auto_channel_error: List[str] = []
            if bead_channels is None and (gate_debris_removal or bead_normalization):
                bead_channels = list(filter(lambda channel: re.match("^bead", channel, re.IGNORECASE), channels)) #type: ignore
                if len(bead_channels) == 0: auto_channel_error.append("bead_channels")
            if DNA_channels is None and gate_intact_cells:
                DNA_channels = list(filter(lambda channel: re.match("dna", channel, re.IGNORECASE), channels)) #type: ignore
                if len(DNA_channels) == 0: auto_channel_error.append("DNA_channels")
            if dead_channel is None and gate_live_cells:
                dead_channel = list(filter(lambda channel: re.match("dead", channel, re.IGNORECASE), channels)) #type: ignore
                if len(dead_channel) == 0: auto_channel_error.append("dead_channel")
            if time_channel is None and bead_normalization:
                time_channel = list(filter(lambda channel: re.match("time", channel, re.IGNORECASE), channels)) #type: ignore
                if len(time_channel) == 0: auto_channel_error.append("time_channel")
            if cor_channels is None and gate_center_offset_residual:
                cor_channels = list(filter(lambda channel: re.match("residual", channel, re.IGNORECASE), channels)) #type: ignore
                cor_channels += list(filter(lambda channel: re.match("center", channel, re.IGNORECASE), channels)) #type: ignore
                cor_channels += list(filter(lambda channel: re.match("offset", channel, re.IGNORECASE), channels)) #type: ignore
                if len(cor_channels) < 3: auto_channel_error.append("cor_channels")
                
            if len(auto_channel_error) > 0:
                raise exceptions.AutoChannelError(auto_channel_error)
            
        indices_temp: np.ndarray
        
        if arcsinh:
            _verbose("Runinng Arcsinh transformation...", verbose=verbose)
            expression_processed = preprocess.arcsinh(expression_processed, self.channels, transform_channels=self.lineage_channels, cofactor=cofactor)
        
        if gate_debris_removal:
            _verbose("Runinng debris remvoal...", verbose=verbose)
            assert bead_channels is not None
            expression_processed, indices_temp = preprocess.gate_debris_removal(expression_processed, self.channels, bead_channels)
            indices = indices[indices_temp]
        
        if gate_intact_cells:
            _verbose("Runinng gating intact cells...", verbose=verbose)
            assert DNA_channels is not None
            expression_processed, indices_temp = preprocess.gate_intact_cells(expression_processed, self.channels, DNA_channels, cutoff_DNA_sd)
            indices = indices[indices_temp]
            
        if gate_live_cells:
            _verbose("Runinng gating live cells...", verbose=verbose)
            assert dead_channel is not None
            expression_processed, indices_temp = preprocess.gate_live_cells(expression_processed, self.channels, dead_channel, dead_cutoff_quantile)
            indices = indices[indices_temp]
            
        if gate_center_offset_residual:
            _verbose("Runinng gating Center, Offset, and Residual...", verbose=verbose)
            assert cor_channels is not None
            expression_processed, indices_temp = preprocess.gate_center_offset_residual(expression_processed, self.channels, cor_channels, cor_cutoff_quantile)
            indices = indices[indices_temp]
            
        if bead_normalization:
            _verbose("Runinng bead normalization...", verbose=verbose)
            assert bead_channels is not None
            assert time_channel is not None
            assert self.lineage_channels is not None
            expression_processed, indices_temp = preprocess.bead_normalization(expression_processed, self.channels, bead_channels, time_channel, self.lineage_channels)
            indices = indices[indices_temp]
        
        self.expression_matrix = expression_processed
        if gate_debris_removal or gate_intact_cells or gate_live_cells or gate_center_offset_residual or bead_normalization:
            self.cell_types = self.cell_types[indices]
            self.sample_index = self.sample_index[indices]
            
            
    def run_dr_methods(self,
                       methods: Union[str, List[str]]="all",
                       out_dims: int=2,
                       n_jobs: int=-1,
                       verbose: bool=True,
                       suppress_error_msg: bool=False
                   ):
        """Run dimension reduction methods.

        This is a one-size-fits-all dispatcher that runs all supported methods in the module. It
        supports running multiple methods at the same time at the sacrifice of some more
        granular control of parameters. If you would like more customization, please use the
        ``CytofDR`` package directly.
        
        :param methods: DR methods to run (not case sensitive).
        :type methods: Union[str, List[str]]
        :param out_dims: Output dimension of DR.
        :type out_dims: int
        :param n_jobs: The number of jobs to run when applicable, defaults to -1.
        :type n_jobs: int
        :param verbose: Whether to print out progress, defaults to ``True``.
        :type verbose: bool
        :param suppress_error_msg: Whether to suppress error messages print outs, defaults to ``False``.
        :type supress_error_msg: bool
        
        :raises ImoportError: ``CytofDR`` is not installed.
        """
        if not OPT_PCK["CytofDR"]:
            raise ImportError("`CytofDR` is not installed. Please install `CytofDR` first.")
        
        self.reductions = dr.run_dr_methods(data=self.expression_matrix, methods=methods, out_dims=out_dims,
                                            n_jobs=n_jobs, verbose=verbose, suppress_error_msg=suppress_error_msg)
        
        self.reductions.add_evaluation_metadata(original_data=self.expression_matrix)
        if np.any(self.cell_types != None):
            self.reductions.add_evaluation_metadata(original_cell_types=self.cell_types)
            
    
    def subset(self, sample: Optional[ArrayLike]=None, cell_types: Optional[ArrayLike]=None, not_in: bool=False, in_place: bool=True) -> Optional[PyCytoData]:
        """Subset the dataset with specific cell types or samples.

        This method allows you to subset and retain only certain samples or cell types of interest.

        :param sample: The names of the samples to perform subset, defaults to None
        :type sample: Optional[ArrayLike], optional
        :param cell_types: The name of the cell types to perform subset, defaults to None
        :type cell_types: Optional[ArrayLike], optional
        :param not_in: Whether to filter out the provided cell types or samples, defaults to False
        :type not_in: bool, optional
        :param in_place: Whether to perform the subset in place. If not, a new object will be created and returned. defaults to True.
        :type in_place: bool, optional
        :return: A new PyCytoData after subsetting
        :rtype: PyCytoData, optional
        
        :raises ValueErro: Filtering out all cells with nothing in the expression matrix, which is unsupported.
        """
        if sample is None and cell_types is None:
            raise TypeError("'sample' and 'cell_types' cannot both be None.")
        
        filter_condition: np.ndarray = np.repeat(True, self.n_cells)
        if sample is not None:
            if not isinstance(sample, np.ndarray):
                sample = np.array(sample)
            filter_condition = np.logical_and(filter_condition, np.isin(self.sample_index, sample)) 
            
        if cell_types is not None:
            if not isinstance(cell_types, np.ndarray):
                cell_types = np.array(cell_types)
            filter_condition = np.logical_and(filter_condition, np.isin(self.cell_types, cell_types))
            
        if not_in:
             filter_condition = np.invert(filter_condition)
             
        if not np.any(filter_condition):
            raise ValueError("Filtering out all cells with nothing in the expression matrix. This is unsupported.")
        
        if not in_place:
            new_exprs: PyCytoData = deepcopy(self)
            new_exprs.expression_matrix = new_exprs.expression_matrix[filter_condition, :]
            new_exprs.sample_index = new_exprs.sample_index[filter_condition]
            new_exprs.cell_types = new_exprs.cell_types[filter_condition]
            return new_exprs
            
        self.expression_matrix = self.expression_matrix[filter_condition, :]
        self.sample_index = self.sample_index[filter_condition]
        self.cell_types = self.cell_types[filter_condition]
         
         
    @property
    def expression_matrix(self) -> np.ndarray:
        """Getter for the expression matrix.
        
        :return: The expression matrix.
        :rtype: np.ndarray
        """
        return self._expression_matrix   
    
    
    @expression_matrix.setter
    def expression_matrix(self, expression_matrix: ArrayLike):
        """Set expression matrix.

        :param expression_matrix: The new expression matrix.
        :type expression_matrix: ArrayLike
        :raises exceptions.ExpressionMatrixDimensionError: The expression matrix is not two-dimensional.
        """
        expression_matrix = np.array(expression_matrix)
        if len(expression_matrix.shape) != 2:
            raise exceptions.ExpressionMatrixDimensionError(expression_matrix.shape)
        self.n_cells = expression_matrix.shape[0]
        self.n_channels = expression_matrix.shape[1]
        self._expression_matrix = expression_matrix
        
  
    @property
    def sample_index(self) -> np.ndarray:
        """Getter for sample_index.
        
        :return: The sample index.
        :rtype: np.ndarray
        """
        return self._sample_index
    
    
    @sample_index.setter
    def sample_index(self, sample_index: ArrayLike):
        """Set sample_index.

        :param sample_index: The sample index for each cell.
        :type sample_index: ArrayLike
        :raises exceptions.DimensionMismatchError: Sampel indices' length does not agree with number of features.
        """
        sample_index = np.array(sample_index)
        if sample_index.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "sample_index")
        self._sample_index = sample_index
        self.n_samples = len(set(self._sample_index))
        
        
    @property
    def cell_types(self) -> np.ndarray:
        """Getter for sample_index.
        
        :return: The cell types.
        :rtype: np.ndarray
        """
        return self._cell_types
    
    
    @cell_types.setter
    def cell_types(self, cell_types: ArrayLike):
        """Set cell_types.

        :param cell_types: The cell types.
        :type cell_types: ArrayLike
        :raises exceptions.DimensionMismatchError: Cell types' length does not agree with number of cells.
        """
        cell_types = np.array(cell_types)
        if cell_types.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "cell_types")
        self._cell_types = cell_types
        self.n_cell_types = len(set(self.cell_types))
        
        
    @property
    def channels(self) -> np.ndarray:
        """Getter for sample_index.
        
        :return: The sample index.
        :rtype: np.ndarray
        """
        return self._channels
    
    
    @channels.setter
    def channels(self, channels: ArrayLike):
        """Set channels.

        :param channels: The channel names.
        :type channels: ArrayLike
        :raises exceptions.DimensionMismatchError: Channels names' length does not agree with number of features.
        """
        channels = np.array(channels)
        if channels.shape[0] != self.n_channels:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "channels")
        self._channels = channels
        
        
    @property
    def n_cells(self) -> int:
        """Getter for n_cells.

        :return: The number of cells.
        :rtype: int
        """
        return self._n_cells
    
    
    @n_cells.setter
    def n_cells(self, n_cells: int):
        """Set n_cells.

        :param n_cells: The total number of cells in the ``expression_matrix``.
        :type n_cells: int
        :raises TypeError: The input is not an ``int``.
        """
        if not isinstance(n_cells, int):
            raise TypeError(f"'n_cells' has to be 'int' instead of {type(n_cells)}")
        self._n_cells = n_cells
        

    @property
    def n_channels(self) -> int:
        """Getter for n_channels.

        :return: The number of channels.
        :rtype: int
        """
        return self._n_channels
    
    
    @n_channels.setter
    def n_channels(self, n_channels: int):
        """Set n_channels.

        :param n_channels: The total number of channels in the ``expression_matrix``.
        :type n_channels: int
        :raises TypeError: The input is not an ``int``.
        """
        if not isinstance(n_channels, int):
            raise TypeError(f"'n_channels' has to be 'int' instead of {type(n_channels)}")
        self._n_channels = n_channels
        
        
    @property
    def n_samples(self) -> int:
        """Getter for n_samples.

        :return: The number of samples.
        :rtype: int
        """
        return self._n_samples
    
    
    @n_samples.setter
    def n_samples(self, n_samples: int):
        """Set n_samples.

        :param n_samples: The total number of samples in the ``expression_matrix``.
        :type n_samples: int
        :raises TypeError: The input is not an ``int``.
        """
        if not isinstance(n_samples, int):
            raise TypeError(f"'n_samples' has to be 'int' instead of {type(n_samples)}")
        self._n_samples = n_samples
        
        
    @property
    def n_cell_types(self) -> int:
        """"Getter for n_cell_types.

        :return: The number of cell types.
        :rtype: int
        """
        return self._n_cell_types
    
    
    @n_cell_types.setter
    def n_cell_types(self, n_cell_types: int):
        """Set n_cell_types.

        :param n_cell_types: The total number of cell types in the ``expression_matrix``.
        :type n_cell_types: int
        :raises TypeError: The input is not an ``int``.
        """
        if not isinstance(n_cell_types, int):
            raise TypeError(f"'n_cell_types' has to be 'int' instead of {type(n_cell_types)}")
        self._n_cell_types = n_cell_types
        
    @property
    def lineage_channels(self) -> Optional[np.ndarray]:
        """Getter for lineage_channels.

        :return: An array of lineage channels or None.
        :rtype: np.ndarray, optional
        """
        return self._lineage_channels
    
    
    @lineage_channels.setter
    def lineage_channels(self, lineage_channels: ArrayLike):
        """Set lineage_channels.

        :param lineage_channels: The names of the lineage channels in the ``channels``.
        :type lineage_channels: int
        :raises ValueError: Some lineage channels are not listed in channel names.
        """
        if not np.all(np.isin(lineage_channels, self._channels)):
            raise ValueError("Some lineage channels are not listed in channel names.")
        self._lineage_channels: Optional[np.ndarray] = lineage_channels if lineage_channels is None else np.array(lineage_channels).flatten()
        
        
    @property
    def reductions(self) -> Optional[dr.Reductions]:
        return self._reductions
    
    
    @reductions.setter
    def reductions(self, reductions: Optional[dr.Reductions]):
        if not isinstance(reductions, dr.Reductions) and reductions is not None:
            raise TypeError("'reductions' has to of type 'CytofDR.dr.Reductions' or None")
        self._reductions = reductions


class DataLoader():

    # Package data directory and Path
    _data_dir = pkg_resources.resource_filename("PyCytoData", "data/")
    _data_path: Dict[str, str] = {"levine13": _data_dir + "levine13/",
                                  "levine32": _data_dir + "levine32/",
                                  "samusik": _data_dir + "samusik/"}
    # Check data status
    _data_status: Dict[str, bool] = {}
    for d in _data_path.keys():
        _data_status[d] = os.path.exists(_data_path[d])

    @classmethod    
    def load_dataset(cls, dataset: Literal["levine13", "levine32", "samusik"], force_download: bool = False) -> PyCytoData:
        """Load benchmark datasets.

        This methods downloads and load benchmark datasets. The dataset is downloaded only once, which is then
        cached for future use. Currently, we support three datasets:
        
        - ``levine13``
        - ``levine32``
        - ``samusik``

        :param dataset: The name of the dataset.
        :type dataset: Literal['levine13', 'levine32', 'samusik']
        :param force_download: Whether to download dataset regardless of previous cache, defaults to False
        :type force_download: bool
        :return: The loaded dataset.
        :rtype: PyCytoData
        """
        
        if not cls._data_status[dataset]:
            cls._download_data(dataset = dataset, force_download = force_download)
            
        data: PyCytoData = FileIO.load_delim(cls._data_path[dataset]+dataset+".txt", col_names = True)
        
        cell_type_path: str = cls._data_path[dataset] + dataset + "_metadata.txt"
        if os.path.exists(cell_type_path):
            metadata: np.ndarray = np.loadtxt(fname=cell_type_path, dtype ="str", delimiter="\t")
            data.cell_types = metadata[:,0].flatten()
            data.sample_index = metadata[:,1].flatten()
            
        return data
            

    @classmethod
    def _download_data(cls,
                      dataset: Literal["levine13", "levine32", "samusik"],
                      force_download: bool=False) -> int:
        
        """Method to download datasets."""
        urls: Dict[str, List[str]] = {"levine13": ["http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Levine_13dim/Levine_13dim_fcs_files.zip"],
                                       "levine32": ["http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Levine_32dim/Levine_32dim_fcs_files.zip"],
                                       "samusik": ["http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Samusik/Samusik_fcs_files.zip",
                                                   "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Samusik/Samusik_population_IDs.zip"]}

        if not force_download:
            value = input(f"Would you like to download {dataset}? [y/n]")
            
            if value.lower() != "y":
                message_1 = f"\nYou have declined to download {dataset}.\n"
                print(message_1)

                return 1

        # Download message
        message_2 = "\nDownload in progress...\n"
        message_2 += "This may take quite a while, "
        message_2 += "go grab a coffee or cytomulate it!\n"
        print(message_2)
        
        for url in urls[dataset]:
            contents = urlopen(url)
            contents = contents.read()
            zip_file = ZipFile(BytesIO(contents))
            print("Here")
            zip_file.extractall(cls._data_path[dataset])
        
        data: PyCytoData = cls._preprocess(dataset)
        path: str = cls._data_path[dataset] + dataset + ".txt"
        cell_types: np.ndarray = data.cell_types.reshape(data.n_cells, 1)
        sample_index: np.ndarray = data.sample_index.reshape(data.n_cells, 1)
        meta_data: np.ndarray = np.concatenate((cell_types, sample_index), axis=1)
        metadata_path: str = cls._data_path[dataset] + dataset + "_metadata.txt"
        FileIO.save_np_array(data.expression_matrix, path, col_names = data.channels)
        FileIO.save_np_array(meta_data, metadata_path, dtype="%s")
        return 0
    
    
    @classmethod
    def _preprocess(cls, dataset: Literal["levine13", "levine32", "samusik"]) -> PyCytoData:
        
        fcss: List[str] = glob.glob(cls._data_path[dataset] + "*.fcs")
        exprs: pd.DataFrame = pd.DataFrame()
        meta: Dict[str, Any] = {}
        sample_length = []
        temp: pd.DataFrame
        fcs: str
        for fcs in fcss:
            meta, temp = fcsparser.parse(fcs, reformat_meta=True)
            sample_length.append(len(temp))
            exprs = pd.concat([exprs, temp])
        
        if dataset == "levine13":
            data: PyCytoData = cls._preprocess_levine13(fcss, exprs, meta, sample_length)
        elif dataset == "levine32":
            data: PyCytoData = cls._preprocess_levine32(fcss, exprs, sample_length)
        elif dataset == "samusik":
            data: PyCytoData = cls._preprocess_samusik(fcss, exprs, sample_length)
        
        return data
    
    
    @classmethod
    def _preprocess_levine13(cls, fcss: List[str], exprs: pd.DataFrame, meta: Dict[str, Any], sample_length: List[int]) -> PyCytoData:
        expression_matrix: np.ndarray = exprs.to_numpy() 
        # Cell Types
        fcs_cell_types: List[str] = [types.split("_")[-2] for types in fcss]
        cell_types: np.ndarray = np.array([])
        for s in range(len(sample_length)):
            cell_types = np.concatenate((cell_types, np.repeat(fcs_cell_types[s], sample_length[s])))
        # channels
        colnames: np.ndarray = meta['_channels_']["$PnN"].to_numpy()
        # Construct data
        data: PyCytoData = PyCytoData(expression_matrix, colnames, cell_types)
        
        return data
    
    
    @classmethod
    def _preprocess_levine32(cls, fcss: List[str], exprs: pd.DataFrame, sample_length: List[int]) -> PyCytoData:
        colnames: np.ndarray = np.array(exprs.columns)
        expression_matrix: np.ndarray = exprs.to_numpy()
        # Cell Types and Sample Index
        fcs_cell_types: List[str] = [types.split("_")[-2] for types in fcss]
        cell_types: np.ndarray = np.array([])
        sample_index: np.ndarray = np.array([])
        sample_names: List[str] = [types.split("-")[3] for types in fcss]
        for s in range(len(sample_length)):
            cell_types = np.concatenate((cell_types, np.repeat(fcs_cell_types[s], sample_length[s])))
            sample_index = np.concatenate((sample_index, np.repeat(sample_names[s], sample_length[s])))
        # Construct data
        data: PyCytoData = PyCytoData(expression_matrix, colnames, cell_types, sample_index)
        
        return data
    
    
    @classmethod
    def _preprocess_samusik(cls, fcss: List[str], exprs: pd.DataFrame, sample_length: List[int]) -> PyCytoData:
        colnames: np.ndarray = np.array(exprs.columns)
        expression_matrix: np.ndarray = exprs.to_numpy()
        # Sample Index
        sample_index: np.ndarray = np.array([])
        sample_names: List[str] = [types.split("_")[3] for types in fcss]
        for s in range(len(sample_length)):
            sample_index = np.concatenate((sample_index, np.repeat(sample_names[s], sample_length[s])))
            
        # Cell Types
        meta_data: pd.DataFrame = pd.read_csv(cls._data_path["samusik"] + "population_assignments.txt", delimiter = "\t", header = None) #type: ignore
        temp: pd.DataFrame = meta_data[0].str.split("_| ", expand = True) #type: ignore
        meta_data = pd.concat([meta_data, temp], axis = 1)
        meta_array: np.ndarray = meta_data.iloc[:,[1,5,8]].to_numpy()
        
        cell_types: List[str] = ["unassigned"]*expression_matrix.shape[0]
        unique_sample: np.ndarray = np.unique(meta_array[:,1])
        counter: int = 0
        for i, s in enumerate(sample_length):
            subset: np.ndarray = meta_array[meta_array[:,1] == unique_sample[i]]
            for j in range(len(subset)):
                try:
                    cell_types[counter + int(subset[j,2])] = subset[j, 0]
                except TypeError: # paragma: no cover
                    continue
            counter += s
        
        # Construct data
        data: PyCytoData = PyCytoData(expression_matrix, colnames, cell_types, sample_index)
        
        return data


class FileIO():
    
    @staticmethod
    def load_delim(files: Union[List[str], str],
                   col_names: bool=True,
                   drop_columns: Optional[Union[int, List[int]]]=None,
                   delim: str="\t",
                   dtype = float
                   ) -> PyCytoData:
        
        """Load a deliminited text file as a PyCytoData object.

        This method loads a deliminited file and returns a PyCytoData object. The file
        has to be a standard text file containing the expression matrix. Rows are cells
        and columns are channel names. If ``col_names`` is ``True``, the first row of
        the file will be treated as channel names. If multiple file paths are present,
        they will be automatically concatenated into one object, but the sample
        indices will be recorded.

        :raises TypeError: The ``files`` is neither a string nor a list of strings.
        :return: A PyCytoData object.
        :rtype: PyCytoData
        """
        
        if not isinstance(files, str) and not isinstance(files, list):
            raise TypeError("'files' has to be str or a list of str as paths.")
        elif not isinstance(files, list):
            files = [files]
        
        skiprows: int = 1 if col_names else 0
        colnames: Optional[np.ndarray] = None
        data: Optional[PyCytoData] = None
        
        for i, file in enumerate(files):                            
            # Load Data
            f: "np.ndarray" = np.loadtxt(fname=file, dtype=dtype, skiprows=skiprows, delimiter=delim)
            sample_index: np.ndarray = np.full(f.shape[0], i)
            if drop_columns is not None:
                f = np.delete(f, drop_columns, axis=1)
            if i==0:
                data = PyCytoData(expression_matrix=f, sample_index=sample_index)
            else:
                assert data is not None
                data.add_sample(expression_matrix=f, sample_index=sample_index)
                
            # Load column names
            if i==0:
                if col_names:
                    colnames = np.loadtxt(fname=file, dtype ="str", max_rows=1, delimiter=delim)
                    if drop_columns is not None:
                        colnames = np.delete(colnames, drop_columns)
                    skiprows = 1
                else:
                    colnames = np.full(data.n_channels, None)
        
        assert data is not None
        assert colnames is not None
        data.channels = colnames
        
        return data
    
    
    @staticmethod
    def save_2d_list_to_csv(data: List[List[Any]], path: str, overwrite: bool = False):
        """Save a nested list to a CSV file.

        :param data: The nested list to be written to disk
        :type data: List[List[Any]]
        :param path: Path to save the CSV file
        :type path: str
        
        ..note:: 
        
            By default, this method does not overwrite existing files. In case a file exists,
            a ``FileExistsError`` is thrown.
        """
        if os.path.exists(path) and not overwrite:
            raise FileExistsError()
        
        i: int
        j: int
        
        with open(path, "w") as f:      
            w: "_csv._writer" = csv.writer(f)
            for i in range(len(data[0])):
                row: List[Any] = []
                for j in range(len(data)):
                    row.append(data[j][i])
                w.writerow(row)
            
            
    @staticmethod
    def save_np_array(array: "np.ndarray",
                      path: str,
                      col_names: Optional["np.ndarray"]=None,
                      dtype: str="%.18e",
                      overwrite: bool = False) -> None:
        """Save a NumPy array to a plain text file

        :param array: The NumPy array to be saved
        :type array: np.ndarray
        :param file: Path to save the plain text file
        :type file: str
        :param col_names: Column names to be save as the first row, defaults to None
        :type col_names: np.ndarray, optional
        :param dtype: NumPy data type, defaults to "%.18e"
        :type dtype: str, optional
        
        .. note:: 
        
            By default, this method does not overwrite existing files. In case a file exists,
            a ``FileExistsError`` is thrown.
        """
        if os.path.exists(path) and not overwrite:
            raise FileExistsError()
            
        with open(path, "w") as f:
            if col_names is not None:
                f.write("\t".join(list(map(str, col_names))))
                f.write("\n")
            np.savetxt(f, array, delimiter="\t", fmt=dtype)
            
    
    @staticmethod
    def make_dir(dir_path: str, add_number_if_dir_exists: bool = False, _counter: int=0) -> str:
        """Create a new directory

        :param dir_path: Path to the new directory to be created
        :type dir_path: str
        :param add_number_if_dir_exists: If the directory already exists, append a number to the
            name until creation of directory is successful, defaults to True
        :type add_number_if_dir_exists: bool, optional
        :return: The path to the new directory
        :rtype: str
        
        .. Warning:: 
        
            The ``add_number_if_dir_exists`` can be dangerous when this method is run inside
            of a loop. This behavior may be deprecated and removed in the future.
        """
        dir_path = dir_path.rstrip("/")
        if _counter==0:
            new_dir_path = dir_path
        else:
            new_dir_path = dir_path + str(_counter)
        
        try:
            os.makedirs(new_dir_path)
        except FileExistsError:
            if add_number_if_dir_exists: 
                new_dir_path = FileIO.make_dir(dir_path, _counter = _counter+1)
            else:
                raise
            
        return new_dir_path
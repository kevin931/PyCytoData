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
from typing import Optional, List, Dict, Any, Tuple, Union

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
        
    :Additional Attributes:
    
    - **reductions**: A ``reductions`` object for dimension reduction using ``CytofDR``.
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
        
        self._lineage_channels_indices: np.ndarray
        if lineage_channels is None:
            self._lineage_channels_indices = np.arange(self.n_channels)
        else:
            self._lineage_channels_indices = np.where(np.isin(lineage_channels, self.channels))
        
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
                   cofactor: int=5,
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
        6. Bead normalization.

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
        :param cofactor: The cofactor for arcsinh transforatrion, default to 5.
        :type cofactor: int, optional
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
            _verbose("Running bead normalization...", verbose=verbose)
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
        
        self.reductions = dr.run_dr_methods(data=self.expression_matrix[:,self._lineage_channels_indices], methods=methods, out_dims=out_dims,
                                            n_jobs=n_jobs, verbose=verbose, suppress_error_msg=suppress_error_msg)
        
        if np.any(self.cell_types != None):
            self.reductions.add_evaluation_metadata(original_cell_types=self.cell_types)
            
    
    def subset(self, channels: Optional[ArrayLike]=None, sample: Optional[ArrayLike]=None, cell_types: Optional[ArrayLike]=None, not_in: bool=False, in_place: bool=True) -> Optional[PyCytoData]:
        """Subset the dataset with specific cell types or samples.

        This method allows you to subset the data using channels, samples,
        or cell types. In terms of the expression matrix, channels subsets
        are operations on columns, whereas sample or cell type subsets
        are operations on rows.
        
        .. tip::
            
            To index specific channels and get the expression matrix instead of a ``PyCtyoData``
            object, use the ``get_channel_expressions`` method.
            
        .. tip::
        
            To subset by indices, use the ``[]`` syntax, which supports indexing similar to that
            of ``numpy``.

        :param channels: The names of the channels to perform subset, defaults to None.
        :type channels: Optional[ArrayLike], optional
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
        
        :raises ValueError: Filtering out all cells with nothing in the expression matrix, which is unsupported.
        """
        if sample is None and cell_types is None and channels is None:
            raise TypeError("'channels', 'sample', and 'cell_types' cannot all be None.")
        
        channel_filter_condition: np.ndarray = np.repeat(True, self.n_channels)
        if channels is not None:
            if not isinstance(channels, np.ndarray):
                channels = np.array(channels)
            channel_filter_condition = np.logical_and(channel_filter_condition, np.isin(self.channels, channels))
        
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
             channel_filter_condition = np.invert(channel_filter_condition)
             
        if not np.any(filter_condition):
            raise ValueError("Filtering out all cells with nothing in the expression matrix. This is unsupported.")
        
        if not in_place:
            new_exprs: PyCytoData = deepcopy(self)
            new_exprs.expression_matrix = new_exprs.expression_matrix[filter_condition]
            new_exprs.expression_matrix = new_exprs.expression_matrix[:, channel_filter_condition]
            new_exprs.channels = new_exprs.channels[channel_filter_condition]
            if new_exprs.lineage_channels is not None:
                new_exprs.lineage_channels = new_exprs.lineage_channels[np.isin(new_exprs.lineage_channels, new_exprs.channels)]
            new_exprs.sample_index = new_exprs.sample_index[filter_condition]
            new_exprs.cell_types = new_exprs.cell_types[filter_condition]
            return new_exprs
            
        self.expression_matrix = self.expression_matrix[filter_condition]
        self.expression_matrix = self.expression_matrix[:, channel_filter_condition]
        self.channels = self.channels[channel_filter_condition]
        if self.lineage_channels is not None:
            self.lineage_channels = self.lineage_channels[np.isin(self.lineage_channels, self.channels)]
        else:
            self._lineage_channels_indices = np.arange(self.n_channels)
        self.sample_index = self.sample_index[filter_condition]
        self.cell_types = self.cell_types[filter_condition]
        
        
    def get_channel_expressions(self, channels: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
        """Get the expressions of specific channels.

        This method subsets the expression matrix with the specific channels
        specified and returns the expression matrix along with the channel
        names. As opposed to ``subset``, this method is more useful for
        investigating the expressions themselves rather than subsetting the
        object as a whole.

        :param channels: The channel names to subset the data.
        :type channels: Union[str, List[str]]
        :raises TypeError: The channels n
        :raises ValueError: The channels specified are not listed in the channel names.
        :return: A tuple of the expressions and the corresponding channel names.
        :rtype: Tuple[np.ndarray, np.ndarray]
        """
        if not isinstance(channels, np.ndarray):
            channels = np.array(channels)
        if not np.all(np.isin(channels, self.channels)):
            raise ValueError("Some channels are not listed in channel names.")
        channel_indices: np.ndarray = np.isin(self.channels, channels)
        return self.expression_matrix[:, channel_indices], self.channels[channel_indices]
        
        
    def __len__(self) -> int:
        """The length of the PyCytoData Class.

        This method implements the ``len`` of the builtin
        python method. It returns the number of total cells
        in the expression matrix.

        :return: The length of the object.
        :rtype: int
        """
        return self.n_cells
    
    
    def __iadd__(self, new_object: PyCytoData) -> PyCytoData:
        """Concatenate a new `PyCytoData` object with the `+=` operator.

        This essentially works the same way the ``add_sample`` method. However,
        instead of the necessity of providing the expression matrices, sample
        indices, and the cell types manually, the concatenation is automatically
        performed from a new `PyCytoData` object.

        :param new_object: The second `PyCytoData` object.
        :type new_object: PyCytoData
        :raises TypeError: The provided object is not a `PyCytoData` object.
        :return: A new `PyCytoData` object after concatenation.
        :rtype: PyCytoData
        """
        if not isinstance(new_object, PyCytoData):
            raise TypeError("The right hand side has to be a 'PyCytoData' object.")
        self.add_sample(new_object.expression_matrix, sample_index=new_object.sample_index, cell_types=new_object.cell_types)
        return self
    
    
    def __add__(self, new_object: PyCytoData) -> PyCytoData:
        """Concatenate two `PyCytoData` objects with the `+` operator.

        This method concatenates two `PyCytOData` objects together by using
        the `add_sample` method internally. A new `PyCytoData` object is returned.

        :param new_object: The second `PyCytoData` object.
        :type new_object: PyCytoData
        :raises TypeError: The provided object is not a `PyCytoData` object.
        :return: A new `PyCytoData` object after concatenation.
        :rtype: PyCytoData
        """
        if not isinstance(new_object, PyCytoData):
            raise TypeError("The right hand side has to be a 'PyCytoData' object.")
        out_object = deepcopy(self)
        out_object.add_sample(new_object.expression_matrix, sample_index=new_object.sample_index, cell_types=new_object.cell_types)
        return out_object
    
    
    def __str__(self) -> str:
        """String representation of the PyCytoData class.

        This method returns a string containing the most basic metadata
        of the class along with the memory address.

        :return: The string representation of the class.
        :rtype: str
        """
        out_str: str = f"A 'PyCytoData' object with {self.n_cells} cells, {self.n_channels} channels, {self.n_cell_types} cell types, and {self.n_samples} samples at {hex(id(self))}."
        return out_str
    
    
    def __getitem__(self, items: Union[slice, List[int], np.ndarray, Tuple[Union[slice, List[int], np.ndarray], Union[slice, List[int], np.ndarray]]]) -> PyCytoData:
        """The method to index elements of the PyCytoData object.

        This method implements the bracket notation to index part of the class.
        The notation is mostly consistent with the numpy indexing notation with a
        few excetions, which is listed below. When indexing specific cells, the
        metadata are appropriately indexed as well.
        
        A few deviations from the numpy notations:
        
        1. Integer indices are currently not supported. This is because indexing by
           integer returns a 1-d array instead of a 2-d array, which can possibly cause confusion.
        2. Indexing by two lists or arrays with different lengths are supported.
           They are treated to index rows and columns, such as ``exprs[[0,1,2], [3,4]]`` is
           perfectly valid to index the first 3 cells with the fourth and fifth channel.
        
        .. tip::
        
            To index columns/channels by name, use the ``subset`` method instead.

        :param items: The indices for items.
        :type items: Union[int, slice, List[int], Tuple[Any, Any]]
        :raises IndexError: Two or more indices present.
        :raises TypeError: Indexing by integer in either or both axes.
        :raises IndexError: An higher dimensional array is used.
        :raises TypeError: Invalid indices type used.
        :return: An appropriately indexed ``PyCytoData`` object.
        :rtype: PyCytoData
        """
        if isinstance(items, tuple):
            if len(items) > 2:
                raise IndexError("Invalid indices: Must be 1 or 2 indices.")
            if (not isinstance(items[0], slice) and not isinstance(items[0], tuple) and
                not isinstance(items[0], list) and not isinstance(items[0], np.ndarray)):
                raise TypeError("Invalid indices: Must be slice, tuple, list, or numpy array.")
            if (not isinstance(items[1], slice) and not isinstance(items[1], tuple) and
                not isinstance(items[1], list) and not isinstance(items[1], np.ndarray)):
                raise TypeError("Invalid indices: Must be slice, tuple, list, or numpy array.")
            
        if isinstance(items, np.ndarray):
            if len(items.shape) != 1:
                raise IndexError("Invalid indices: Must be a 1d array.")
            
        if (not isinstance(items, slice) and not isinstance(items, tuple) and
            not isinstance(items, list) and not isinstance(items, np.ndarray)):
            raise TypeError("Invalid indices: Must be slice, tuple, list, or numpy array.")
        
        out_object = deepcopy(self)
        if isinstance(items, tuple):
            out_object.expression_matrix = out_object.expression_matrix[items[0],:][:,items[1]]
            out_object.cell_types = out_object.cell_types[items[0]]
            out_object.sample_index = out_object.sample_index[items[0]]
            out_object.channels = out_object.channels[items[1]]
            if out_object.lineage_channels is not None:
                out_object.lineage_channels = out_object.lineage_channels[np.isin(out_object.lineage_channels, out_object.channels)]
            return out_object
        
        out_object.expression_matrix = out_object.expression_matrix[items]
        out_object.cell_types = out_object.cell_types[items]
        out_object.sample_index = out_object.sample_index[items]
        return out_object
    
         
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

        :return: An array of lineage channels or ``None``.
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
        self._lineage_channels_indices = np.where(np.isin(self.lineage_channels, self.channels))
        
        
    @property
    def reductions(self) -> Optional[dr.Reductions]:
        """Getter for reductions.

        :return: A ``Reductions`` object or ``None``.
        :rtype: CytofDR.dr.Reductions, optional
        """
        return self._reductions
    
    
    @reductions.setter
    def reductions(self, reductions: Optional[dr.Reductions]):
        """Set reductions.

        :param reductions: A ``Reductions`` object from the ``CytofDR`` package.
        :type reductions: int
        :raises TypeError: The object is not a ``Reductions`` object.
        """
        if not isinstance(reductions, dr.Reductions) and reductions is not None:
            raise TypeError("'reductions' has to of type 'CytofDR.dr.Reductions' or None")
        self._reductions = reductions


class DataLoader():
    """The class with utility functions to load datasets.

    This class offers one public utility function to load datasets, ``load_dataset``,
    which loads and preprocesses existing benchmark datasets. All other methods are
    private methods. Instantiation is not necessary.
    """

    # Package data directory and Path
    _data_dir = pkg_resources.resource_filename("PyCytoData", "data/")
    _data_path: Dict[str, str] = {"levine13": _data_dir + "levine13/",
                                  "levine32": _data_dir + "levine32/",
                                  "samusik": _data_dir + "samusik/"}


    @classmethod    
    def load_dataset(cls, dataset: str, sample: Optional[ArrayLike]=None, force_download: bool = False, preprocess: bool=False) -> PyCytoData:
        """Load benchmark datasets.

        This methods downloads and load benchmark datasets. The dataset is downloaded only once, which is then
        cached for future use. Currently, we support three datasets:
        
        - ``levine13``
        - ``levine32``
        - ``samusik``
        
        This method also supports specifying a specific sample instead of loading the entire dataset. Below is a list
        of samples available:
        
        - ``levine13``: ``0`` (There is only one sample in this case)
        - ``levine32``: ``AML08`` and ``AML09``.
        - ``samusik``: ``01``, ``02``, ..., ``09``, ``10``

        :param dataset: The name of the dataset.
        :type dataset: str
        :param sample: The specific sample to load from the dataset, defaults to None.
        :type sample: ArrayLike, optional
        :param force_download: Whether to download dataset regardless of previous cache, defaults to False
        :type force_download: bool
        :param preprocess: Whether to automatically perform all the necessary preocessing, defaults to false. 
            In the case of the existing three datasets, preprocessing includes just arcsinh transformation
            with cofactor of 5.
        :type preprocess: bool, optional
        :return: The loaded dataset.
        :rtype: PyCytoData
        """
        
        dataset = dataset.lower()
        if dataset not in ["levine13", "levine32", "samusik"]:
            raise ValueError("Unsupported dataset: Have to be 'levine13', 'levine32', or 'samusik'.")
                
        if not os.path.exists(cls._data_path[dataset]):
            cls._download_data(dataset = dataset, force_download = force_download)   
            
        data: PyCytoData = cls._preprocess(dataset)
        
        if sample is not None and not isinstance(sample, np.ndarray):
            sample = np.array(sample).flatten()
            data.subset(sample = sample)
                    
        if preprocess:
            data.preprocess(arcsinh=True, verbose=False)
            
        return data
            

    @classmethod
    def _download_data(cls,
                      dataset: str,
                      force_download: bool=False) -> int:
        """Method to download datasets."""
        urls: Dict[str, str] = {"levine13": "https://github.com/kevin931/PyCytoData/releases/download/datasets.rev.1/levine13.zip",
                                "levine32": "https://github.com/kevin931/PyCytoData/releases/download/datasets.rev.1/levine32.zip",
                                "samusik": "https://github.com/kevin931/PyCytoData/releases/download/datasets.rev.1/samusik.zip"}

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
        
        # Download
        contents = urlopen(urls[dataset])
        contents = contents.read()
        zip_file = ZipFile(BytesIO(contents))
        zip_file.extractall(cls._data_path[dataset])
        
        return 0
    
    
    @classmethod
    def _preprocess(cls, dataset: str) -> PyCytoData:
        """Preprocess the Samusik dataset."""
        fcs: str
        metadata: str
        data: PyCytoData
        
        if dataset == "levine13":
            fcs = cls._data_path[dataset] + "Levine_13dim_notransform.fcs"
            metadata = cls._data_path[dataset] + "Levine_13dim_cell_types.txt"
            data = cls._preprocess_levine13(fcs, metadata)
        elif dataset == "levine32":
            fcs = cls._data_path[dataset] + "Levine_32dim_notransform.fcs"
            metadata = cls._data_path[dataset] + "Levine_32dim_cell_types.txt"
            data = cls._preprocess_levine32(fcs, metadata)
        elif dataset == "samusik":
            fcs = cls._data_path[dataset] + "Samusik_all_notransform.fcs"
            metadata = cls._data_path[dataset] + "Samusik_cell_types.txt"
            data = cls._preprocess_samusik(fcs, metadata)
        else:
            raise ValueError("Unsupported dataset: Have to be 'levine13', 'levine32', or 'samusik'.")
            
        return data
    
    
    @classmethod
    def _preprocess_levine13(cls, fcs: str, metadata: str) -> PyCytoData:
        """Preprocess the Levine13 dataset."""
        df: pd.DataFrame
        _, df = fcsparser.parse(fcs, reformat_meta=True)
        labels: np.ndarray = np.loadtxt(metadata, delimiter="\t", dtype=str)
        
        labels = df["label"].apply(lambda x: "unassigned" if np.isnan(x) else labels[int(x),1]).to_numpy()
        df = df.drop(["label"], axis = 1)
        
        data: PyCytoData  = PyCytoData(expression_matrix=df.to_numpy(),
                                       cell_types=labels,
                                       channels=df.columns.to_numpy())
        data.lineage_channels = data.channels
        return data
    
    
    @classmethod
    def _preprocess_levine32(cls, fcs: str, metadata: str) -> PyCytoData:
        """Preprocess the Levine32 dataset and use the old formatting."""
        df: pd.DataFrame
        _, df = fcsparser.parse(fcs, reformat_meta=True)
        labels: np.ndarray = np.loadtxt(metadata, delimiter="\t", dtype=str)
        
        labels = df["label"].apply(lambda x: "unassigned" if np.isnan(x) else labels[int(x),1]).to_numpy()
        samples: np.ndarray = df["individual"].apply(lambda x: "AML08" if x == 1 else "AML09").to_numpy() #type: ignore
        df = df.drop(["individual", "label"], axis = 1)
        channels: np.ndarray = np.array(['Time', 'Cell_length', 'DNA1(Ir191)Di', 'DNA2(Ir193)Di', 'CD45RA(La139)Di',
                                         'CD133(Pr141)Di', 'CD19(Nd142)Di', 'CD22(Nd143)Di', 'CD11b(Nd144)Di',
                                         'CD4(Nd145)Di', 'CD8(Nd146)Di', 'CD34(Nd148)Di', 'Flt3(Nd150)Di', 'CD20(Sm147)Di',
                                         'CXCR4(Sm149)Di', 'CD235ab(Sm152)Di', 'CD45(Sm154)Di', 'CD123(Eu151)Di',
                                         'CD321(Eu153)Di', 'CD14(Gd156)Di', 'CD33(Gd158)Di', 'CD47(Gd160)Di', 'CD11c(Tb159)Di',
                                         'CD7(Dy162)Di', 'CD15(Dy164)Di', 'CD16(Ho165)Di', 'CD44(Er166)Di', 'CD38(Er167)Di',
                                         'CD13(Er168)Di', 'CD3(Er170)Di', 'CD61(Tm169)Di', 'CD117(Yb171)Di', 'CD49d(Yb172)Di',
                                         'HLA-DR(Yb174)Di', 'CD64(Yb176)Di', 'CD41(Lu175)Di', 'Viability(Pt195)Di', 'file_number', 'event_number'])
        
        data: PyCytoData  = PyCytoData(expression_matrix=df.to_numpy(),
                                       cell_types=labels,
                                       sample_index=samples,
                                       channels=channels)
        data.lineage_channels = data.channels[4:36]
        return data
    
    
    @classmethod
    def _preprocess_samusik(cls, fcs: str, metadata: str) -> PyCytoData:
        
        df: pd.DataFrame
        _, df = fcsparser.parse(fcs, reformat_meta=True)
        labels: np.ndarray = np.loadtxt(metadata, delimiter="\t", dtype=str)
        
        labels = df["label"].apply(lambda x: "unassigned" if np.isnan(x) else labels[int(x)-1]).to_numpy()
        samples: np.ndarray = df["sample"].apply(lambda x: str(int(x)) if x == 10 else f"0{str(int(x))}").to_numpy() #type: ignore
        df = df.drop(["sample", "event", "label"], axis = 1)
        
        data: PyCytoData  = PyCytoData(expression_matrix=df.to_numpy(),
                                       cell_types=labels,
                                       sample_index=samples,
                                       channels=df.columns.to_numpy())
        data.lineage_channels = data.channels[8:47]
        return data


class FileIO():
    """A utility class to handle common IO workflows for CyTOF data.

    This class includes a few utility static methods to load and save
    CyTOF data. Currently, it includes the following methods:
    
    - load_delim
    - load_expression
    - save_2d_list_to_csv
    - save_np_array
    
    Most of the methods are wrappers, but we offer a few advantages, such
    as returning ``PyCytoData`` data and saving ``numpy`` array along
    with channel names. For detailed documentations, read the docstring
    for each method.
    """
    
    @staticmethod
    def load_delim(files: Union[List[str], str],
                   skiprows: int=0,
                   drop_columns: Optional[Union[int, List[int]]]=None,
                   delim: str="\t",
                   dtype: type = float,
                   return_sample_indices: bool=False
                   ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        
        """Load deliminated file(s) as a numpy array.

        This method loads a deliminited file and returns a numpy array. The file
        has to be a standard text file. It is essentially a wrapper for the
        ``np.loadtxt`` function, but we offer the functionality of loading a list
        of files all at once, which are automatically concatenated.
    
        :param files: The path (or a list of paths) to the files to be loaded.
        :type files: Union[List[str], str]
        :param skiprows: The number of rows to skip, default to 0.
        :type skiprows: int, optional
        :param drop_colums: The columns indices for those that need to be dropped, defaults to None.
        :type drop_columns: Union[int, List[int]], optional.
        :param delim: The delimiter to use, defaults to ``\\t``
        :type delim: str, optional.
        :param dtype: The data type for the arrays, defaults to ``float``.
        :type dtype: type, optional

        :raises TypeError: The ``files`` is neither a string nor a list of strings.
        :return: An array or an array along with the sample indices.
        :rtype: Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]
        """
        
        if not isinstance(files, str) and not isinstance(files, list):
            raise TypeError("'files' has to be str or a list of str as paths.")
        elif not isinstance(files, list):
            files = [files]
        
        f: np.ndarray = np.array([])
        indices: np.ndarray = np.array([])
        for i, file in enumerate(files):                            
            # Load Data
            temp_f: np.ndarray = np.loadtxt(fname=file, dtype=dtype, skiprows=skiprows, delimiter=delim)
            if i==0:
                f = temp_f
                indices = np.repeat(i, temp_f.shape[0])
            else:
                f = np.vstack((f, temp_f))
                indices = np.hstack((indices, np.repeat(i, temp_f.shape[0])))
                
        if drop_columns is not None:
            f = np.delete(f, drop_columns, axis=1)
                
        if return_sample_indices:
            return f, indices
        else:
            return f
    
    
    @staticmethod
    def load_expression(files: Union[List[str], str],
                        col_names: bool=True,
                        drop_columns: Optional[Union[int, List[int]]]=None,
                        delim: str="\t",
                        dtype = float
                        ) -> PyCytoData:
        
        """Load a deliminited text file as a PyCytoData object.

        This method loads deliminited file(s) and returns a PyCytoData object. The file
        has to be a standard text file containing the expression matrix. Rows are cells
        and columns are channels. If ``col_names`` is ``True``, the first row of
        the file will be treated as channel names. If multiple file paths are present,
        they will be automatically concatenated into one object, but the sample
        indices will be recorded.
        
        :param files: The path (or a list of paths) to the files to be loaded.
        :type files: Union[List[str], str]
        :param col_names: Whether the first row is channel names, default to False.
        :type col_names: bool, optional
        :param drop_columns: The columns indices for those that need to be dropped, defaults to None.
        :type drop_columns: Union[int, List[int]], optional.
        :param delim: The delimiter to use, defaults to ``\\t``
        :type delim: str, optional.
        :param dtype: The data type for the arrays, defaults to ``float``.
        :type dtype: type, optional

        :raises TypeError: The ``files`` is neither a string nor a list of strings.
        :raises ValueError: The expression matrices' channels are mismatched or misaligned.
        :return: A PyCytoData object.
        :rtype: PyCytoData
        """
        
        if not isinstance(files, str) and not isinstance(files, list):
            raise TypeError("'files' has to be str or a list of str as paths.")
        elif not isinstance(files, list):
            files = [files]
        
        skiprows: int = 1 if col_names else 0
        exprs: np.ndarray
        indices: np.ndarray
        colnames: np.ndarray
        
        exprs, indices = FileIO.load_delim(files, skiprows, drop_columns, delim, dtype, return_sample_indices=True)
        data: PyCytoData = PyCytoData(expression_matrix=exprs, sample_index=indices)
        
        if col_names:
            colnames = np.loadtxt(fname=files[0], dtype ="str", max_rows=1, delimiter=delim)
            if drop_columns is not None:
                colnames = np.delete(colnames, drop_columns)
                
            # Check whether channels are aligned properly
            if len(files) > 1:
                f: int
                for f in range(1, len(files)):
                    temp_colnames: np.ndarray = np.loadtxt(fname=files[f], dtype ="str", max_rows=1, delimiter=delim)
                    if drop_columns is not None:
                        temp_colnames = np.delete(temp_colnames, drop_columns)
                    if not np.all(temp_colnames == colnames):
                        msg: str = f"The channels of expression matrices the first and {f+1}-th are not the same. "
                        msg += "Please ensure expression matrices' channels are in the same order with the same channels."
                        raise ValueError(msg)
                    
        else:
            colnames = np.full(data.n_channels, None)

        data.channels = colnames
        return data
    
    
    @staticmethod
    def save_2d_list_to_csv(data: List[List[Any]], path: str, overwrite: bool = False):
        """Save a nested list to a CSV file.

        :param data: The nested list to be written to disk
        :type data: List[List[Any]]
        :param path: Path to save the CSV file
        :type path: str
        
        .. note:: 
        
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
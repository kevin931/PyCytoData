from PyCytoData import FileIO, PyCytoData, DataLoader, exceptions
import numpy as np
import pandas as pd

import pytest
import os
import sys
from io import StringIO
import shutil
import csv
import _csv

from typing import List, Any, Dict, Tuple, Literal, Union, Optional
from numpy.typing import ArrayLike

OPT_PCK: Dict[str, bool] = {"CytofDR": True}

try:
    from CytofDR import dr
except ImportError:
    OPT_PCK["CytofDR"] = False


class TestCytoData():
    
    @classmethod
    def setup_class(cls):
        expression_matrix: np.ndarray = np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])
        channels: np.ndarray = np.array(["feature1", "feature2", "feature3"])
        cell_types: np.ndarray = np.array(["TypeA", "TypeB"])
        sample_index: np.ndarray = np.array(["SampleA", "SampleA"])
        lineage_channels: np.ndarray = np.array(["feature1", "feature2"])
        cls.dataset = PyCytoData(expression_matrix=expression_matrix, cell_types=cell_types, channels=channels,
                                 sample_index=sample_index, lineage_channels=lineage_channels)
        
        
    @pytest.mark.parametrize("attr,expected",
            [("expression_matrix", np.ndarray),
             ("channels", np.ndarray),
             ("cell_types", np.ndarray),
             ("sample_index", np.ndarray),
             ("n_channels", int),
             ("n_cells", int),
             ("n_samples", int),
             ("n_cell_types", int),
             ("lineage_channels", np.ndarray)])
    def test_attributes_type(self, attr: str, expected: Any):
        assert isinstance(getattr(self.dataset, attr), expected)
        
        
    @pytest.mark.parametrize("attr,expected",
            [("n_channels", 3),
             ("n_cells", 2),
             ("n_samples", 1),
             ("n_cell_types", 2)])
    def test_metadata_attributes_values(self, attr: str, expected: Any):
        assert getattr(self.dataset, attr) == expected
        
        
    def test_preprocess(self):
        expression_matrix: np.ndarray = np.abs(np.random.rand(500, 11))
        channels: np.ndarray = np.array(["feature1", "feature2", "Bead1", "Bead2", "Center", "Offset", "Residual", "Dead", "DNA1", "DNA2", "Time"])
        lineage_channels: np.ndarray = np.array(["feature1", "feature2"])
        cell_types: np.ndarray = np.repeat(["A", "B"], 250)
        sample_index: np.ndarray = np.repeat(["C", "D"], 250)
        dataset = PyCytoData(expression_matrix=expression_matrix, channels=channels,lineage_channels=lineage_channels,
                             sample_index=sample_index, cell_types=cell_types)
        dataset.preprocess(gate_debris_removal=True,
                           gate_center_offset_residual=True,
                           gate_live_cells=True,
                           gate_intact_cells=True,
                           bead_normalization=True,
                           arcsinh=True)
        assert dataset.expression_matrix.shape[0] <= 1000
        assert dataset.expression_matrix.shape[0] == dataset.n_cells
        assert dataset.sample_index.shape[0] == dataset.n_cells
        assert dataset.cell_types.shape[0] == dataset.n_cells
        
    
    @pytest.mark.parametrize("auto_channels", [True, False])
    def test_preprocess_manual_channels(self, auto_channels: bool):
        expression_matrix: np.ndarray = np.abs(np.random.rand(500, 11))
        channels: np.ndarray = np.array(["feature1", "feature2", "1bead", "2bead", "cen", "off", "res", "live", "dna1", "dna2", "clock"])
        lineage_channels: np.ndarray = np.array(["feature1", "feature2"])
        cell_types: np.ndarray = np.repeat(["A", "B"], 250)
        sample_index: np.ndarray = np.repeat(["C", "D"], 250)
        dataset = PyCytoData(expression_matrix=expression_matrix, channels=channels,lineage_channels=lineage_channels,
                             sample_index=sample_index, cell_types=cell_types)
        dataset.preprocess(gate_debris_removal=True,
                           gate_center_offset_residual=True,
                           gate_live_cells=True,
                           gate_intact_cells=True,
                           bead_normalization=True,
                           arcsinh=True,
                           cor_channels=["cen", "off", "res"],
                           DNA_channels=["dna1", "dna2"],
                           dead_channel=['live'],
                           bead_channels=["1bead", "2bead"],
                           time_channel=["clock"],
                           auto_channels=auto_channels)
        assert dataset.expression_matrix.shape[0] <= 1000
        assert dataset.expression_matrix.shape[0] == dataset.n_cells
        assert dataset.sample_index.shape[0] == dataset.n_cells
        assert dataset.cell_types.shape[0] == dataset.n_cells
    
    
    @pytest.mark.parametrize("channel_names,error_channel",
                             [(["feature1", "feature2", "Bead1", "Bead2", "cen", "Offset", "Residual", "Dead", "DNA1", "DNA2", "Time"], "cor_channels"),
                              (["feature1", "feature2", "1bead", "2bead", "Center", "Offset", "Residual", "Dead", "DNA1", "DNA2", "Time"], "bead_channels"),
                              (["feature1", "feature2", "Bead1", "Bead2", "Center", "Offset", "Residual", "live", "DNA1", "DNA2", "Time"], "dead_channel"),
                              (["feature1", "feature2", "Bead1", "Bead2", "Center", "Offset", "Residual", "Dead", "rna1", "rna2", "Time"], "DNA_channels"),
                              (["feature1", "feature2", "Bead1", "Bead2", "Center", "Offset", "Residual", "Dead", "DNA1", "DNA2", "clock"], "time_channel")])
    def test_preprocess_autochannel_error(self, channel_names: List[str], error_channel: str):
        expression_matrix: np.ndarray = np.abs(np.random.rand(500, 11))
        channels: np.ndarray = np.array(channel_names)
        lineage_channels: np.ndarray = np.array(["feature1", "feature2"])
        dataset = PyCytoData(expression_matrix=expression_matrix, channels=channels,lineage_channels=lineage_channels)
        try:
            dataset.preprocess(gate_debris_removal=True,
                               gate_center_offset_residual=True,
                               gate_live_cells=True,
                               gate_intact_cells=True,
                               bead_normalization=True)
        except exceptions.AutoChannelError as e:
            assert f"Auto channel detection failed for the following channels: {error_channel}." in str(e)
        else:
            assert False
    
    
    @pytest.mark.parametrize("attr,expected",
        [("expression_matrix", np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])),
         ("channels", np.array(["feature1", "feature2", "feature3"])),
         ("cell_types", np.array(["TypeA", "TypeB"])),
         ("sample_index", np.array(["SampleA", "SampleA"])),
         ("lineage_channels", np.array(["feature1", "feature2"]))])
    def test_array_attributes_values(self, attr: str, expected: Any):
        assert np.all(getattr(self.dataset, attr) == expected)
        
        
    def test_cytof_dr(self):
        if OPT_PCK["CytofDR"]:
            data: PyCytoData = PyCytoData(np.random.rand(100, 10))
            data.run_dr_methods(["PCA", "ICA"])
            assert isinstance(data.reductions, dr.Reductions)
            assert data.reductions.reductions["PCA"].shape == (100, 2)
            assert data.reductions.original_data is not None
            
    
    def test_cytof_dr_cell_types(self):
        if OPT_PCK["CytofDR"]:
            data: PyCytoData = PyCytoData(np.random.rand(100, 10),
                                          cell_types=np.repeat(["a", "b"], 50))
            data.run_dr_methods(["PCA", "ICA"])
            assert isinstance(data.reductions, dr.Reductions)
            assert data.reductions.reductions["PCA"].shape == (100, 2)
            assert data.reductions.original_cell_types is not None
            
            
    def test_reductions_setter(self):
        if OPT_PCK["CytofDR"]:
            self.dataset.reductions = dr.Reductions()
            assert isinstance(self.dataset.reductions, dr.Reductions)
            self.dataset.reductions = None
            assert self.dataset.reductions is None
            
            wrong_reductions: Any = []
            try:
                self.dataset.reductions = wrong_reductions
            except TypeError as e:
                assert "'reductions' has to of type 'CytofDR.dr.Reductions' or None" in str(e)
            else: 
                assert False            
        
        
    def test_add_sample(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        cell_types: np.ndarray = np.array(["TypeC", "TypeD"])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB"])
        self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index, cell_types=cell_types)
        assert self.dataset.expression_matrix.shape == (4,3)
        assert self.dataset.cell_types.shape[0] == 4
        assert self.dataset.sample_index.shape[0] == 4
        assert "SampleB" in self.dataset.sample_index
        assert "TypeC" in self.dataset.cell_types
        assert "TypeD" in self.dataset.cell_types
        
        
    def test_add_sample_sample_index_error(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        cell_types: np.ndarray = np.array(["TypeC", "TypeD"])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index, cell_types=cell_types)
        except exceptions.DimensionMismatchError as e:
            message = str(e)
            assert "The `sample_index` attribute has to be of length 2." in message
    
    
    def test_add_sample_cell_types_error(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        cell_types: np.ndarray = np.array(["TypeC", "TypeD", "TypeD"])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index, cell_types=cell_types)
        except exceptions.DimensionMismatchError as e:
            message = str(e)
            assert "The `cell_types` attribute has to be of length 2." in message
            
            
    def test_add_sample_expression_matrix_dimension_error(self):
        new_sample: np.ndarray = np.array([1.2, 2.3, 3.4])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index)
        except exceptions.ExpressionMatrixDimensionError as e:
            assert "The shape (3,) is unsupported. Please reshape it and ensure that the number of channels match." in str(e)
        else:
            raise
            
            
    def test_add_sample_expression_feature_number_error(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3], [3.4, 5.5]])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index)
        except exceptions.ExpressionMatrixDimensionError as e:
            assert "The shape (2, 2) is unsupported. Please reshape it and ensure that the number of channels match." in str(e)
        else:
            raise
            
            
    def test_add_sample_no_cell_types(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        sample_index: np.ndarray = np.array(["SampleC", "SampleC"])
        self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index)
        assert None in self.dataset.cell_types
        
        
    def test_setters(self):
        expression_matrix = np.array([[1.2, 2.3, 3.4, 4.5], [4.5, 5.6, 6.7, 7.8]])
        self.dataset.expression_matrix = expression_matrix
        self.dataset.channels = np.array(["feature1", "feature2", "feature3", "feature4"])
        self.dataset.sample_index = np.array(["SampleA", "SampleB"])
        self.dataset.cell_types = np.array(["TypeA", "TypeB"])
        self.dataset.lineage_channels = np.array(["feature3", "feature4"])
        assert self.dataset.n_cells == 2
        assert self.dataset.n_channels == 4
        assert self.dataset.n_samples == 2
        assert self.dataset.n_cell_types == 2
        assert np.all(np.isin(np.array(["feature1", "feature2", "feature3", "feature4"]), self.dataset.channels))
        assert np.all(np.isin(np.array(["feature3", "feature4"]), self.dataset.lineage_channels))
        
        
    def test_setters_n_values(self):
        self.dataset.n_cells = 2
        self.dataset.n_channels = 4
        self.dataset.n_samples = 2
        self.dataset.n_cell_types = 2
        assert self.dataset.n_cells == 2
        assert self.dataset.n_channels == 4
        assert self.dataset.n_samples == 2
        assert self.dataset.n_cell_types == 2
    
    
    @pytest.mark.parametrize("attr",
        ["n_channels","n_cells","n_samples", "n_cell_types"])
    def test_setter_type_error(self, attr: str):
        try:
            setattr(self.dataset, attr, "test")
        except TypeError as e:
            assert f"'{attr}' has to be 'int'" in str(e)
        else:
            assert False
        
        
    @pytest.mark.parametrize("attr,input_array,length",
        [("sample_index", np.array(["a", "a", "a"]), 2),
         ("cell_types", np.array(["a", "a", "a"]), 2),
         ("channels", np.array(["a", "a", "a"]), 2)])
    def test_setter_dimension_mismatch_error(self, attr: str, input_array: np.ndarray, length: int):
        try:
            setattr(self.dataset, attr, input_array)
        except exceptions.DimensionMismatchError as e:
            assert f"The `{attr}` attribute has to be of length {length}." in str(e)
        else:
            assert False
            
            
    def test_setter_expression_dimension_erro(self):
        try:
            self.dataset.expression_matrix = np.array([1, 2, 3])
        except exceptions.ExpressionMatrixDimensionError as e:
            assert "The shape (3,) is unsupported. Please reshape it and ensure that the number of channels match." in str(e)
        else:
            assert False
            
            
    def test_lineage_channel_setter_error(self):
        try:
            self.dataset.lineage_channels = np.array(["feature5"])
        except ValueError as e:
            assert "Some lineage channels are not listed in channel names." in str(e)
        else:
            raise
            
            
    def test_constructor_expression_matrix_dimmenion_error(self):
        try:
            PyCytoData(np.array([1,2,3]))
        except exceptions.ExpressionMatrixDimensionError as e:
            assert "The shape (3,) is unsupported. Please reshape it and ensure that the number of channels match." in str(e)
        else:
            assert False
            
    
    @pytest.mark.parametrize("channels,sample_index,cell_types,n,var",
                             [(np.array(["1", "2"]), None, None, 3, "channels"),
                              (None, np.array(["1", "1", "1"]), None, 2, "sample_index"),
                              (None, None, np.array(["1", "1", "1"]), 2, "cell_types")])
    def test_constructor_dimension_mismatch_error(self, channels: np.ndarray, sample_index: np.ndarray, cell_types: np.ndarray, n: int, var: str):
        expression: np.ndarray = np.array([[1,2,3],[2,1,3]])
        try:
            PyCytoData(expression, channels=channels, cell_types=cell_types, sample_index=sample_index)
        except exceptions.DimensionMismatchError as e:
            assert f"The `{var}` attribute has to be of length {n}." in str(e)
        else:
            assert False
            
            
    def test_constructor_lineage_channel_error(self):
        expression: np.ndarray = np.array([[1,2,3],[2,1,3]])
        channels: np.ndarray = np.array(["1", "2", "3"])
        lineage_channels: np.ndarray = np.array(["1", "2", "4"])
        try:
            PyCytoData(expression, channels=channels, lineage_channels=lineage_channels)
        except ValueError as e:
            assert "Some lineage channels are not listed in channel names." in str(e)
        else:
            assert False
            
            
    def test_constructor_channel_name_value_error(self):
        expression: np.ndarray = np.array([[1,2,3],[2,1,3]])
        channels: np.ndarray = np.array(["1", "1", "3"])
        try:
            PyCytoData(expression, channels=channels)
        except ValueError as e:
            assert "Channel names not unique: This can result in ambiguities." in str(e)
        else:
            assert False
            
            
    def test_subset_channels(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        lineage_channels: np.ndarray = np.array(["Channel0", "Channel1", "Channel2"])
        exprs = PyCytoData(exprs_matrix, lineage_channels=lineage_channels)
        exprs.subset(channels=["Channel1", "Channel2"])
        assert exprs.n_cells == 100
        assert exprs.n_channels == 2
        assert exprs.lineage_channels is not None
        assert exprs.lineage_channels.shape[0] == 2
        assert np.all(np.equal(exprs._lineage_channels_indices, np.array([0, 1])))
        assert not np.isin("Channel0", exprs.lineage_channels)
        
        
    def test_subset_channels_lineage_indices(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        channels: np.ndarray = np.arange(10).astype(str)
        exprs = PyCytoData(exprs_matrix, channels=channels)
        exprs.subset(channels=["1", "2"])
        assert exprs.n_cells == 100
        assert exprs.n_channels == 2
        assert exprs.lineage_channels is None
        assert exprs._lineage_channels_indices.shape[0] == 2
        assert np.all(np.equal(exprs._lineage_channels_indices, np.array([0, 1])))
        assert not np.isin("0", exprs.lineage_channels)
            
            
    def test_subset_cell_types(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB"], 50)

        exprs = PyCytoData(exprs_matrix, cell_types=cell_types)
        exprs.subset(cell_types="TypeA")
        assert exprs.n_cell_types == 1
        assert exprs.n_cells == 50
            
    def test_subset_sample(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        sample_index: np.ndarray = np.repeat(["SampleA", "SampleB"], 50)

        exprs = PyCytoData(exprs_matrix, sample_index=sample_index)
        exprs.subset(sample="SampleA")
        assert exprs.n_cell_types == 1
        assert exprs.n_cells == 50
        
        
    def test_subset_all(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        sample_index: np.ndarray = np.repeat(["SampleA", "SampleB"], 50)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 25)

        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        exprs.subset(channels=["Channel1", "Channel2"], sample="SampleA", cell_types="TypeA")
        assert exprs.n_cell_types == 1
        assert exprs.n_cells == 25
        assert exprs.n_channels == 2
    
    
    def test_subset_not_in(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 25)

        exprs = PyCytoData(exprs_matrix, cell_types=cell_types)
        exprs.subset(cell_types="TypeA", not_in=True)
        assert exprs.n_cell_types == 3
        assert exprs.n_cells == 75
        
        
    def test_subset_not_in_place(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 25)
        lineage_channels: np.ndarray = np.array(["Channel0", "Channel1", "Channel2"])

        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, lineage_channels=lineage_channels)
        new_exprs = exprs.subset(cell_types="TypeA", channels ="Channel0", not_in=True, in_place=False)
        assert isinstance(new_exprs, PyCytoData)
        assert new_exprs is not exprs
        assert new_exprs.n_cell_types == 3
        assert new_exprs.n_cells == 75
        assert new_exprs.n_channels == 9
        assert not np.isin("Channel0", new_exprs.lineage_channels)
        assert not np.isin("Channel0", new_exprs.channels)
        
    
    def test_subset_value_error_filter_all(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 25)

        exprs = PyCytoData(exprs_matrix, cell_types=cell_types)
        try:
            exprs.subset(cell_types=["TypeA", "TypeB", "TypeC", "TypeD"], not_in=True)
        except ValueError as e:
            assert "Filtering out all cells with nothing in the expression matrix. This is unsupported." in str(e)
        else:
            assert False
            
            
    def test_subset_type_error(self):
        exprs_matrix: np.ndarray = np.random.rand(100, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 25)

        exprs = PyCytoData(exprs_matrix, cell_types=cell_types)
        try:
            exprs.subset()
        except TypeError as e:
            assert "'channels', 'sample', and 'cell_types' cannot all be None." in str(e)
        else:
            assert False
        
        
    @pytest.mark.parametrize("features,n_channels",
                             [("feature1", 1), (["feature1", "feature2"], 2)]) 
    def test_get_channel_expressions(self, features: ArrayLike, n_channels: int):
        expression, channels = self.dataset.get_channel_expressions(features)
        assert isinstance(expression, np.ndarray)
        assert isinstance(channels, np.ndarray)
        assert expression.shape == (2, n_channels)
        assert channels.shape[0] == n_channels
        
        
    def test_get_channel_expressions_value_error(self):
        try:
            self.dataset.get_channel_expressions("channel100")
        except ValueError as e:
            assert "Some channels are not listed in channel names." in str(e)
        else:
            assert False
            
            
    def test_len(self):
        assert len(self.dataset) == self.dataset.n_cells
        
        
    def test_str(self):
        obj_str: str = f"A 'PyCytoData' object with {self.dataset.n_cells} cells, {self.dataset.n_channels} channels, {self.dataset.n_cell_types} cell types, and {self.dataset.n_samples} samples at"
        assert obj_str in str(self.dataset)
        
        
    def test_add(self):
        exprs1: PyCytoData = PyCytoData(np.random.rand(20, 10),
                                        sample_index=np.repeat("Sample1", 20),
                                        cell_types=np.repeat(["Type1", "Type2"], 10))
        exprs2: PyCytoData = PyCytoData(np.random.rand(30, 10),
                                        sample_index=np.repeat("Sample2", 30),
                                        cell_types=np.repeat(["Type3", "Type4"], 15))
        exprs3: PyCytoData = exprs1 + exprs2
        assert isinstance(exprs3, PyCytoData)
        assert exprs3.n_cells == 50
        assert exprs3.n_cell_types == 4
        assert exprs3.n_samples == 2
        
        
    def test_add_type_error(self):
        exprs1: PyCytoData = PyCytoData(np.random.rand(20, 10),
                                        sample_index=np.repeat("Sample1", 20),
                                        cell_types=np.repeat(["Type1", "Type2"], 10))
        try:
            exprs1 + 1 #type: ignore
        except TypeError as e:
            assert "The right hand side has to be a 'PyCytoData' object." in str(e)
        else:
            assert False
        
    
    def test_iadd(self):
        exprs1: PyCytoData = PyCytoData(np.random.rand(20, 10),
                                        sample_index=np.repeat("Sample1", 20),
                                        cell_types=np.repeat(["Type1", "Type2"], 10))
        exprs2: PyCytoData = PyCytoData(np.random.rand(30, 10),
                                        sample_index=np.repeat("Sample2", 30),
                                        cell_types=np.repeat(["Type3", "Type4"], 15))
        exprs1 += exprs2
        assert isinstance(exprs1, PyCytoData)
        assert exprs1.n_cells == 50
        assert exprs1.n_cell_types == 4
        assert exprs1.n_samples == 2
        
        
    def test_iadd_type_error(self):
        exprs1: PyCytoData = PyCytoData(np.random.rand(20, 10),
                                        sample_index=np.repeat("Sample1", 20),
                                        cell_types=np.repeat(["Type1", "Type2"], 10))
        try:
            exprs1 += 1 #type: ignore
        except TypeError as e:
            assert "The right hand side has to be a 'PyCytoData' object." in str(e)
        else:
            assert False
            
    
    @pytest.mark.parametrize("indices,n_cells,n_cell_types,n_samples",
                             [([0,1,2], 3,1,1),
                              (np.array([0,1,2]),3,1,1),
                              ([10], 1,1,1)]
                            )
    def test_getitem_1d(self, indices: Union[List[int], np.ndarray], n_cells: int, n_cell_types: int, n_samples: int):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        exprs = exprs[indices]
        
        assert exprs.n_channels == 10
        assert exprs.n_cells == n_cells
        assert exprs.n_cell_types == n_cell_types
        assert exprs.n_samples == n_samples
        
    
    @pytest.mark.parametrize("slice1,slice2,n_cells,n_cell_types,n_samples",
                             [(None,10,10,2,2),
                              (10,None,10,2,2),
                              (None,None,20,4,4),
                              (5,10,5,1,1)]
                            )   
    def test_getitem_1d_slice(self, slice1: Optional[int], slice2: Optional[int], n_cells: int, n_cell_types: int, n_samples: int):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        exprs = exprs[slice1:slice2]
        
        assert exprs.n_channels == 10
        assert exprs.n_cells == n_cells
        assert exprs.n_cell_types == n_cell_types
        assert exprs.n_samples == n_samples
        
    
    @pytest.mark.parametrize("index1,index2,n_cells,n_cell_types,n_samples,n_channels",
                             [([0,1,2],[0,1], 3,1,1,2),
                              (np.array([0,1,2]),np.array([0,1,2]),3,1,1,3),
                              ([10], [0], 1,1,1,1)]
                            )
    def test_getitem_2d(self,
                        index1: Union[slice, List[int], np.ndarray],
                        index2: Union[slice, List[int], np.ndarray],
                        n_cells: int,
                        n_cell_types: int,
                        n_samples: int,
                        n_channels: int):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        exprs = exprs[index1, index2]
        
        assert exprs.n_channels == n_channels
        assert exprs.n_cells == n_cells
        assert exprs.n_cell_types == n_cell_types
        assert exprs.n_samples == n_samples
    
    
    def test_getitem_2d_slice(self):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        exprs = exprs[:5,3:5]
        
        assert exprs.n_channels == 2
        assert exprs.n_cells == 5
        assert exprs.n_cell_types == 1
        assert exprs.n_samples == 1
        
        
    def test_getitem_2d_lineage_channels(self):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        lineage_channels: np.ndarray = np.array(["Channel1", "Channel2"])
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index, lineage_channels=lineage_channels)
        exprs = exprs[:5,:2]
        
        assert exprs.n_channels == 2
        assert exprs.n_cells == 5
        assert exprs.n_cell_types == 1
        assert exprs.n_samples == 1
        assert exprs.lineage_channels is not None
        assert not np.isin("Channel2", exprs.lineage_channels)
        assert exprs.lineage_channels.shape[0] == 1
        
        
    def test_getitem_1d_type_error(self):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        try:
            exprs[1] #type: ignore
        except TypeError as e:
            assert "Invalid indices: Must be slice, tuple, list, or numpy array." in str(e)
        else:
            assert False
            
    
    @pytest.mark.parametrize("index1,index2", 
                             [(1, [0,1]),
                              ([0,1], 1)])
    def test_getitem_2d_type_error(self, index1: Any, index2: Any):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        try:
            exprs[index1, index2]
        except TypeError as e:
            assert "Invalid indices: Must be slice, tuple, list, or numpy array." in str(e)
        else:
            assert False
            
            
    def test_getitem_2d_tuple_length(self):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        exprs = PyCytoData(exprs_matrix)
        try:
            exprs[:,:,:] #type: ignore
        except IndexError as e:
            assert "Invalid indices: Must be 1 or 2 indices." in str(e)
        else:
            assert False
            
    
    @pytest.mark.parametrize("index", 
                             [np.array(0),
                              np.array([[0,1],[0,1]])
                             ])
    def test_getitem_array_shape_error(self, index: np.ndarray):
        exprs_matrix: np.ndarray = np.random.rand(20, 10)
        cell_types: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        sample_index: np.ndarray = np.repeat(["TypeA", "TypeB", "TypeC", "TypeD"], 5)
        exprs = PyCytoData(exprs_matrix, cell_types=cell_types, sample_index=sample_index)
        try:
            exprs[index]
        except IndexError as e:
            assert "Invalid indices: Must be a 1d array." in str(e)
        else:
            assert False
            
            
    @classmethod
    def teardown_class(cls):
        pass
        

class TestDataLoader():
    
    @classmethod
    def setup_class(cls):
        os.makedirs("./tmp_pytest/data/levine13")
        os.makedirs("./tmp_pytest/data/levine32")
        os.makedirs("./tmp_pytest/data/samusik")
        
        for dataset in ["levine13", "levine32", "samusik"]:
            path: str
            if dataset == "levine13":
                path: str = "./tmp_pytest/data/"+ dataset + "/" + "Levine_13dim" + "_cell_types.txt"
                cell_types: np.ndarray = np.array([[1, "TypeA"], [2, "TypeB"]])
                FileIO.save_np_array(cell_types, path, dtype = "%s",
                                     col_names=np.array(["label", "population"]))
            elif dataset == "levine32":
                path: str = "./tmp_pytest/data/"+ dataset + "/" + "Levine_32dim" + "_cell_types.txt"
                cell_types: np.ndarray = np.array([[1, "TypeA"], [2, "TypeB"]])
                FileIO.save_np_array(cell_types, path, dtype = "%s",
                                     col_names=np.array(["label", "population"]))
            else:
                path: str = "./tmp_pytest/data/"+ dataset + "/" + "Samusik" + "_cell_types.txt"
                cell_types: np.ndarray = np.array(["TypeA", "TypeB"])
                FileIO.save_np_array(cell_types, path, dtype = "%s")
                
        
    def test_load_dataset_value_error(self):        
        try:
            DataLoader.load_dataset(dataset="levine14")
        except ValueError as e:
            assert "Unsupported dataset: Have to be 'levine13', 'levine32', or 'samusik'." in str(e)
        else:
            assert False
            
            
    def test_preprocess_value_eror(self):
        try:
            DataLoader._preprocess(dataset="levine14")
        except ValueError as e:
            assert "Unsupported dataset: Have to be 'levine13', 'levine32', or 'samusik'." in str(e)
        else:
            assert False
        
    
    @pytest.mark.parametrize("dataset",
                             ["levine13",
                              "levine32","samusik"]
                             )   
    def test_preprocess(self, mocker, dataset: Literal["levine13", "levine32", "samusik"]):
        preprocess_mock = mocker.MagicMock()
        mocker.patch("PyCytoData.DataLoader._preprocess_"+dataset, preprocess_mock)
        DataLoader._preprocess(dataset=dataset)
        
        assert preprocess_mock.called
        
        
    def test_preprocess_levine13(self, mocker):
        df = pd.DataFrame({"label": [1, 2, np.nan],
                           "CD1": [1.5, 1.5, 1.5],
                           "CD2": [2.5, 3.5, 4.5]})
        df = ("", df)
        mocker.patch("PyCytoData.data.fcsparser.parse", return_value = df)
        data: PyCytoData = DataLoader._preprocess_levine13(fcs = "",
                                                           metadata = "./tmp_pytest/data/levine13/Levine_13dim_cell_types.txt")
        assert isinstance(data, PyCytoData)
        assert np.all(np.equal(np.array(["TypeA", "TypeB", "unassigned"]), data.cell_types))
        
        
    def test_preprocess_levine32(self, mocker):
        df = pd.DataFrame({'Time': [1.5, 1.5, 1.5],
                           'Cell_length': [1.5, 1.5, 1.5],
                           'DNA1(Ir191)Di': [1.5, 1.5, 1.5],
                           'DNA2(Ir193)Di': [1.5, 1.5, 1.5],
                           'CD45RA(La139)Di': [1.5, 1.5, 1.5],
                           'CD133(Pr141)Di': [1.5, 1.5, 1.5],
                           'CD19(Nd142)Di': [1.5, 1.5, 1.5],
                           'CD22(Nd143)Di': [1.5, 1.5, 1.5],
                           'CD11b(Nd144)Di': [1.5, 1.5, 1.5],
                           'CD4(Nd145)Di': [1.5, 1.5, 1.5],
                           'CD8(Nd146)Di': [1.5, 1.5, 1.5],
                           'CD34(Nd148)Di': [1.5, 1.5, 1.5],
                           'Flt3(Nd150)Di': [1.5, 1.5, 1.5],
                           'CD20(Sm147)Di': [1.5, 1.5, 1.5],
                           'CXCR4(Sm149)Di': [1.5, 1.5, 1.5],
                           'CD235ab(Sm152)Di': [1.5, 1.5, 1.5],
                           'CD45(Sm154)Di': [1.5, 1.5, 1.5],
                           'CD123(Eu151)Di': [1.5, 1.5, 1.5],
                           'CD321(Eu153)Di': [1.5, 1.5, 1.5],
                           'CD14(Gd156)Di': [1.5, 1.5, 1.5],
                           'CD33(Gd158)Di': [1.5, 1.5, 1.5],
                           'CD47(Gd160)Di': [1.5, 1.5, 1.5],
                           'CD11c(Tb159)Di': [1.5, 1.5, 1.5],
                           'CD7(Dy162)Di': [1.5, 1.5, 1.5],
                           'CD15(Dy164)Di': [1.5, 1.5, 1.5],
                           'CD16(Ho165)Di': [1.5, 1.5, 1.5],
                           'CD44(Er166)Di': [1.5, 1.5, 1.5],
                           'CD38(Er167)Di': [1.5, 1.5, 1.5],
                           'CD13(Er168)Di': [1.5, 1.5, 1.5],
                           'CD3(Er170)Di': [1.5, 1.5, 1.5],
                           'CD61(Tm169)Di': [1.5, 1.5, 1.5],
                           'CD117(Yb171)Di': [1.5, 1.5, 1.5],
                           'CD49d(Yb172)Di': [1.5, 1.5, 1.5],
                           'HLA-DR(Yb174)Di': [1.5, 1.5, 1.5],
                           'CD64(Yb176)Di': [1.5, 1.5, 1.5],
                           'CD41(Lu175)Di': [1.5, 1.5, 1.5],
                           'Viability(Pt195)Di': [1.5, 1.5, 1.5],
                           'file_number': [1.5, 1.5, 1.5],
                           'event_number': [1.5, 1.5, 1.5],
                           "label": [1, 2, np.nan],
                           "individual": [1, 2, 1]})
        df = ("", df)
        mocker.patch("PyCytoData.data.fcsparser.parse", return_value = df)
        data: PyCytoData = DataLoader._preprocess_levine32(fcs = "",
                                                           metadata = "./tmp_pytest/data/levine32/Levine_32dim_cell_types.txt")
        channels: np.ndarray = np.array(['Time', 'Cell_length', 'DNA1(Ir191)Di', 'DNA2(Ir193)Di', 'CD45RA(La139)Di',
                                         'CD133(Pr141)Di', 'CD19(Nd142)Di', 'CD22(Nd143)Di', 'CD11b(Nd144)Di',
                                         'CD4(Nd145)Di', 'CD8(Nd146)Di', 'CD34(Nd148)Di', 'Flt3(Nd150)Di', 'CD20(Sm147)Di',
                                         'CXCR4(Sm149)Di', 'CD235ab(Sm152)Di', 'CD45(Sm154)Di', 'CD123(Eu151)Di',
                                         'CD321(Eu153)Di', 'CD14(Gd156)Di', 'CD33(Gd158)Di', 'CD47(Gd160)Di', 'CD11c(Tb159)Di',
                                         'CD7(Dy162)Di', 'CD15(Dy164)Di', 'CD16(Ho165)Di', 'CD44(Er166)Di', 'CD38(Er167)Di',
                                         'CD13(Er168)Di', 'CD3(Er170)Di', 'CD61(Tm169)Di', 'CD117(Yb171)Di', 'CD49d(Yb172)Di',
                                         'HLA-DR(Yb174)Di', 'CD64(Yb176)Di', 'CD41(Lu175)Di', 'Viability(Pt195)Di', 'file_number', 'event_number'])
        
        assert isinstance(data, PyCytoData)
        assert np.all(np.array(["TypeA", "TypeB", "unassigned"]) == data.cell_types)
        assert np.all(np.array(["AML08", "AML09", "AML08"]) == data.sample_index)
        assert np.all(channels == data.channels)
        
        
    def test_preprocess_samusik(self, mocker):
        df = pd.DataFrame({"label": [1, 2, np.nan],
                           "sample": [1.0, 2.0, 10],
                           "event": [1, 1, 1],
                           "CD1": [1.5, 1.5, 1.5],
                           "CD2": [2.5, 3.5, 4.5]})
        df = ("", df)
        mocker.patch("PyCytoData.data.fcsparser.parse", return_value = df)
        data: PyCytoData = DataLoader._preprocess_samusik(fcs = "",
                                                          metadata = "./tmp_pytest/data/samusik/Samusik_cell_types.txt")
        assert isinstance(data, PyCytoData)
        assert np.all(np.equal(np.array(["TypeA", "TypeB", "unassigned"]), data.cell_types))
        assert np.all(data.sample_index == np.array(["01", "02", "10"]))
        
        
    def test_load_dataset_sample_subset(self, mocker):
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/" + "levine13/",
                                                          "levine32": "./tmp_pytest/data/" + "levine32/",
                                                          "samusik": "./tmp_pytest/data/" + "samusik/"})
        
        df: PyCytoData = PyCytoData(expression_matrix=np.array([[0.2, 1], [1, 0.2], [2, 1.1]]),
                                    sample_index=["A", "B", "A"])
        mocker.patch("PyCytoData.DataLoader._preprocess", return_value = df)

        data: PyCytoData = DataLoader.load_dataset(dataset="levine32", sample="A")
        assert isinstance(data, PyCytoData)
        assert data.n_cells == 2
        
        
    def test_load_dataset_preprocess(self, mocker):
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/" + "levine13/",
                                                          "levine32": "./tmp_pytest/data/" + "levine32/",
                                                          "samusik": "./tmp_pytest/data/" + "samusik/"})
        
        df: PyCytoData = PyCytoData(expression_matrix=np.array([[0.2, 1], [1, 0.2], [2, 1.1]]))
        mocker.patch("PyCytoData.DataLoader._preprocess", return_value = df)

        data: PyCytoData = DataLoader.load_dataset(dataset="levine32", preprocess=True)
        assert isinstance(data, PyCytoData)
        assert data.expression_matrix[0,0] == np.arcsinh(0.2/5)
        
            
    @pytest.mark.parametrize("input_value",
                             ["n", "N"]
                             )     
    def test_download_assets_force_and_input_no(self, mocker, input_value):
        
        screen_stdout = sys.stdout
        string_stdout = StringIO()
        sys.stdout = string_stdout
        
        mocker.patch("builtins.input", return_value=input_value)
        out:int = DataLoader._download_data(dataset="levine13", force_download=False)
        
        output = string_stdout.getvalue()
        sys.stdout = screen_stdout
        
        assert "You have declined to download" in output    
        assert "Download in progress" not in output
        assert out == 1   
            
    
    @pytest.mark.parametrize("input_value,force_download",
                            [("y", False),
                             ("Y", False),
                             ("y", True)])   
    def test_download_assets_url_download_extract_mock(self, mocker, input_value: str, force_download: bool):
        contents_mock = mocker.MagicMock()
        contents_mock.read.return_value = b"This is a test."
        zip_mock = mocker.MagicMock()
        zip_mock.extractall.return_value = None
        
        mocker.patch("builtins.input", return_value=input_value)
        mocker.patch("PyCytoData.data.urlopen", return_value = contents_mock)
        mocker.patch("PyCytoData.data.ZipFile", return_value = zip_mock)
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        out: int = DataLoader._download_data(dataset="levine13", force_download=force_download)
        
        zip_mock.extractall.assert_called_once()
        contents_mock.read.assert_called_once()
        assert out == 0
        
        
    def test_load_dataset_download(self, mocker):
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/levine13/",
                                                          "levine32": "./tmp_pytest/data/levine32/",
                                                          "samusik": "./tmp_pytest/data/samusik/"})
        df: PyCytoData = PyCytoData(expression_matrix=np.array([[0.2, 1], [1, 0.2], [2, 1.1]]),
                                    sample_index=["A", "B", "A"])
        mocker.patch("PyCytoData.DataLoader._preprocess", return_value = df)
        download_mock = mocker.MagicMock()
        mocker.patch("PyCytoData.DataLoader._download_data", download_mock)
        
        shutil.rmtree("./tmp_pytest/data/levine32/")

        data: PyCytoData = DataLoader.load_dataset(dataset="levine32")
        download_mock.assert_called_once()
        
        
    def test_purge_dataset_cache(self, mocker):
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/levine13/",
                                                          "levine32": "./tmp_pytest/data/levine32/",
                                                          "samusik": "./tmp_pytest/data/samusik/"})
        DataLoader.purge_dataset_cache(dataset="levine13")
        assert not os.path.exists("./tmp_pytest/data/levine13")
        
        
    def test_purge_dataset_cache_all(self, mocker):
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/levine13/",
                                                          "levine32": "./tmp_pytest/data/levine32/",
                                                          "samusik": "./tmp_pytest/data/samusik/"})
        DataLoader.purge_dataset_cache_all()
        assert not os.path.exists("./tmp_pytest/data/samusik")
        
        
    def test_purge_dataset_cache_value_error(self):
        try:
            DataLoader.purge_dataset_cache(dataset="WrongDataset")
        except ValueError as e:
            assert str(e) == "Unsupported dataset: Have to be 'levine13', 'levine32', or 'samusik'."
        else:
            assert False
        
    
    @classmethod
    def teardown_class(cls):
        shutil.rmtree("./tmp_pytest/")
        

class TestFileIO():
    
    @classmethod
    def setup_class(cls):
        os.mkdir("./tmp_pytest/")
        
        with open("./tmp_pytest/file_read_test_csv.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter=",")
            tsv_writer.writerow(["col1", "col2", "col3"])
            tsv_writer.writerow([1, 2, 3])
            tsv_writer.writerow([4, 5, 6])
            
        with open("./tmp_pytest/file_read_test_csv2.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter=",")
            tsv_writer.writerow(["col1", "col2", "col3"])
            tsv_writer.writerow([7, 8, 9])
            tsv_writer.writerow([10, 11, 12])
            
        with open("./tmp_pytest/file_read_test_csv3.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter=",")
            tsv_writer.writerow(["col3", "col2", "col1"])
            tsv_writer.writerow([7, 8, 9])
            tsv_writer.writerow([10, 11, 12])
            
        with open("./tmp_pytest/file_read_test_tsv.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter="\t")
            tsv_writer.writerow([1.1, 2.2, 3.3])
            tsv_writer.writerow([4.4, 5.5, 6.6])
        
    
    @pytest.mark.parametrize("path,delim,col_names,dtype",
                [("./tmp_pytest/file_read_test_csv.txt", ",", True, int),
                 ("./tmp_pytest/file_read_test_tsv.txt", "\t", False, float)]
                )
    def test_load_expression_filetype(self, path: str, delim: str, col_names: bool, dtype):
        out_file: PyCytoData = FileIO.load_expression(path, col_names, delim=delim, dtype = dtype)
        assert isinstance(out_file, PyCytoData)
        assert out_file.n_cells == 2
        assert out_file.n_channels == 3
        
        if path == "./tmp_pytest/file_read_test_csv.txt":
            assert "col1" in out_file.channels
            assert out_file.expression_matrix.dtype == np.dtype("int64")
        else:
            assert None in out_file.channels
            assert out_file.expression_matrix.dtype == np.dtype("float64")
            
            
    def test_load_expression_concat(self):
        out_file: PyCytoData = FileIO.load_expression(["./tmp_pytest/file_read_test_csv.txt", "./tmp_pytest/file_read_test_csv2.txt"],
                                                 col_names=True, delim=",", dtype = int)
        assert isinstance(out_file, PyCytoData)
        assert out_file.n_cells == 4
        assert out_file.n_channels == 3
        assert np.all(np.isin([0, 1], out_file.sample_index))
        
        
    def test_laod_expression_mismatch_error(self):
        try:
            out_file: PyCytoData = FileIO.load_expression(["./tmp_pytest/file_read_test_csv2.txt", "./tmp_pytest/file_read_test_csv3.txt"],
                                                          col_names=True, delim=",", dtype = int)
        except ValueError as e:
            msg: str = "The channels of expression matrices the first and 2-th are not the same. "
            msg += "Please ensure expression matrices' channels are in the same order with the same channels."
            assert msg in str(e)
        else:
            assert False
            
            
    def test_laod_expression_mismatch_drop(self):
        out_file: PyCytoData = FileIO.load_expression(["./tmp_pytest/file_read_test_csv2.txt", "./tmp_pytest/file_read_test_csv3.txt"],
                                                      col_names=True, delim=",", dtype = int, drop_columns=[0,2])
        assert out_file.n_channels == 1
        assert out_file.n_cells == 4
    
    
    @pytest.mark.parametrize("drop_cols,expected_shape",
            [([0, 1], (2, 1)), (1, (2,2))]
            )      
    def test_load_expression_drop_col(self, drop_cols, expected_shape):
        path: str = "./tmp_pytest/file_read_test_csv.txt"
        out_file: PyCytoData = FileIO.load_expression(path, True, drop_columns=drop_cols, delim=",", dtype = int)
        assert out_file.expression_matrix.shape == expected_shape
        
        
    def test_load_delim_concat(self):
        out_file: Union[np.ndarray, Tuple[np.ndarray, np.ndarray]] = FileIO.load_delim(["./tmp_pytest/file_read_test_csv.txt", "./tmp_pytest/file_read_test_csv2.txt"],
                                                                                       skiprows=1, delim=",", dtype = int, return_sample_indices=False)
        assert isinstance(out_file, np.ndarray)
        assert out_file.shape == (4, 3)
            
    
    @pytest.mark.parametrize("drop_cols,expected_shape",
            [([0, 1], (2, 1)), (1, (2,2))]
            )      
    def test_load_delim_drop_col(self, drop_cols, expected_shape):
        path: str = "./tmp_pytest/file_read_test_csv.txt"
        out_file: Union[np.ndarray, Tuple[np.ndarray, np.ndarray]] = FileIO.load_delim(path,
                                                                                       1,
                                                                                       drop_columns=drop_cols,
                                                                                       delim=",",
                                                                                       dtype = int,
                                                                                       return_sample_indices=False)
        assert isinstance(out_file, np.ndarray)
        assert out_file.shape == expected_shape
        
    
    @pytest.mark.parametrize("method", ["load_delim", "load_expression"])
    def test_load_expression_path_type_error(self, method: str):
        try:
            path: Any = {"a": "./tmp_pytest/file_read_test_csv.txt"}
            getattr(FileIO, method)(path)
        except TypeError as e:
            assert "'files' has to be str or a list of str as paths." in str(e)
        else:
            assert False
        
        
    @pytest.mark.parametrize("save_list,path",
                    [([[1, 2, 3], [1, 2, 3]], "./tmp_pytest/2d_list_1.csv"),
                    ([[1.4, 2.2], ["a", "b"]], "./tmp_pytest/2d_list_2.csv"),
                    ([[1.4, 1], [1, "a"]], "./tmp_pytest/2d_list_3.csv")]
                    )   
    def test_save_2d_list_to_csv(self, save_list: List[List[Any]], path: str):
        FileIO.save_2d_list_to_csv(save_list, path)
        assert os.path.exists(path)


    def test_save_2d_list_to_csv_exception(self):
        path: str = "./tmp_pytest/2d_list_3.csv"
        save_list: List[List[Any]] = [[1.4, 1], [1, "a"]]
        try:
            FileIO.save_2d_list_to_csv(save_list, path)
        except FileExistsError:
            assert True
            

    def test_save_np_array(self):
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        path: str = "./tmp_pytest/nparr.txt"
        FileIO.save_np_array(arr, path)
        assert os.path.exists("./tmp_pytest/nparr.txt")
        
        
    def test_save_np_array_colnames(self): 
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        header: "np.ndarray" = np.array(["A", "B"])
        path: str = "./tmp_pytest/nparr_colnames.txt"
        FileIO.save_np_array(arr, path, col_names=header)
        assert os.path.exists("./tmp_pytest/nparr_colnames.txt")
        
    
    def test_save_np_array_exception(self):
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        path: str = "./tmp_pytest/nparr.txt"
        try:
            FileIO.save_np_array(arr, path)
        except FileExistsError:
            assert True
    
    
    @classmethod
    def teardown_class(cls):
        shutil.rmtree("./tmp_pytest/")
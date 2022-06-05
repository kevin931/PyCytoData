from PyCytoData import preprocess
import numpy as np

import pytest

class TestPreprocess():
    
    @classmethod
    def setup_class(cls):
        cls.data: np.ndarray = np.abs(np.random.rand(100, 11))
        cls.channels: np.ndarray = np.array(["feature1", "feature2", "Bead1", "Bead2", "Center", "Offset", "Residual", "Dead", "DNA1", "DNA2", "Time"])
        cls.lineage_channels: np.ndarray = np.array(["feature1", "feature2"])
        cls.bead_channels: np.ndarray = np.array(["Bead1", "Bead2"])
        cls.cor_channels: np.ndarray = np.array(["Center", "Offset", "Residual"])
        cls.dead_channel: np.ndarray = np.array(["Dead"])
        cls.DNA_channels: np.ndarray = np.array(["DNA1", "DNA2"])
        cls.time_channel: np.ndarray = np.array(["Time"])
        
    
    @pytest.mark.parametrize("cofactor", [5, 10])
    def test_arcsinh(self, cofactor: int):
        out: np.ndarray = preprocess.arcsinh(self.data, self.channels, self.lineage_channels, cofactor=cofactor)
        assert isinstance(out, np.ndarray)
        assert out.shape == (100, 11)
        
        
    def test_arcsinh_arraylike(self):
        out: np.ndarray = preprocess.arcsinh(self.data.tolist(), self.channels.tolist(), self.lineage_channels.tolist())
        assert isinstance(out, np.ndarray)
        assert out.shape == (100, 11)
        
        
    def test_arcsinh_transform_channel_none(self):
        out: np.ndarray = preprocess.arcsinh(self.data.tolist(), self.channels.tolist())
        assert isinstance(out, np.ndarray)
        assert out.shape == (100, 11)
        
        
    def test_arcsinh_channel_none(self):
        out: np.ndarray = preprocess.arcsinh(self.data.tolist())
        assert isinstance(out, np.ndarray)
        assert out.shape == (100, 11)
        
        
    def test_gate_debris_removal(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_debris_removal(self.data, self.channels, self.bead_channels)
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
    
    def test_gate_debris_removal_arraylike(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_debris_removal(self.data.tolist(), self.channels.tolist(), self.bead_channels.tolist())
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
    
    
    @pytest.mark.parametrize("cutoff", [1, 2])
    def test_gate_intact_cells(self, cutoff: float):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_intact_cells(self.data, self.channels, self.DNA_channels, cutoff_DNA_sd=cutoff)
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    def test_gate_intact_cells_arraylike(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_intact_cells(self.data.tolist(), self.channels.tolist(), self.DNA_channels.tolist())
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    @pytest.mark.parametrize("cutoff", [0.01, 0.03])
    def test_gate_live_cells(self, cutoff: float):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_live_cells(self.data, self.channels, self.dead_channel, cutoff_quantile=cutoff)
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    def test_gate_live_cells_arraylike(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_live_cells(self.data.tolist(), self.channels.tolist(), self.dead_channel.tolist())
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
    
    @pytest.mark.parametrize("cutoff", [0.01, 0.03])
    def test_gate_center_offset_residual(self, cutoff: float):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_center_offset_residual(self.data, self.channels, self.cor_channels, cutoff_quantile=cutoff)
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    def test_gate_center_offset_residual_arraylike(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.gate_center_offset_residual(self.data.tolist(), self.channels.tolist(), self.cor_channels.tolist())
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    def test_bead_normalization(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.bead_normalization(self.data, self.channels, self.bead_channels, self.time_channel, self.lineage_channels)
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    def test_bead_normalization_arraylike(self):
        out: np.ndarray
        indices: np.ndarray
        out, indices = preprocess.bead_normalization(self.data.tolist(), self.channels.tolist(), self.bead_channels.tolist(),
                                                     self.time_channel.tolist(), self.lineage_channels.tolist())
        assert isinstance(out, np.ndarray)
        assert isinstance(indices, np.ndarray)
        assert out.shape[0] == indices.shape[0]
        assert out.shape[0] <= 100
        
        
    def test_bead_normalization_time_channel_error(self):
        time_channel = np.array(["Time", "DNA1"])
        try:
            preprocess.bead_normalization(self.data.tolist(), self.channels.tolist(), self.bead_channels.tolist(),
                                          time_channel.tolist(), self.lineage_channels.tolist())
        except ValueError as e:
            assert "Only 1 'Time' channel allowed." in str(e)
        else:
            raise
        
        
    def test_gate_live_cell_dead_channel_error(self):
        dead_channel = np.array(["Dead", "DNA1"])
        try:
            preprocess.gate_live_cells(self.data, self.channels, dead_channel)
        except ValueError as e:
            assert "Only 1 'Dead' channel allowed." in str(e)
        else:
            raise
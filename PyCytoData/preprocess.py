import numpy as np

from numpy.typing import ArrayLike
from typing import Optional, List


def arcsinh(data: ArrayLike, channels: Optional[ArrayLike]=None, transform_channels: Optional[ArrayLike]=None, cofactor: int=5) -> np.ndarray:
    """Arcsinh transformation for CyTOF data.

    Arcsinh transformation is often the first step to preprocessing data. This function flexibly allows users
    to transform their data at a cofactor of their choice and to specify a transformation of their choice.

    :param data: The expression matrix array of two dimensions.
    :param channels: The channel names of the expression matrix in the order of the columns, defaults to None
    :param transform_channels: The channels to transformed as specify by name, defaults to None
    :param cofactor: The cofactor, defaults to 5
    
    :return: The arcsinh transformed expression matrix.
    """
    
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    if channels is None or (channels is not None and transform_channels is None):
        return np.arcsinh(data/cofactor)
    
    if not isinstance(channels, np.ndarray):
        channels = np.array(channels)
    if not isinstance(transform_channels, np.ndarray):
        transform_channels = np.array(transform_channels)
    
    channels = channels.flatten()
    transform_channels = transform_channels.flatten()
        
    data[:,np.isin(channels, transform_channels)] = np.arcsinh(data[:,np.isin(channels, transform_channels)]/cofactor)
    return data


def gate_debris_removal(data: ArrayLike, channels: ArrayLike, bead_channels: ArrayLike) -> np.ndarray:
    """Pre-gating to remove debris.

    This is a first step of the gating procedures to remove debris. Channels names and bead channel names are needed.

    :param data: The expression matrix array of two dimensions.
    :param channels: The channel names of the expression matrix in the order of the columns, defaults to None
    :param bead_channels: The bead channels as specify by name, defaults to None
    
    :return: The gated expression matrix.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    if not isinstance(channels, np.ndarray):
        channels = np.array(channels)
    if not isinstance(bead_channels, np.ndarray):
        bead_channels = np.array(bead_channels)
        
    channels = channels.flatten()
    bead_channels = bead_channels.flatten()
    
    bead_exp: np.ndarray = np.mean(data[:, np.isin(channels, bead_channels)], axis=1)
    cutoff_bead: float = np.mean(bead_exp) + np.std(bead_exp)*3
    
    return data[bead_exp < cutoff_bead,:]
    
    
def gate_intact_cells(data: ArrayLike, channels: ArrayLike, DNA_channels: ArrayLike, cutoff_DNA_sd: float=2) -> np.ndarray:
    """Gating for intact cells.

    This gating procedure gates for intact cells following the debris removal procedure.
    All channel names and DNA channel names needed.

    :param data: The expression matrix array of two dimensions.
    :param channels: The channel names of the expression matrix in the order of the columns, defaults to None
    :param DNA_channels: The DNA channels as specify by name, defaults to None
    :param cutoff_DNA_sd: The number of standard deviations away from the mean to use as a cutoff for DNA channels.
    
    :return: The gated expression matrix.
    """
    
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    if not isinstance(channels, np.ndarray):
        channels = np.array(channels)
    if not isinstance(DNA_channels, np.ndarray):
        DNA_channels = np.array(DNA_channels)
    
    channels = channels.flatten()
    DNA_channels = DNA_channels.flatten()

    cutoff_DNA1_se: float = cutoff_DNA_sd*np.std(data[:,np.isin(channels, DNA_channels[0])])/3
    cutoff_DNA2_se: float = cutoff_DNA_sd*np.std(data[:,np.isin(channels, DNA_channels[1])])/3
    
    cutoff_DNA1_upper: float = np.mean(data[:,np.isin(channels, DNA_channels[0])]) + cutoff_DNA1_se
    cutoff_DNA1_lower: float = np.mean(data[:,np.isin(channels, DNA_channels[0])]) - cutoff_DNA1_se
    cutoff_DNA2_upper: float = np.mean(data[:,np.isin(channels, DNA_channels[1])]) + cutoff_DNA2_se
    cutoff_DNA2_lower: float = np.mean(data[:,np.isin(channels, DNA_channels[1])]) - cutoff_DNA2_se
    
    return data[np.where((data[:,np.isin(channels, DNA_channels[0])] > cutoff_DNA1_lower).flatten() & 
                         (data[:,np.isin(channels, DNA_channels[0])] < cutoff_DNA1_upper).flatten() & 
                         (data[:,np.isin(channels, DNA_channels[1])] > cutoff_DNA2_lower).flatten() &
                         (data[:,np.isin(channels, DNA_channels[1])] < cutoff_DNA2_upper).flatten())[0],:]
    
    
def gate_live_cells(data: ArrayLike, channels: ArrayLike, dead_channel: ArrayLike, cutoff_quantile: float=0.03) -> np.ndarray:
    """Gating for live cells.

    This gating procedure gates for living cells following the gating procedure for intact cells.
    All channel names and 'Dead' channel names needed.

    :param data: The expression matrix array of two dimensions.
    :param channels: The channel names of the expression matrix in the order of the columns, defaults to None
    :param dead_channel: The dead channels as specify by name, defaults to None
    :param cutoff_quantile: The top quantile to be excluded.
    
    :return: The gated expression matrix.
    
    :raises ValueError: More than 1 "Dead" channel provided.
    """
    
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    if not isinstance(channels, np.ndarray):
        channels = np.array(channels)
    if not isinstance(dead_channel, np.ndarray):
        dead_channel = np.array(dead_channel)
        
    channels = channels.flatten()   
    dead_channel = dead_channel.flatten()
    
    if dead_channel.shape[0] > 1:
        raise ValueError("Only 1 'Dead' channel allowed.")
    
    cutoff_dead: float = np.quantile(data[:,np.isin(channels, dead_channel)], 1-cutoff_quantile)
    return data[(data[:, np.isin(channels, dead_channel)] < cutoff_dead).flatten(),:]
    
    
def gate_center_offset_residual(data: ArrayLike, channels: ArrayLike, cor_channels: ArrayLike, cutoff_quantile: float=0.03):
    """Gating for live cells.

    This gating procedure gates for living cells following the gating procedure for intact cells.
    All channel names and 'Dead' channel names needed.

    :param data: The expression matrix array of two dimensions.
    :param channels: The channel names of the expression matrix in the order of the columns, defaults to None
    :param cor_channels: The center, offset, and residual channels as specify by name, defaults to None
    :param cutoff_quantile: The top and bottom quantile to be excluded.
    
    :return: The gated expression matrix.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    if not isinstance(channels, np.ndarray):
        channels = np.array(channels)
    if not isinstance(cor_channels, np.ndarray):
        cor_channels = np.array(cor_channels)
        
    channels = channels.flatten()
    cor_channels = cor_channels.flatten()
    
    cutoff0_upper: float = np.quantile(data[:,np.isin(channels, cor_channels[0])], 1-cutoff_quantile)
    cutoff0_lower: float = np.quantile(data[:,np.isin(channels, cor_channels[0])], cutoff_quantile)
    cutoff1_upper: float = np.quantile(data[:,np.isin(channels, cor_channels[1])], 1-cutoff_quantile)
    cutoff1_lower: float = np.quantile(data[:,np.isin(channels, cor_channels[1])], cutoff_quantile)
    cutoff2_upper: float = np.quantile(data[:,np.isin(channels, cor_channels[2])], 1-cutoff_quantile)
    cutoff2_lower: float = np.quantile(data[:,np.isin(channels, cor_channels[2])], cutoff_quantile)
    return data[np.where((data[:,np.isin(channels, cor_channels[0])] > cutoff0_lower).flatten() & 
                         (data[:,np.isin(channels, cor_channels[0])] < cutoff0_upper).flatten() &
                         (data[:,np.isin(channels, cor_channels[1])] > cutoff1_lower).flatten() & 
                         (data[:,np.isin(channels, cor_channels[1])] < cutoff1_upper).flatten() &
                         (data[:,np.isin(channels, cor_channels[2])] > cutoff2_lower).flatten() & 
                         (data[:,np.isin(channels, cor_channels[2])] < cutoff2_upper).flatten())[0],:]
    
    
def bead_normalization(data: ArrayLike, channels: ArrayLike, bead_channels: ArrayLike, time_channel: ArrayLike, transform_channels: ArrayLike):
    """Bead normalization to correct time-shift throughout an experiment.

    Sometimes, CyTOF has a phenomenon known as time-dependent shift, meaning that the over expression
    becomes biased over time. To correct this, bead normaliztion is used.

    :param data: The expression matrix array of two dimensions.
    :param channels: The channel names of the expression matrix in the order of the columns
    :param bead_channels: The bead channels as specify by name
    :param time_channel: The time channels as specify by name
    :param transform_quantile: The transform channels to apply the normalization as specify by name
    
    :return: The normalized expression matrix.
    
    :raises ValueError: More than 1 "Time" channel provided.
    """
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    if not isinstance(channels, np.ndarray):
        channels = np.array(channels)
    if not isinstance(bead_channels, np.ndarray):
        bead_channels = np.array(bead_channels)
    if not isinstance(time_channel, np.ndarray):
        time_channel = np.array(time_channel)
    if not isinstance(transform_channels, np.ndarray):
        transform_channels = np.array(transform_channels)
        
    channels = channels.flatten()
    time_channel = time_channel.flatten()
    bead_channels = bead_channels.flatten()
    transform_channels = transform_channels.flatten()
    
    if time_channel.shape[0] > 1:
        raise ValueError("Only 1 'Time' channel allowed.")
    
    intervals: np.ndarray = np.round(data[:,np.isin(channels, time_channel)]/np.max(data[:,np.isin(channels, time_channel)])*200).flatten()
    unique_intervals: np.ndarray = np.unique(intervals)
    
    # Remove aberrant data
    i: int
    interval_mean: np.ndarray = np.empty((unique_intervals.shape[0], data.shape[1]), dtype=float)
    for i in range(0, unique_intervals.shape[0]):
        interval_mean[i] = np.mean(data[np.isin(intervals, unique_intervals[i]),:], axis=0)
    
    keep: np.ndarray = np.full((unique_intervals.shape[0],), 0)
    for c in bead_channels:
        keep += (interval_mean[:,np.isin(channels, c)] > np.quantile(interval_mean[:,np.isin(channels, c)], 0.05)).astype(int).flatten()
    keep = keep > (bead_channels.shape[0]/2)
    unique_intervals = unique_intervals[keep]
    data = data[np.isin(intervals, unique_intervals),:]
    
    # Normalization
    intervals: np.ndarray = np.round(data[:,np.isin(channels, time_channel)]/np.max(data[:,np.isin(channels, time_channel)])*200).flatten() # type: ignore
    unique_intervals: np.ndarray = np.unique(intervals)
    interval_mean: np.ndarray = np.empty((unique_intervals.shape[0], data.shape[1]), dtype=float) # type: ignore
    
    for c in bead_channels:
        data[:,np.isin(channels, c)] = data[:, np.isin(channels, c)]/np.mean(data[:, np.isin(channels, c)]) #type: ignore
        
    if bead_channels.shape[0] > 1:
        cormat: np.ndarray = np.corrcoef(interval_mean[:, np.isin(channels, bead_channels)], rowvar=False)
        np.fill_diagonal(cormat, -1)
        cormat[np.tril_indices_from(cormat)] = -1
        max_cor: int = np.argmax(cormat) #type: ignore
        indices: List[int] = list(np.unravel_index(max_cor, [cormat.shape[0], cormat.shape[1]])) #type: ignore
        bead_channels = bead_channels[indices]
        
        
    correct_f: np.ndarray = np.mean(interval_mean[:, np.isin(channels, bead_channels)], axis=1)
    correct_f_full: np.ndarray = np.empty(intervals.shape[0])
    for i in range(0, unique_intervals.shape[0]):
        correct_f_full[np.isin(intervals, unique_intervals[i])] = correct_f[i]
    correct_f_full = correct_f_full/np.mean(correct_f_full)
    keep = correct_f_full > 0
    
    for c in transform_channels:
        data[:, np.isin(channels, c)] = data[:, np.isin(channels, c)]/(correct_f_full.reshape(data.shape[0], 1)) #type: ignore
    return data[keep,:] #type: ignore

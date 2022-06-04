from typing import Tuple, List

class ExpressionMatrixDimensionError(Exception):
    """Error for Expression Matrix Dimension Error.

    This error is thrown when the given expression matrix is not or cannot
    be cast into a two dimensional array.

    :param shape: The actual shape of the expression matrix.
    :type shape: Tuple[int, ...]
    """
    
    def __init__(self, shape: Tuple[int, ...]):
        self.shape = shape
        super().__init__()
        
    def __str__(self):
        return f"The shape {self.shape} is unsupported. Please reshape it and ensure that the number of channels match."
    

class DimensionMismatchError(Exception):
    """Error for mismatched dimensions.

    This error is thrown when the given attribute does not match the length
    of the data in the class.

    :param n: The desired length of the attribute.
    :type n: int
    :param var: The name of the attribute.
    :type var: str
    """
    def __init__(self, n: int, var: str):
        self.n = n
        self.var = var
        super().__init__()
        
    def __str__(self):
        return f"The `{self.var}` attribute has to be of length {self.n}."
    
    
class AutoChannelError(Exception):
    """Auto Channel Failure.

    This error is thrown when the ``auto_channels`` option has failed due
    to regex unable to match necessary channels.

    :param channels: The necessary channels.
    :type channels: List[str]
    """
    def __init__(self, channel: List[str]):
        self.channel: str = ", ".join(channel)
        super().__init__()
        
    def __str__(self):
        return f"Auto channel detection failed for the following channels: {self.channel}."
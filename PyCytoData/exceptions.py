from typing import Tuple

class ExpressionMatrixDimensionError(Exception):
    
    def __init__(self, shape: Tuple[int, ...]):
        self.shape = shape
        super().__init__()
        
    def __str__(self):
        return f"The shape {self.shape} is unsupported. Please reshape it and ensure that the number of channels match."
    

class DimensionMismatchError(Exception):
    def __init__(self, n: int, var: str):
        self.n = n
        self.var = var
        super().__init__()
        
    def __str__(self):
        return f"The `{self.var}` attribute has to be of length {self.n}."
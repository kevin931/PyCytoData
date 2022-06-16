import pytest
from PyCytoData import __main__
from PyCytoData.__init__ import __VERSION__

import sys
from io import StringIO

def test_main():
    
    # Redirect StdOut
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    # Capture Output
    __main__.main()
    output = string_stdout.getvalue()
    # Set stdout back
    sys.stdout = screen_stdout
    #Check
    assert f"PyCytoData Version: v{__VERSION__}" in output
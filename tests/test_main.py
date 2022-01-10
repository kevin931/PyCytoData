import pytest
from PythonTemplate import __main__

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
    assert "Hello World" in output
"""
There are a lot of different versions of Clipper-Python going around at the
moment, so this submodule is a way of trying to adapt to whichever version
might be available in any particular Python environment.
"""

try:
    # Import CCP4 CP
    from clipper import *
    mode = 0
except ImportError:
    print('WARNING: failed to import Clipper-Python, some functions will be unavailable')
    mode = -1

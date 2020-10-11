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
    try:
        # If on Windows, import local CP
        import os
        if os.name == 'nt':
            import sys
            if sys.version_info[0] == 2:
                from py2.clipper import *
            elif sys.version_info[0] == 3:
                raise ImportError #from py3.clipper import *
            mode = 0
        # Else, import pip-installed CP
        else:
            from clipper_python._clipper import *
            mode = 1
    except ImportError:
        print('WARNING: failed to import Clipper-Python, some functions will be unavailable')
        mode = -1

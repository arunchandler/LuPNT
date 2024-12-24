import logging

try:
    from ._blender import *
except ImportError as e:
    print(e)
    try:
        from _blender import *
    except ImportError:
        pass

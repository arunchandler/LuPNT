try:
    from ._chebyshev import *
    from ._polynomial import *
    from ._sinusodial import *
    from ._angles_fitting import *
except ImportError:
    from _chebyshev import *
    from _polynomial import *
    from _sinusodial import *
    from _angles_fitting import *

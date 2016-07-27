try:
    from ._fastdtw import fastdtw
except ImportError:
    from .fastdtw import fastdtw
    # user has been warned on installation
    pass    
